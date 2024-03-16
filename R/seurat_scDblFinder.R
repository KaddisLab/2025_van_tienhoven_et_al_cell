#' Run scDblFinder for doublet detection on a Seurat object
#'
#' @title seurat_scDblFinder
#' @param seurat_object A Seurat object.
#' @param clusters A column name from the Seurat object's metadata to use for clustering information.
#'                 If NULL, clusters will be computed using scDblFinder::fastcluster.
#' @param cmo_column For CellPlex assays, the name of the column that indicates CMO or sample, or NULL (default).
#' 
#' @return The path to a modified Seurat object with updated metadata containing the scDblFinder results.
#' @author Denis O'Meally
#' @export
seurat_scDblFinder <- function(seurat_object, cluster_col = NULL, cmo_column = NULL) {

    seurat_object <- load_seurat(seurat_object)

    bp <- BiocParallel::MulticoreParam(hprcc::slurm_allocation()$CPUs, RNGseed = 42)

    Seurat::DefaultAssay(seurat_object) <- "RNA"

    if (!isNormalized(seurat_object)) {
        seurat_object <- Seurat::NormalizeData(seurat_object)
    }

    sce <- Seurat::as.SingleCellExperiment(seurat_object)

    if (is.null(cluster_col)) {
        sce$clusters <- scDblFinder::fastcluster(sce, BPPARAM = bp)
    } else {
        sce$clusters <- seurat_object[[cluster_col]]
    }

    if (!is.null(cmo_column) && !(cmo_column %in% colnames(SingleCellExperiment::colData(sce)))) {
        stop("The specified cmo_column does not exist in the Seurat object's metadata.")
    }
    
    sce <- try({
            sce <- sce |> scran::computeSumFactors(cluster = sce$clusters, BPPARAM = bp) 
            sce <- sce |> scuttle::logNormCounts(BPPARAM = bp) 
            sce <- sce |> scDblFinder::scDblFinder(clusters = "clusters", samples = cmo_column, verbose = TRUE, BPPARAM = bp)
        })

    if (class(sce) == "try-error") {
        message("scDblFinder failed")
        seurat_object@meta.data$scDblFinder.class <- "scDblFinder failed"
    } else {
        seurat_object@meta.data <- SingleCellExperiment::colData(sce) |> as.data.frame()
    }

    seurat_metadata <- seurat_object@meta.data %>% 
        tibble::rownames_to_column("cell") %>% 
        dplyr::select("cell", dplyr::any_of(dplyr::contains("scDblFinder"))) %>% 
        dplyr::rename_with(~ stringr::str_remove(., "scDblFinder\\."))
    
    sample_id <- seurat_object[[]]$orig.ident[1]
    
    seurat_metadata_path <- glue::glue("{analysis_cache}/scDblFinder_out/{sample_id}_scDdblFinder.csv")
    dir.create(dirname(seurat_metadata_path), showWarnings = FALSE, recursive = TRUE)

    seurat_metadata |> write.csv(seurat_metadata_path, row.names = FALSE, quote = FALSE)

    return(seurat_metadata_path)
}
