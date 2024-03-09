#' Run scDblFinder for doublet detection on a Seurat object
#'
#' @title seurat_scDblFinder
#' @param seurat_object A Seurat object.
#' @param clusters A column name from the Seurat object's metadata to use for clustering information.
#'                 If NULL, clusters will be computed using scDblFinder::fastcluster.
#' @param sample_column The name of the column to use for sample identifiers in scDblFinder. 
#'                      If NULL, defaults to 'hash.ID'.
#' @return The modified Seurat object with updated metadata.
#' @author Denis O'Meally
#' @export
seurat_scDblFinder <- function(seurat_object, cluster_col = NULL, sample_column = NULL) {

    # Initialisation for parallel processing ------------------------------------
    init_multisession()
    resources <- slurm_allocation()
    bp <- BiocParallel::MulticoreParam(resources$CPUs, RNGseed = 42)

    Seurat::DefaultAssay(seurat_object) <- "RNA"
    sce <- Seurat::as.SingleCellExperiment(seurat_object)

    if (is.null(cluster_col)) {
        sce$clusters <- scDblFinder::fastcluster(sce, BPPARAM = bp)
    } else {
        sce$clusters <- seurat_object[[cluster_col]]
    }

    if (!is.null(sample_column) && !(sample_column %in% colnames(SingleCellExperiment::colData(sce)))) {
        stop("The specified sample_column does not exist in the Seurat object's metadata.")
    }
    
    sce <- try(
        {
            sce <- sce |> scran::computeSumFactors(cluster = sce$clusters, BPPARAM = bp) 
            sce <- sce |> scuttle::logNormCounts(BPPARAM = bp) 
            sce <- sce |> scDblFinder::scDblFinder(clusters = "clusters", samples = sample_column, verbose = TRUE, BPPARAM = bp)
        }
        #silent = TRUE
    )

    if (class(sce) == "try-error") {
        message("scDblFinder failed")
        seurat_object@meta.data$scDblFinder.class <- "scDblFinder failed"
    } else {
        seurat_object@meta.data <- SingleCellExperiment::colData(sce) |> as.data.frame()
    }

    return(seurat_object)
}
