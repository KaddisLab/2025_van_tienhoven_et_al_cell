#' Detect Cell Doublets Using scDblFinder
#'
#' Performs doublet detection on a Seurat object using the scDblFinder algorithm.
#' The function handles data normalization if needed, supports both automated and
#' pre-existing clustering, and is compatible with CellPlex experiments.
#'
#' @param seurat_object A Seurat object or path to a Seurat object file
#' @param cluster_col Name of an existing cluster column in metadata. If NULL,
#'        clusters will be computed using scDblFinder::fastcluster.
#' @param cmo_column For CellPlex assays, the name of the column containing CMO/sample
#'        information. Set to NULL for non-CellPlex data.
#'
#' @return Path to a CSV file containing doublet detection results, with columns:
#'   - cell: Cell barcode
#'   - class: Classification (singlet/doublet)
#'   - score: Doublet score
#'   Additional scDblFinder metrics if available
#'
#' @details
#' The function:
#' 1. Normalizes data if not already normalized
#' 2. Converts to SingleCellExperiment
#' 3. Performs clustering (if needed)
#' 4. Runs scDblFinder algorithm
#' 5. Saves results to "{analysis_cache}/scDblFinder_out/{sample_id}_scDdblFinder.csv"
#'
#' @import Seurat SingleCellExperiment
#' @importFrom scDblFinder fastcluster scDblFinder
#' @importFrom scran computeSumFactors
#' @importFrom scuttle logNormCounts
#' @importFrom BiocParallel MulticoreParam
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select any_of contains rename_with
#' @importFrom stringr str_remove
#' @author Denis O'Meally
#' @export
seurat_scDblFinder <- function(seurat_object, cluster_col = NULL, cmo_column = NULL) {

    seurat_object <- load_seurat(seurat_object)

#    bp <- BiocParallel::MulticoreParam(hprcc::slurm_allocation()$CPUs, RNGseed = 42)
ncpus <- if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
    as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else if (Sys.getenv("SLURM_CPUS_ON_NODE") != "") {
    as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
} else {
    1 # default if not in SLURM environment
}
bp <- BiocParallel::MulticoreParam(workers = ncpus)

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
