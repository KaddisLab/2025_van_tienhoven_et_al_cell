#' Transfer cell labels from reference to query Seurat Object
#'
#' This function applies cell type annotations from a reference Seurat object
#' to a query Seurat object using SingleR. It supports optional normalization
#' of the reference object and uses differential expression methods to assist
#' with the annotation transfer.
#'
#' @param query_seurat_object a seurat object or path to the QS file for the query Seurat object.
#' @param ref_seurat_object a seurat object or path to the QS file for the reference Seurat object.
#' @param cell_type_col The name of the metadata column in the reference Seurat object
#'        that contains cell type annotations.
#' @param out_folder The name of the folder to save the output files.
#' @param de_method The differential expression method to be used by SingleR. Possible
#'        values are "wilcox", "bimod", "t", "negbinom", "poisson", or "LR".
#'        Default is "wilcox".
#'
#' @return The path to the CSV file containing the transferred cell type annotations
#'         for the query Seurat object.
#'
#' @details The function first performs normalization on the query Seurat object. It
#'          optionally performs normalization on the reference Seurat object as well.
#'          SingleR is then used to transfer cell type annotations from the reference
#'          to the query dataset. The function is designed to work in a high-performance
#'          computing environment, leveraging parallel processing capabilities.
#'
#' @examples
#' query_seurat_object <- "path/to/query_seurat_object.qs"
#' ref_seurat_object <- "path/to/reference_seurat_object.qs"
#' cell_type_col <- "cell_type"
#' 
#' # Transfer cell type annotation without normalizing the reference
#' result_path <- seurat_transfer_cell_type_annotation(
#'   query_seurat_object,
#'   ref_seurat_object
#'   cell_type_col,
#'   norm_ref = FALSE,
#'   de_method = "wilcox"
#' )
#' 
#' # Check the result
#' print(result_path)
#'
#' @importFrom Seurat NormalizeData GetAssayData
#' @importFrom qs qread
#' @importFrom BiocParallel register MulticoreParam bpparam
#' @importFrom SingleR SingleR
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom glue glue
#' @export
seurat_singleR_transfer_label <- function(query_seurat_object, ref_seurat_object, cell_type_col, de_method = "wilcox", out_folder = "cell_type_out") {
    query_seurat_object <- load_seurat(query_seurat_object)
    ref_seurat_object <- load_seurat(ref_seurat_object)

    sample_id <- query_seurat_object|> Seurat::Project()
    ref_id <- ref_seurat_object |> Seurat::Project()

    if (!isNormalized(query_seurat_object)) {
        query_seurat_object <- Seurat::NormalizeData(query_seurat_object)
    }

    if (!isNormalized(ref_seurat_object)) {
        ref_seurat_object <- Seurat::NormalizeData(ref_seurat_object)
    }

    # # # Setup parallel processing
    # BiocParallel::register(
    #     BiocParallel::MulticoreParam(workers = hprcc::slurm_allocation()$CPUs, progressbar = TRUE)
    # )
    ncpus <- if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
        as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    } else if (Sys.getenv("SLURM_CPUS_ON_NODE") != "") {
        as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
    } else {
        1 # default if not in SLURM environment
    }

    # Make it reproducible
    set.seed(42)

    # Run SingleR
    cell_labels <- SingleR::SingleR(
        test = Seurat::GetAssayData(query_seurat_object, assay = "RNA"),
        ref =  Seurat::GetAssayData(ref_seurat_object, assay = "RNA"),
        labels = ref_seurat_object@meta.data[[cell_type_col]],
        de.method = de_method,
        BPPARAM = BiocParallel::MulticoreParam(workers = ncpus)
    ) 
    # drop the scores: it's a matrix of scores for each cell type vs each cell type. We just want 
    # the label and the delta.next, which is the difference between the score for the next best cell type
    cell_labels$scores<-NULL
    
    # rename
    cell_type_table <- cell_labels |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "cell") |>
        tibble::as_tibble() |>
        dplyr::select(cell,
            cell_type = pruned.labels,
            cell_type_label = labels,
            delta_next = delta.next) |>
        dplyr::rename_with(~ paste0(ref_id, .), -cell)
    
    cell_type_table_path <- glue::glue("{analysis_cache}/{out_folder}/{sample_id}_{ref_id}_cell_type.csv")

    dir.create(dirname(cell_type_table_path), showWarnings = FALSE, recursive = TRUE)

    cell_type_table |> write.csv(cell_type_table_path, row.names = FALSE, quote = FALSE)

    return(cell_type_table_path)
}
