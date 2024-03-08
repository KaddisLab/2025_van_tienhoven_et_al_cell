#' Transfer Cell Type Annotation from Reference to Test Seurat Object
#'
#' This function applies cell type annotations from a reference Seurat object
#' to a test Seurat object using SingleR. It supports optional normalization
#' of the reference object and uses differential expression methods to assist
#' with the annotation transfer.
#'
#' @param test_seurat_object a seurat object or path to the QS file for the test Seurat object.
#' @param ref_seurat_object a seurat object or path to the QS file for the reference Seurat object.
#' @param cell_type_col The name of the metadata column in the reference Seurat object
#'        that contains cell type annotations.
#' @param de_method The differential expression method to be used by SingleR. Possible
#'        values are "wilcox", "bimod", "t", "negbinom", "poisson", or "LR".
#'        Default is "wilcox".
#'
#' @return The path to the CSV file containing the transferred cell type annotations
#'         for the test Seurat object.
#'
#' @details The function first performs initial quality control (using an assumed
#'          `initialQC` function) and normalization on the test Seurat object. It
#'          optionally performs normalization on the reference Seurat object as well.
#'          SingleR is then used to transfer cell type annotations from the reference
#'          to the test dataset. The function is designed to work in a high-performance
#'          computing environment, leveraging parallel processing capabilities.
#'
#' @examples
#' test_seurat_object <- "path/to/test_seurat_object.qs"
#' ref_seurat_object <- "path/to/reference_seurat_object.qs"
#' cell_type_col <- "cell_type"
#' 
#' # Transfer cell type annotation without normalizing the reference
#' result_path <- seurat_transfer_cell_type_annotation(
#'   test_seurat_object,
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
seurat_transfer_cell_type_annotation <- function(test_seurat_object, ref_seurat_object, cell_type_col, de_method = "wilcox") {

    if (!inherits(test_seurat_object, "Seurat")) {
        test_seurat_object <- qs::qread(test_seurat_object)
    }
    if (!inherits(ref_seurat_object, "Seurat")) {
        ref_seurat_object <- qs::qread(ref_seurat_object)
    }

    test_seurat_object <- test_seurat_object |> initialQC()
    
    ref_seurat_object <- ref_seurat_object |> initialQC() 

    if (!isNormalized(test_seurat_object)) {
        test_seurat_object <- Seurat::NormalizeData(test_seurat_object)
    }

    if (!isNormalized(ref_seurat_object)) {
        ref_seurat_object <- Seurat::NormalizeData(ref_seurat_object)
    }

    # Setup parallel processing
    BiocParallel::register(
        BiocParallel::MulticoreParam(workers = hprcc::slurm_allocation()$CPUs, progressbar = TRUE)
    )

    # Make it reproducible
    set.seed(42)

    # Run SingleR
    cell_labels <- SingleR::SingleR(
        test = Seurat::GetAssayData(test_seurat_object, assay = "RNA"),
        ref =  Seurat::GetAssayData(ref_seurat_object, assay = "RNA"),
        labels = ref_seurat_object@meta.data[[cell_type_col]],
        de.method = de_method,
        BPPARAM = BiocParallel::bpparam()
    ) 
    # drop the scores: it's a matrix of scores for each cell type vs each cell type. We just want 
    # the label and the delta.next, which is the difference between the score for the next best cell type
    cell_labels$scores<-NULL
    
    # rename
    cell_type_table <- cell_labels |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = "cell") |> 
        tibble::as_tibble()

    sample_id <- test_seurat_object@meta.data$orig.ident[1]
    ref_id <- ref_seurat_object@project.name
    
    cell_type_table_path <- glue::glue("{analysis_cache}/cell_type_out/{sample_id}_{ref_id}_cell_type.csv")

    dir.create(dirname(cell_type_table_path), showWarnings = FALSE, recursive = TRUE)

    cell_type_table |> write.csv(cell_type_table_path)

    return(cell_type_table_path)
}
