#' Sketch Seurat Object
#'
#' Subsamples a Seurat object to a specified number of cells, normalizes the data,
#' finds variable features, and performs sketching based on leverage scores.
#' The sketched data is then saved to a specified path.
#'
#' @param seurat_object The path to the Seurat object file to be loaded.
#' @param n_cells The number of cells to subsample to; defaults to 750.
#' @return The path to the saved sketched Seurat object.
#' @export
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     seurat_path <- "path/to/seurat_object.rds"
#'     sketched_path <- seurat_sketch(seurat_path, n_cells = 500)
#' }
#' }
#' @importFrom Seurat NormalizeData FindVariableFeatures
#' @importFrom qs qsave
#' @importFrom glue glue
seurat_sketch <- function(seurat_object, n_cells = 750) {
    object <- load_seurat(seurat_object)
    sample_id <- Seurat::Project(object)
    object <- Seurat::NormalizeData(object)
    object <- Seurat::FindVariableFeatures(object, verbose = FALSE)
    #  Sample representative cells from each dataset
    object <- Seurat::SketchData(
        object = object,
        ncells = n_cells,
        method = "LeverageScore",
        sketched.assay = "sketch",
        seed = 42
    )
    
    object_path <- glue::glue("{analysis_cache}/seurat_sketch_out/{sample_id}_{n_cells}_cells.qs")
    save_results(object, object_path)
    return(object_path)
}

