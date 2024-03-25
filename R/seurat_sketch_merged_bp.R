#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param merged_seurat_path
#' @return
#' @author Denis O'Meally
#' @export
seurat_sketch_merged_bp <- function(seurat_object, n_cells = 1000) {

    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    object <- load_seurat(seurat_object)
    object <- Seurat::NormalizeData(object)
    object <- Seurat::FindVariableFeatures(object, verbose = FALSE)

    #  Sample representative cells from each dataset
    object <- Seurat::SketchData(object = object, ncells = n_cells, method = "LeverageScore", sketched.assay = "sketch")
    object_path <- glue::glue("{analysis_cache}/data/merged_seurat_bp_sketch.rds")

    # Save seurat object
    saveRDS(object, file = object_path)
}
