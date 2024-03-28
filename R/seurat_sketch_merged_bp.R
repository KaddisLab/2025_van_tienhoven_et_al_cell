#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param merged_seurat_path
#' @return
#' @author Denis O'Meally
#' @export
seurat_sketch_merged_bp <- function(seurat_object, n_cells = 750) {

    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()
    message("Loading Seurat object", seurat_object)
    object <- load_seurat(seurat_object)
    message("Joining layers")
    object <- SeuratObject::JoinLayers(object)
    #TODO Use scTransform here to regress out unwanted sources of variation?
    message("Normalizing data")
    object <- Seurat::NormalizeData(object)
    message("Splitting by sample")
    object[["RNA"]] <- split(object[["RNA"]], f = object$orig.ident)
    message("Finding variable features")
    object <- Seurat::FindVariableFeatures(object, verbose = FALSE)
    message("Sketching data")
    #  Sample representative cells from each dataset
    object <- Seurat::SketchData(
        object = object,
        ncells = n_cells,
        method = "LeverageScore",
        sketched.assay = "sketch")

    object_path <- glue::glue("{analysis_cache}/data/merged_seurat_bp_sketch_{n_cells}cells.rds")
    message("Saving object to ", object_path)
    # Save seurat object
    saveRDS(object, file = object_path)
    message("Done")
    return(object_path)
}
