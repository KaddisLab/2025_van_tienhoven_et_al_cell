#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seurat_object(s)
#' @return
#' @author Denis O'Meally
#' @export
seurat_merge_and_sketch <- function(seurat_object, n_cells = 750) {

    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    if (length(seurat_object) > 1) {
        message("Loading objects...")
        object_list<-list()
        for (object in seurat_object) {
            object <- load_seurat(object)
            sample_id <- object[["orig.ident"]][1, ]
            object_list[[sample_id]] <- object
            message("Loaded ", sample_id)
        }
        message("Merging objects...")
        object <- scCustomize::Merge_Seurat_List(object_list)
    } else {
        message("Loading Seurat object", seurat_object)
        object <- load_seurat(seurat_object)
    }
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
