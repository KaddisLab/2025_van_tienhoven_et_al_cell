#' Merge Multiple Seurat Objects
#'
#' Merges multiple Seurat objects into one, applying a specified project name to the merged object.
#' Each object is loaded and merged in sequence. Memory usage is reported after each merge.
#'
#' @param seurat_object A character vector of paths to Seurat object files to be merged.
#' @param project_name A character string specifying the project name for the merged Seurat object.
#' @return The path to the saved merged Seurat object file.
#' @importFrom SeuratObject JoinLayers
#' @importFrom glue glue
#' @importFrom qs qsave
#' @seealso \code{\link[SeuratObject:JoinLayers]{JoinLayers}} for joining layers,
#'   \code{\link[glue:glue]{glue}} for constructing paths,
#'   and \code{\link[qs:qsave]{qsave}} for saving objects efficiently.
#' @examples
#' \dontrun{
#' seurat_paths <- c("path/to/seurat1.rds", "path/to/seurat2.rds")
#' project_name <- "MyMergedProject"
#' merged_path <- seurat_merge(seurat_paths, project_name)
#' }
seurat_merge2 <- function(seurat_objects, project_name) {
    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    message("Loading objects to list...")
    # Parallel loading of Seurat objects
    seurat_list <- future.apply::future_lapply(seurat_objects, function(obj_path) {
        loaded_object <- load_seurat(obj_path)
        return(loaded_object)
    }, future.seed = 42)
    gc()
    # Extract sample IDs and set them as names of the list elements
    names(seurat_list) <- sapply(seurat_list, function(obj) {
        sample_id <- obj[["orig.ident"]][1, ]
        return(sample_id)
    })

    message("Loaded ", length(seurat_list), " objects")
    message("List size: ", format(object.size(seurat_list), units = "Gb"))
    message("Merging objects...")
    object <- merge(seurat_list[[1]], seurat_list[-1])
    rm(seurat_list)
    gc()
    message("Merged size: ", format(object.size(object), units = "Gb"))
    object <- SeuratObject::JoinLayers(object)
    message("Joined size: ", format(object.size(object), units = "Gb"))
    Seurat::Project(object) <- project_name
    object_path <- glue::glue("{analysis_cache}/seurat_merged/{project_name}.qs")
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    qs::qsave(object, file = object_path)
    return(object_path)
}
