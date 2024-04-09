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
seurat_merge <- function(seurat_object, project_name) {
    format_gb <- function(bytes) {
        gb <- bytes / 2^30
        formatted_gb <- sprintf("%.2f Gb", gb)
        return(formatted_gb)
    }
    if (length(seurat_object) > 1) {
        message("Loading objects...")
        object <- load_seurat(seurat_object[1])
        sample_id <- object[["orig.ident"]][1, ]
        message("Loaded ", sample_id)
        for (i in 2:length(seurat_object)) {
            new_object <- load_seurat(seurat_object[i])
            sample_id <- new_object[["orig.ident"]][1, ]
            object <- merge(object, new_object)
            rm(new_object)
            gc()
            message("Merged ", sample_id, " | ", i, " of ", length(seurat_object), " | ", format_gb(lobstr::mem_used()))
        }
    } else {
        message("Loading Seurat object...")
        object <- load_seurat(seurat_object)
    }
    object <- SeuratObject::JoinLayers(object)
    Seurat::Project(object) <- project_name
    object_path <- glue::glue("{analysis_cache}/seurat_merged/{project_name}.qs")
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    qs::qsave(object, file = object_path)
    return(object_path)
}
