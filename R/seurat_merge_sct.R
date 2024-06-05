#' Merge RNA assay from Multiple Seurat Objects
#'
#' Merges RNA assays from multiple Seurat objects into one, applying a specified project name to the merged object.
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
seurat_merge_sct <- function(seurat_objects, project_name, do.umap = TRUE, do.norm.rna = TRUE) {
    hprcc::init_multisession()
    future::plan("multisession")

    message("Loading objects to list...\n")
    # Parallel loading of Seurat objects
    seurat_list <- future.apply::future_lapply(seurat_objects, function(obj_path) {
        message("Loading ", obj_path, "...")
        loaded_object <- load_seurat(obj_path)
        return(loaded_object)
    }, future.seed = 42)
    gc()
    # Extract sample IDs and set them as names of the list elements
    names(seurat_list) <- sapply(seurat_list, function(obj) {
        sample_id <- Seurat::Project(obj)
        return(sample_id)
    })

    message("Loaded ", length(seurat_list), " objects")
    message("List size: ", format(object.size(seurat_list), units = "Gb"))
    message("Merging objects...\n")
    object <- merge(seurat_list[[1]], seurat_list[-1])
    rm(seurat_list)
    gc()
    message("Merged size: ", format(object.size(object), units = "Gb"))

    message("Normalising merged object with PrepSCTFindMarkers()...\n")
    # normalise to median UMI
    future::plan("multisession", workers = 4)
    object<-Seurat::PrepSCTFindMarkers(object)

    message("Setting variable features to scale.data...\n")
    # set variable features
    # https://github.com/satijalab/seurat/issues/2814
    Seurat::VariableFeatures(object[["SCT"]]) <- rownames(object[["SCT"]]@scale.data)

    if (do.umap) {
        object <- object |>
            Seurat::RunPCA(verbose = TRUE) |>
            harmony::RunHarmony(group.by.vars = c("orig.ident"), verbose = TRUE) |>
            Seurat::FindNeighbors(reduction = "harmony", verbose = TRUE) |>
            Seurat::FindClusters(resolution = 0.15, method = "igraph", cluster.name = "sct_clusters", verbose = TRUE) |>
            Seurat::RunUMAP(reduction = "harmony", dims = 1:50, verbose = TRUE)
    }
    
    if (do.norm.rna) {
        message("Normalising RNA assay counts with NormalizeData()...\n")
        object <- object |>
            Seurat::NormalizeData(assay = "RNA", verbose = TRUE)
    }
    Seurat::DefaultAssay(object) <- "RNA"
    Seurat::Project(object) <- project_name
    object_path <- glue::glue("{analysis_cache}/seurat_merged/{project_name}.qs")
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    message("Saving Seurat object...")
    qs::qsave(object, file = object_path)
    message("Done.")
    return(object_path)
}
