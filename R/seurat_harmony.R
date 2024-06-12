
#' Perform Harmony Integration on a Seurat Object
#'
#' This function performs PCA, Harmony integration, neighbor finding, clustering, and UMAP on a Seurat object.
#' The Seurat object must have variable features and some form of normalization applied.
#'
#' @param seurat_object A Seurat object to be processed.
#' @param assay The assay to be used for the analysis.
#' @param group_by_vars A character vector specifying the variables to group by for Harmony integration. Default is "batch".
#' @param res The resolution parameter for clustering. Default is 0.3.
#' @param dims Dimensions to use for UMAP. Default is 1:50.
#' @return The file path where the processed Seurat object is saved.
#' @examples
#' \dontrun{
#' seurat_harmony(seurat_object, assay = "RNA", group_by_vars = "batch", res = 0.3, dims = 1:50)
#' }
#' @importFrom Seurat Project DefaultAssay RunPCA FindNeighbors FindClusters RunUMAP
#' @importFrom harmony RunHarmony
#' @importFrom glue glue
#' @importFrom qs qsave
#' @importFrom future plan
#' @importFrom hprcc init_multisession
#' @importFrom tibble rownames_to_column
#' @export
seurat_harmony <- function(seurat_object, assay, group_by_vars = "batch", res = 0.3, dims = 1:50) {
    future::plan("multisession")
    options(future.globals.maxSize = hprcc::slurm_allocation()$Memory_GB * 1024^3)

    message("This is seurat_harmony()...\n")
    seurat_object <- load_seurat(seurat_object)
    
    project_name <- Seurat::Project(seurat_object)

    Seurat::DefaultAssay(seurat_object) <- assay

    seurat_object <- seurat_object |>
            Seurat::RunPCA(verbose = TRUE) |>
            harmony::RunHarmony(group.by.vars = group_by_vars, verbose = TRUE) |>
            Seurat::FindNeighbors(reduction = "harmony", verbose = TRUE) |>
            Seurat::FindClusters(resolution = res, method = "igraph", cluster.name = "harmony_clusters", verbose = TRUE) |>
            Seurat::RunUMAP(reduction = "harmony", dims = dims, verbose = TRUE)

    new_project_name <- glue::glue("harmony__{paste(group_by_vars, collapse = '_')}__{project_name}")
    
    Seurat::Project(seurat_object) <- new_project_name
    #--------------------------------------------------------------------------------
    # Save results
    #object_path <- (glue::glue("{analysis_cache}/harmony_out/{new_project_name}.qs"))
    #save_results(seurat_object, object_path)
    return(seurat_object)
}
