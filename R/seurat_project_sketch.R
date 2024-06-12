#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param integrated_seurat_sketch_750
#' @return
#' @author Denis O'Meally
#' @export
seurat_project_sketch <- function(
    seurat_object,
    assay = "RNA",
    reduction = "harmony",
    full_reduction = "harmony.full",
    sketched_assay = "sketch",
    sketched_reduction = "harmony.full",
    umap_model = "umap",
    dims = 1:50,
    reduction_key = "harmony_full_UMAP",
    reduction_name = "harmony_full_UMAP",
    refdata = list(cluster_full = "seurat_clusters")) {

    seurat_object <- load_seurat(seurat_object)
    new_project_name <- Seurat::Project(seurat_object) |> str_replace("_sketch_", "_full_")
    
    seurat_object <- ProjectIntegration(
        seurat_object,
        sketched.assay = sketched_assay,
        assay = assay,
        reduction = "harmony"
    )

    # catch errors with try()
    tryCatch({
        seurat_object <- ProjectData(
            object = seurat_object,
            sketched.assay = sketched_assay,
            assay = assay,
            sketched.reduction = full_reduction,
            full.reduction = full_reduction,
            dims = dims,
            refdata = refdata
        )
    })
    
    DefaultAssay(seurat_object) <- "RNA"

    seurat_object <- RunUMAP(seurat_object,
        reduction = full_reduction, dims = dims, reduction.name = reduction_name,
        reduction.key = reduction_key
    )

    Seurat::Project(seurat_object) <- new_project_name
    
    return(seurat_object)

}
