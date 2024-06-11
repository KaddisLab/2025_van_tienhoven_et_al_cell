#' Transfer Cell Type Annotation from Reference to Query Seurat Object
#'
#' This function applies cell type annotations from a reference Seurat object
#' to a test Seurat object using SingleR. It supports optional normalization
#' of the reference object and uses differential expression methods to assist
#' with the annotation transfer.
#'
#' @param query_seurat_object a seurat object or path to the QS file for the test Seurat object.
#' @param ref_seurat_object a seurat object or path to the QS file for the reference Seurat object.
#' @param reduction_model the dimensionality reduction model to use for the projection.
#' @param ref_data a list of reference data to project onto the query object.

#' @export
seurat_project_into_ref <- function(query_seurat_object, ref_seurat_object, reduction_model = "umap", ref_data = list(cell_type = "Cluster")) {
    # https://satijalab.org/seurat/articles/integration_mapping
    set.seed(42)
    require(Seurat)
    options(parallelly.availableCores.methods = "Slurm")
    future::plan("multisession")
    options(future.globals.maxSize = hprcc::slurm_allocation()$Memory_GB / hprcc::slurm_allocation()$CPUs * 1024^3)

    query_seurat_object <- load_seurat(query_seurat_object) 

    ref_seurat_object <- load_seurat(ref_seurat_object) 

    if (!isNormalized(query_seurat_object)) {
        query_seurat_object <- Seurat::NormalizeData(query_seurat_object)
    }

    anchors <- Seurat::FindTransferAnchors(
        reference = ref_seurat_object,
        query = query_seurat_object,
        reference.reduction = "pca")

    updated_query <- Seurat::MapQuery(
        anchorset = anchors,
        reference = ref_seurat_object,
        query = query_seurat_object,
        refdata = ref_data,
        reference.reduction = "pca",
        reduction.model = reduction_model)

    sample_id <- query_seurat_object@meta.data$orig.ident[1]

    ref_id <- ref_seurat_object@project.name

    plot_updated_query <- DimPlot(updated_query, group.by = "predicted.cell_type", cols = custom_palette, label = TRUE, label.size = 4, repel = TRUE, reduction = "ref.umap", shuffle = TRUE) +
        NoLegend() + 
        labs(title = glue::glue("{sample_id} projected into {ref_id} {reduction_model}")) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    
    plot_path <- glue::glue("{analysis_cache}/ref_projection_out/{sample_id}_{reduction_model}_{ref_id}.png")
    
    dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(plot_path, plot = plot_updated_query, width = 10, height = 10, dpi = 300)

    predicted_ref_data <- updated_query[[]] |> as.data.frame() |> tibble::rownames_to_column(var = "cell") |> select("cell", matches("predicted")) |> tibble::as_tibble()
    predicted_ref_data_path <- glue::glue("{analysis_cache}/ref_projection_out/{sample_id}_{reduction_model}_{ref_id}_predictions.csv")
    
    predicted_ref_data |> write.csv(predicted_ref_data_path, row.names = FALSE, quote = FALSE)

    updated_query_path <- glue::glue("{analysis_cache}/ref_projection_out/{sample_id}_{reduction_model}_{ref_id}.ps")

    qs::qsave(updated_query, updated_query_path)

    return(updated_query_path)
}
