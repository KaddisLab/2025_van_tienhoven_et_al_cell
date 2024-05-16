#' Transfer Cell Type Annotation from Reference to Query Seurat Object
#'
#' This function applies cell type annotations from a reference Seurat object
#' to a query Seurat object using the Azimuth package. It saves the updated
#' query object with transferred annotations, as well as a UMAP plot and a
#' CSV file containing the predicted annotations and mapping scores.
#'
#' @param query_seurat_object A Seurat object or a path to an RDS file containing the query Seurat object.
#' @param azimuth_reference_path A path to the directory containing the reference Seurat object for Azimuth.
#'
#' @return The file path of the updated query Seurat object with transferred annotations.
#' @export
seurat_azimuth <- function(query_seurat_object, azimuth_reference_path) {
    # https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
    set.seed(42)
    require(Seurat)
    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    query_seurat_object <- load_seurat(query_seurat_object) 

    azimuth_reference <- dirname(azimuth_reference_path[1])
    
    azimuth <- Azimuth::RunAzimuth(query_seurat_object, reference = azimuth_reference)

    sample_id <- query_seurat_object@meta.data$orig.ident[1]

    plot_azimuth <- DimPlot(azimuth, group.by = "predicted.annotation.l1", cols = cell_type_palette, label = TRUE, label.size = 4, repel = TRUE, reduction = "ref.umap", shuffle = TRUE) +
        NoLegend() + 
        labs(title = glue::glue("{sample_id} projected into Azimuth ref")) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    
    plot_path <- glue::glue("{analysis_cache}/azimuth_out/{sample_id}_azimuth.png")
    
    dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(plot_path, plot = plot_azimuth, width = 10, height = 10, dpi = 300)

    predicted_azimuth_data <- azimuth[[]] |> tibble::rownames_to_column(var = "cell") |> select("cell", matches("predicted"), "mapping.score")
    predicted_azimuth_data_path <- glue::glue("{analysis_cache}/azimuth_out/{sample_id}_azimuth.csv")
    
    predicted_azimuth_data |> write.csv(predicted_azimuth_data_path, row.names = FALSE, quote = FALSE)

    azimuth_path <- glue::glue("{analysis_cache}/azimuth_out/{sample_id}_azimuth.qs")

    qs::qsave(azimuth, azimuth_path)

    return(azimuth_path)
}
