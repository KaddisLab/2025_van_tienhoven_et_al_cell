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
    cell_type_palette_orig <- c(
  # Endocrine
  "Alpha" = "#2ECC71",              # Green
  "alpha" = "#2ECC71",              # Green
  "Beta" = "#3498DB",               # Blue
  "Beta-like" = "#15bfd9",               # Blue
  "beta" = "#3498DB",               # Blue
  "Beta-like" = "#15bfd9",               # Blue
  "Delta" = "#1ABC9C",              # Teal
  "delta" = "#1ABC9C",              # Teal
  "Gamma" = "#16A085",              # Dark Teal
  "gamma" = "#16A085",              # Dark Teal
  "PP_Gamma" = "#16A085",              # Dark Teal
  "epsilon" = "#27AE60",            # Emerald
  "Epsilon" = "#27AE60",            # Emerald
  
  # Exocrine
  "Acinar-s" = "#E74C3C",           # Red
  "Acinar-i" = "#E67E22",           # Orange
  "Acinar-REG+" = "#F39C12",        # Amber
  "acinar" = "#E74C3C",             # Red
  "Acinar" = "#E74C3C",             # Red
  "Ductal" = "#9B59B6",             # Purple
  "ductal" = "#9B59B6",             # Purple
  "MUC5B+ Ductal" = "#8E44AD",      # Dark Purple
  
  # Immune
  "Macrophage" = "brown",         # Dark Grey
  "macrophage" = "brown",         # Dark Grey
  "immune" = "brown",             # Pink
  "Immune" = "brown",             # Pink
  
  # Other
  "Other" = "#314c4e",        # 
  "Endothelial" = "#314c4e",        # 
  "endothelial" = "#314c4e",        # 
  "Activated Stellate" = "#F1C40F", # Yellow
  "Stellates_Mesenchymal" = "#F1C40F", # Yellow
  "activated_stellate" = "#F1C40F", # Yellow
  "Quiescent Stellate" = "#FDFD96", # Light Yellow
  "quiescent_stellate" = "#FDFD96", # Light Yellow
  "Schwann" = "#2C3E50",            # Dark blue
  "schwann" = "#2C3E50",            # Dark blue
  "cycling" <- "#FF7F50",
  "Unknown" <- "cornsilk2" 
)

    hprcc::init_multisession()

    query_seurat_object <- load_seurat(query_seurat_object) 

    azimuth_reference <- dirname(azimuth_reference_path[1])
    
    azimuth <- Azimuth::RunAzimuth(query_seurat_object, reference = azimuth_reference)

    sample_id <- query_seurat_object|> Seurat::Project()

    plot_azimuth <- DimPlot(azimuth, group.by = "predicted.annotation.l1", cols = cell_type_palette_orig, label = TRUE, label.size = 4, repel = TRUE, reduction = "ref.umap", shuffle = TRUE) +
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
