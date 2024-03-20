#' Plot CellBender analysis results
#'
#' This function generates a QC plot of CellBender analysis results using Seurat object, including UMAPs showing the gene
#' whose expression differs most after correction, a plot of all genes whose expression by more that the specified threshold,
#' and a ClusterProfilter enrichment plot for that list of genes, using the specified cell marker data.
#' 
#' @param seurat_object The Seurat object or path to the Seurat object containing the scRNA-seq data.
#' @param pct_diff_threshold The percentage difference threshold for selecting genes that differ most after correction (default: 5).
#' @param gene_of_interest A single gene of interest to highlight in the plot (default: "GAPDH").
#' @param cell_marker_data The path to the cell marker data file (default: "{analysis_cache}/data/cell_marker_data.txt").
#'
#' @return the file path of the generated plot
#'
#' @import scCustomize
#' @import ggplot2
#' @import patchwork
#' @import clusterProfiler
#'
#' @examples
#' plot_cellbender(seurat_object = my_seurat_object, pct_diff_threshold = 5, gene_of_interest = "GAPDH", cell_marker_data = "path/to/cell_marker_data.txt")
#'
seurat_plot_cellbender <- function(seurat_object, gene_of_interest = "GAPDH", pct_diff_threshold = 5, cell_marker_data = glue::glue("{analysis_cache}/data/cell_marker_data.txt")) {
    suppressPackageStartupMessages({       
    require(scCustomize)
    require(ggplot2)
    require(patchwork)
    require(clusterProfiler)
    require(enrichplot)
    })

    plots <- list()
    seurat_object<-load_seurat(seurat_object)
    sample_id <- seurat_object$orig.ident[1]
    seurat_object <- NormalizeData(seurat_object, assay = "RNA", verbose = FALSE) |>
        NormalizeData(assay = "RAW", verbose = FALSE) |>
        .clusterData() # from seurat_ddqc.R
    if (!gene_of_interest %in% rownames(seurat_object)) {
        message <- glue::glue("{gene_of_interest}, not found. Plotting GAPDH instead")
        gene_of_interest <- "GAPDH"
        message(message) } else {message <- NULL}
    
    feature_diff <- CellBender_Feature_Diff(seurat_object, raw_assay = "RAW", cell_bender_assay = "RNA")

    plots$max_contaminant_plot <- FeaturePlot_DualAssay(
        seurat_object, features = rownames(feature_diff)[1], assay1 = "RAW", assay2 = "RNA", colors_use = viridis_plasma_dark_high, na_color = "gray", num_col = 2) &
            theme_classic() &
            theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    plots$gene_of_interest_plot <- FeaturePlot_DualAssay(
        seurat_object, features = gene_of_interest, assay1 = "RAW", assay2 = "RNA", colors_use = viridis_plasma_dark_high, na_color = "gray", num_col = 2) & NoLegend() &
            theme_classic() & labs(subtitle = message) &
            theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    top_genes <- feature_diff %>%
      as_tibble(rownames = "Gene") %>%
      filter(Pct_Diff >= 5) %>%
      mutate(Weighted_Score = Count_Diff * Pct_Diff) %>%
      arrange(desc(Weighted_Score)) %>%
      dplyr::select(Gene, Weighted_Score) %>%
      {setNames(.$Weighted_Score, .$Gene)}
    if (!is.data.frame(cell_marker_data)) {cell_marker_data <- vroom::vroom(cell_marker_data, show_col_types = FALSE)}
    cell_types <- cell_marker_data |> dplyr::select(cellName, geneSymbol) |> mutate(geneSymbol = strsplit(geneSymbol, ', ')) |> tidyr::unnest(cols = c(geneSymbol))
    enriched_cell_types <- clusterProfiler::GSEA(top_genes, TERM2GENE = cell_types, scoreType = "pos")
    diff_plot <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = pct_diff_threshold) + theme_classic()
    enriched_plot <- if (nrow(enriched_cell_types@result) == 0) grid::textGrob('No enriched cell types found') else barplot(enriched_cell_types, showCategory=10) 
    plots$bottom_row <-patchwork::wrap_plots(diff_plot, enriched_plot, ncol = 2)
    combined_plot <- patchwork::wrap_plots(plots, nrow = 3) &
        patchwork::plot_annotation(title = glue::glue("{sample_id} cellbender plot"), theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
    combined_plot_path <- glue::glue("{analysis_cache}/cellbender_out/{sample_id}_plot.png")
    dir.create(dirname(combined_plot_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = combined_plot_path, plot = combined_plot, width = 10, height = 20)
    return(combined_plot_path)

# http://xteam.xbio.top/CellMarker/download/all_cell_markers.txt
# http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt
# http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt
# http://xteam.xbio.top/CellMarker/download/Single_cell_markers.txt

}
