#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param seurat_object
#' @param cell_metadata
#' @param sample_metadata
#' @return
#' @author Denis O'Meally
#' @export
seurat_annotate_cells <- function(seurat_object, cell_metadata) {
    # if (is.null(seurat_object)) seurat_object <- tar_read(integrated_seurat_sketch_750)
    # if (is.null(cell_metadata)) cell_metadata <- tar_read(aggregated_cell_annot_csv)
    # if (is.null(sample_metadata)) sample_metadata <- tar_read(pancdb_metadata_gt)

    #    seurat_object <- targets::tar_read(integrated_seurat_sketch_750)
    #    cell_metadata <- targets::tar_read(aggregated_cell_annot_csv)
    #    sample_metadata <- targets::tar_read(pancdb_metadata_agg)

    seurat_object <- load_seurat(seurat_object)
    assay <- Seurat::DefaultAssay(seurat_object)
    cell_metadata <- cell_metadata %>%
        dplyr::mutate(
            cell_cycle = Phase,
            tosti_cell_type = tosti_etalcell_type,
            elgamal_cell_type = HPAPcell_type,
            cell_type = dplyr::case_when(
                HPAPcell_type %in% c("Alpha", "Cycling Alpha", "Beta", "Alpha+Beta", "Delta", "Gamma+Epsilon", "Epsilon") ~ HPAPcell_type,
                HPAPcell_type %in% c("Acinar") ~ HPAPcell_type,
                HPAPcell_type %in% c("Ductal", "MUC5B+ Ductal") ~ HPAPcell_type,
                TRUE ~ "Other"),
            hpap_cell_type = hpap_celltype,
            hpap_clusters = cluster_labels,
            seurat_clusters = manual_clusters,
            xbp1u_psi = xbp1_psi,
            spliced_ratio_INS = gene_ratio_INS,
            spliced_ratio_XBP1 = gene_ratio_XBP1,
            spliced_ratio_GAPDH = gene_ratio_GAPDH) %>%
        dplyr::select(
            cell,
            cell_cycle,
            tosti_cell_type,
            elgamal_cell_type,
            cell_type,
            hpap_cell_type,
            hpap_clusters,
            seurat_clusters,
            xbp1u_psi,
            INS_hk,
            spliced_ratio_INS,
            spliced_ratio_XBP1,
            spliced_ratio_GAPDH,
            contains("_UCell")) %>%
        tibble::column_to_rownames(var = "cell")

    seurat_object <- Seurat::AddMetaData(
        object = seurat_object,
        metadata = cell_metadata
    )

    
    # Add INS+ cell_type annotations
    require(tidyseurat)
    INS_hk_threshold <- find_expression_valley(seurat_object, "INS_hk")$valley
    seurat_object <- seurat_object |>
        mutate(
            cell_type_extra = ifelse(INS_hk >= INS_hk_threshold, "INS+", "")
        ) 

    return(seurat_object)
}
