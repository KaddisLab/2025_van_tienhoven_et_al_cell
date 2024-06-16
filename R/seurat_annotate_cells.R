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

    # Add Beta_like & epsilon cell_type annotations
    require(tidyseurat)
    marker_genes <- c("GHRL", "SOX15", "FEV", "INS", "GCG", "SST", "PPY")

    seurat_object <- seurat_object |>
        join_features(
            features = marker_genes,
            shape = "wide",
            assay = assay
        ) |>
        mutate(
            epsilon_score =
                scales::rescale(GHRL + SOX15 + FEV, to = c(0, 1)) - scales::rescale(INS + GCG + SST + PPY, to = c(0, 1)),
            cell_type_extra = case_when(
                INS >= 5.8 ~ "INS+",
                epsilon_score > 0.3 ~ "Epsilon",
                TRUE ~ ""
            )
        )
    # # add Stress scores
    # # Add upr_score and er_stress_score to the Seurat object
    # seurat_object <- seurat_object |>
    #     join_features(
    #         features = c(upr_genes, er_stress_genes),
    #         shape = "wide",
    #         assay = assay
    #     ) |>
    #     mutate(
    #         upr_score = scales::rescale(rowSums(across(all_of(upr_genes))), to = c(0, 1)),
    #         er_stress_score = scales::rescale(rowSums(across(all_of(er_stress_genes))), to = c(0, 1))
    #     ) |>
    #     # clean up
    #     select(-all_of(c(marker_genes, upr_genes, er_stress_genes)))

    # project_name <- Seurat::Project(object = seurat_object)
    # seurat_object_path <- glue::glue("{analysis_cache}/data/{project_name}_annotated.qs")
    # qs::qsave(seurat_object, seurat_object_path)
    # message("Saved seurat object")

    return(seurat_object)
}
