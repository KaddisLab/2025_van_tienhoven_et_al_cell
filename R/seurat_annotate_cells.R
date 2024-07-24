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
            contains("gene_ratio_INS"),
            contains("spliced_counts"),
            contains("unspliced_counts"),
            contains("_UCell")) %>%
        tibble::column_to_rownames(var = "cell")

    seurat_object <- Seurat::AddMetaData(
        object = seurat_object,
        metadata = cell_metadata
    )

    # Add INS+ cell_type_extra ------------------------------
    require(tidyseurat)
    INS_hk_threshold <- find_valley(seurat_object[[]]$INS_hk)$valley
    seurat_object <- seurat_object |>
        mutate(
            cell_type_extra = ifelse(INS_hk >= INS_hk_threshold, "INS+", "")
        )
    ## Add simple Gene Set scores ------------------------------
    # Function to calculate signature scores
    calculate_signature_score <- function(seurat_object, signature_genes) {
        require(tidyseurat)
        require(dplyr)
        require(scales)

        positive_genes <- gsub("\\+$", "", signature_genes[!grepl("-$", signature_genes)])
        negative_genes <- gsub("-$", "", signature_genes[grepl("-$", signature_genes)])

        seurat_object <- seurat_object %>%
            join_features(
                features = c(positive_genes, negative_genes),
                shape = "wide",
                slot = "data",
                assay = assay
            )

        score <- seurat_object %>%
            mutate(
                pos_score = rowSums(across(any_of(positive_genes))),
                neg_score = if (length(negative_genes) > 0) rowSums(across(any_of(negative_genes))) else 0,
                scaled_pos_score = scales::rescale(pos_score, to = c(0, 1)),
                scaled_neg_score = if (length(negative_genes) > 0) scales::rescale(neg_score, to = c(0, 1)) else 0,
                signature_score = scaled_pos_score - scaled_neg_score
            ) %>%
            pull(signature_score)

        return(score)
    }

    # Calculate and add stress scores
    stress_signatures <- c(
        "chronic_er_stress", "active_er_stress", "islet_er_stress",
        "islet_stress", "cellular_stress", "core_upr_stress", "msigdb_upr_stress"
    )

    for (sig_name in stress_signatures) {
        message(paste("Calculating signature score for", sig_name))
        seurat_object[[paste0(sig_name, "_score")]] <- calculate_signature_score(seurat_object, signatures[[sig_name]])
    }

    return(seurat_object)
}
