#' Annotate cells in a Seurat object with metadata and signature scores
#'
#' Integrates cell metadata and calculates gene signature scores for stress-related pathways.
#' Handles cell type classification standardization and adds splicing metrics.
#'
#' @param seurat_object A Seurat object or file path to Seurat object
#' @param cell_metadata data.frame containing cell annotations with columns:
#'   - Phase (cell cycle phase)
#'   - HPAPcell_type (original cell type labels)
#'   - xbp1_psi (XBP1 splicing ratio)
#'   - counts_ratio_* (splicing metrics)
#'   - INS_hk (INS housekeeping counts)
#' @param signatures Named list of gene signatures containing named character vectors of genes
#'
#' @return Modified Seurat object with added metadata:
#'   - Standardized cell_type classifications
#'   - Original elgamal_cell_type labels
#'   - Splicing ratios (xbp1u_psi, spliced_ratio_*)
#'   - Stress signature scores (*_stress_score)
#' 
#'
#' @importFrom Seurat AddMetaData DefaultAssay
#' @importFrom dplyr mutate case_when select contains
#' @importFrom tibble column_to_rownames
#' @export
seurat_annotate_cells <- function(seurat_object, cell_metadata) {
    seurat_object <- load_seurat(seurat_object)
    assay <- Seurat::DefaultAssay(seurat_object)
    cell_metadata <- cell_metadata |>
        dplyr::mutate(
            cell_cycle = Phase,
            elgamal_cell_type = HPAPcell_type,
            cell_type = dplyr::case_when(
                HPAPcell_type %in% c("Alpha", "Cycling Alpha", "Beta", "Alpha+Beta", "Delta", "Gamma+Epsilon", "Epsilon") ~ HPAPcell_type,
                HPAPcell_type %in% c("Acinar") ~ HPAPcell_type,
                HPAPcell_type %in% c("Ductal", "MUC5B+ Ductal") ~ HPAPcell_type,
                TRUE ~ "Other"
            ),
            xbp1u_psi = xbp1_psi,
            spliced_ratio_INS = counts_ratio_INS,
            spliced_ratio_XBP1 = counts_ratio_XBP1,
            spliced_ratio_GAPDH = counts_ratio_GAPDH
        ) |>
        dplyr::select(
            cell,
            cell_cycle,
            elgamal_cell_type,
            cell_type,
            xbp1u_psi,
            INS_hk,
            contains("counts_ratio"),
            contains("spliced_counts"),
            contains("unspliced_counts"),
            contains("_UCell")
        ) |>
        tibble::column_to_rownames(var = "cell")

    seurat_object <- Seurat::AddMetaData(
        object = seurat_object,
        metadata = cell_metadata
    )

    ## Add simple Gene Set scores ------------------------------
    # Calculate and add stress scores
    stress_signatures <- names(signatures)

    for (signature_name in stress_signatures) {
        message(paste("Calculating signature score for", signature_name))
        seurat_object[[paste0(signature_name, "_score")]] <- calculate_signature_score(seurat_object, signatures[[signature_name]])
    }

    return(seurat_object)
}
