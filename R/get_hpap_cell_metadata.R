get_hpap_cell_metadata <- function() {
    hpap <- load_seurat(glue::glue("{analysis_cache}/data/T1D_T2D_20220428.rds"))
    hpap_metadata_csv_path <- glue::glue("{analysis_cache}/data/T1D_T2D_20220428_cell_metadata.csv")

    hpap_metadata <- hpap[[]] |>
        select(hpap_celltype = "cell_type", hpap_id) |>
        tibble::rownames_to_column("cell") |>
        dplyr::mutate(
            hpap_id = stringr::str_replace(hpap_id, "HPAP", "HPAP-"),
            cell = stringr::str_replace(cell, "(?<=-)[0-9]+", "1"),
            cell = glue::glue("{hpap_id}_{cell}")
        ) |>
        select(-hpap_id)

    write.csv(hpap_metadata, file = hpap_metadata_csv_path, row.names = FALSE, quote = FALSE)
    return(hpap_metadata_csv_path)
}

import_hpap_seurat <- function() {
    hpap <- load_seurat(glue::glue("{analysis_cache}/data/T1D_T2D_20220428.rds"))
    hpap_metadata_csv_path <- glue::glue("{analysis_cache}/data/T1D_T2D_20220428_cell_metadata.csv")

return(hpap)
}