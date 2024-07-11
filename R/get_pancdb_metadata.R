get_pancdb_metadata <- function() {
    #read excel file
    donor_metadata <- readxl::read_excel(glue::glue("{analysis_cache}/data/metadata/PancDB_Donors.xlsx")) |>
        left_join(readxl::read_excel(glue::glue("{analysis_cache}/data/metadata/hpap_donor_types.xlsx")), by = join_by(DonorID))
    scrnaseq_metadata <- readxl::read_excel(glue::glue("{analysis_cache}/data/metadata/PancDB_scRNA-seq_metadata_2023-12-22.xlsx"))
    metadata <- left_join(donor_metadata, scrnaseq_metadata, by = "DonorID") |>
        janitor::clean_names() |>
        mutate(
            sample_age = readr::parse_number(sample_age),
            diabetes_status = case_when(
                stringr::str_detect(disease_status, "T1DM") ~ "T1DM",
                stringr::str_detect(disease_status, "T2DM") ~ "T2DM",
                stringr::str_detect(ab_positive, "Yes") ~ "AABP",
                TRUE ~ "NODM"
            ),
            technology = case_when(
                stringr::str_detect(reagent_kit, "v2") ~ "10XV2",
                stringr::str_detect(reagent_kit, "v3.1") ~ "10XV3.1",
                stringr::str_detect(reagent_kit, "v3") ~ "10XV3",
                stringr::str_detect(reagent_kit, "C1") ~ "SMARTSEQ3",
                TRUE ~ NA_character_
            ),
            sample_race = dplyr::case_when(
                stringr::str_detect(sample_ethnicity, "African-Am") ~ "AFA",
                stringr::str_detect(sample_ethnicity, "Cauc") ~ "CAU",
                stringr::str_detect(sample_ethnicity, "Hisp") ~ "HIS",
                TRUE ~ NA_character_
            )
        ) |>
        dplyr::filter(str_detect(reagent_kit, "10X")) |>
        dplyr::mutate(
            batch = as.integer(as.factor(reagent_kit)),
            sample_id = donor_id,
            sample_name = donor_id,
            tissue_source = if_else(str_detect(tissue_source, "Penn"), "UPenn", "nPod")
        )
    
    return(metadata)
}

