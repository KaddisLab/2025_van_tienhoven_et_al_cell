get_pancdb_metadata <- function() {
    #read excel file
    donor_metadata <- readxl::read_excel(glue::glue("{analysis_cache}/data/metadata/PancDB_Donors.xlsx"))
    scrnaseq_metadata <- readxl::read_excel(glue::glue("{analysis_cache}/data/metadata/PancDB_scRNA-seq_metadata_2023-12-22.xlsx"))
    metadata <- left_join(donor_metadata, scrnaseq_metadata, by = "DonorID") |>
        janitor::clean_names() |>
        mutate(sample_age = readr::parse_number(sample_age),
             diabetes_status = case_when(
               stringr::str_detect(disease_status, "T1DM") ~ "T1DM",
               stringr::str_detect(disease_status, "T2DM") ~ "T2DM",
               stringr::str_detect(disease_status, "No HX DIAB") ~ "NODM",
               TRUE ~ NA_character_
             ),
             technology = case_when(
               stringr::str_detect(reagent_kit, "V2") ~ "10XV2",
               stringr::str_detect(reagent_kit, "V3") ~ "10XV3",
               stringr::str_detect(reagent_kit, "C1") ~ "SMARTSEQ3",
               TRUE ~ NA_character_
             ),
             sample_race = dplyr::case_when(
                stringr::str_detect(sample_ethnicity, "African-Am") ~ "AFA",
                stringr::str_detect(sample_ethnicity, "Cauc") ~ "CAU",
                stringr::str_detect(sample_ethnicity, "Hisp") ~ "HIS",
                TRUE ~ NA_character_
            )
        )
    return(metadata)
}
