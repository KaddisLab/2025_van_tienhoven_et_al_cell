#' Merge Seurat Objects for BPCells
#'
#' This function takes a set of QCed Seurat objects, PancDB metadata, and additional metadata files
#' (protected cohort, azimuth mapped Seurat objects, cell cycle CSV, and Tosti cell type CSV),
#' merges them into a single Seurat object with enriched metadata, and saves the merged object to disk.
#' The merging process includes metadata integration, sample metadata addition, and cell metadata enrichment.
#'
#' @param ddqc_seurat_objects Vector of file paths to Seurat objects that passed data-driven quality control.
#' @param pancdb_metadata DataFrame containing PancDB metadata.
#' @param protected_cohort DataFrame specifying the protected cohort information.
#' @param azimuth_mapped_seurat_objects Vector of file paths to CSV files with Azimuth mapped Seurat object data.
#' @param cell_cycle_csv Vector of file paths to CSV files containing cell cycle phase data.
#' @param tosti_cell_type_csv Vector of file paths to CSV files containing Tosti cell type annotations.
#'
#' @return A string with the file path to the saved RDS file containing the merged Seurat object.
#'
#' @examples
#' merged_seurat_path <- merge_seurat_bpcells(
#'   ddqc_seurat_objects = c("path/to/ddqc_seurat_object1.qs", "path/to/ddqc_seurat_object2.qs"),
#'   pancdb_metadata = pancdb_metadata_df,
#'   protected_cohort = protected_cohort_df,
#'   azimuth_mapped_seurat_objects = c("path/to/azimuth1.csv", "path/to/azimuth2.csv"),
#'   cell_cycle_csv = c("path/to/cell_cycle1.csv", "path/to/cell_cycle2.csv"),
#'   tosti_cell_type_csv = c("path/to/tosti_cell_type1.csv", "path/to/tosti_cell_type2.csv")
#' )
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom dplyr filter mutate left_join select
#' @importFrom readr read_csv
#' @importFrom progress progress_bar
#' @export
#' @import tidyverse Seurat scCustomize BPCells progress
merge_seurat_bpcells <- function(ddqc_seurat_objects, pancdb_metadata, protected_cohort, azimuth_mapped_seurat_objects, cell_cycle_csv, tosti_cell_type_csv) {
    suppressPackageStartupMessages({
        require(tidyverse)
        require(Seurat)
        require(scCustomize)
        require(BPCells)
        require(progress)
    })

    # Sample metadata
    # rs3842752 annotation
    pancdb_metadata$protected <- pancdb_metadata$donor_id %in% protected_cohort$sample_id
    # 10X libraries that passed QC
    pancdb_metadata <- pancdb_metadata |>
        dplyr::filter(str_detect(reagent_kit, "10X") & !str_detect(donor_id, failed_qc_donor_ids)) |>
        dplyr::mutate(batch = as.integer(as.factor(reagent_kit)))
    ddqc_seurat_object_paths <- data.frame(
        "qs_path" = ddqc_seurat_objects,
        "donor_id" = gsub(".*/(HPAP-\\d+)_ddqc\\.qs", "\\1", ddqc_seurat_objects)
    )
    pancdb_metadata <- pancdb_metadata |> left_join(ddqc_seurat_object_paths, by = "donor_id")

    seurat_paths <- pancdb_metadata |> pull(qs_path)
    meta_data <- select(pancdb_metadata, c(donor_id, protected, batch, sample_sex, sample_age))

    # Initialise progress bar
    pb <- progress_bar$new(
        format = "  Loading [:bar] :percent eta: :eta (Loading :current of :total)",
        total = length(seurat_paths), clear = FALSE, width = 60
    )

    data_list <- c()
    metadata_list <- c()

    for (path in seurat_paths) {
        seurat_object <- load_seurat(path)
        # add sample metadata
        seurat_object <- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, join_by_seurat = "orig.ident", join_by_meta = "donor_id")
        sample_id <- seurat_object$orig.ident[1]
        # add cell metadata 
        azimuth_data <- sub("\\.qs", "\\.csv", grep(sample_id, azimuth_mapped_seurat_objects, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        cell_cycle_data <- sub("\\.qs", "\\.csv", grep(sample_id, cell_cycle_csv, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        tosti_cell_type_data <- sub("\\.qs", "\\.csv", grep(sample_id, tosti_cell_type_csv, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        cell_metadata <- azimuth_data |>
            left_join(cell_cycle_data, by = "cell") |>
            left_join(tosti_cell_type_data, by = "cell") |>
            select(c(cell, azimuth_label = "predicted.annotation.l1", cell_cycle_phase = "Phase", tosti_cell_type = "labels")) |>
            column_to_rownames(var = "cell")
        seurat_object <- Seurat::AddMetaData(seurat_object, cell_metadata)

        bp_dir <- glue::glue("{analysis_cache}/bpcells_out/{sample_id}")

        iterable_matrix <- as(seurat_object[["RNA"]]$counts, "IterableMatrix")
        # Convert counts to BPCells / load the BPCells folder
        if (!dir.exists(bp_dir)) {
            write_matrix_dir(mat = iterable_matrix, dir = bp_dir, compress = TRUE)
        }
        mat <- open_matrix_dir(dir = bp_dir)

        # Get mat and metadata for later merging
        data_list[[sample_id]] <- mat
        metadata_list[[sample_id]] <- seurat_object[[]]
        # Increment progress bar after each iteration
        pb$tick()
    }

    # Build the merged seurat object
    metadata <- Reduce(rbind, metadata_list)

    merged_seurat <- Seurat::CreateSeuratObject(counts = data_list, meta.data = metadata)
    
    merged_seurat_path <- glue::glue("{analysis_cache}/data/merged_seurat_bp.rds")

    saveRDS(merged_seurat, file = merged_seurat_path)

    return(merged_seurat_path)
}
