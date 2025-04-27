#' Rename and symlink for Cellranger fastq files
#'
#' This function renames the samples for Cellranger and creates symlinks for the fastq files
#'
#' @param original_paths A data frame containing the original file paths for the fastq files
#' @return A data.frame with paths to fastq symlinks, formatted for a nfcore sample_sheet
#'
#' @examples
#' rename_sample_sheet(original_sample_sheet)
#'
#' @export
#'
#' @seealso
#' \code{\link{create_symlinks}}
#' @examples
#' # Load the original paths data frame
#' original_paths <- read.csv("original_paths.csv")
#'
#' # Rename the sample sheet and create symlinks
#' updated_paths <- rename_sample_sheet(original_paths)
#'
#' # View the updated paths data frame
#' updated_paths
#'
#' # View the symlinks in the destination directory
#' list.files("/labs/sbranciamore/projects/rna_nadia_sICKLE/combined_analysis_2023/data/fastq/links")
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr select
rename_sample_sheet <- function(original_sample_sheet) {
    dest_dir <- glue::glue("{analysis_cache}/data/hpapdata/fqlinks")

    # unlink(dest_dir, recursive = TRUE, force = TRUE)
    dir.create(dest_dir, recursive = TRUE)

    # Group the original paths by sample and add new columns for the updated file paths
    df <- original_sample_sheet %>%
        group_by(sample_id) %>%
        mutate(
            Lane_Number = row_number(),
            fastq_1 = file.path(dest_dir, paste0(sample_id, "_S1_L00", Lane_Number, "_R1_001.fastq.gz")),
            fastq_2 = file.path(dest_dir, paste0(sample_id, "_S1_L00", Lane_Number, "_R2_001.fastq.gz")),
        ) |>
        select(sample_id, fastq_1, fastq_2, everything(), -Lane_Number)

    # Function to create symlinks
    create_symlinks <- function(original_df, updated_df) {
        mapply(function(old_1, old_2, new_1, new_2) {
            if (!is.na(old_1) && !is.na(new_1)) {
                file.symlink(old_1, new_1)
            }
            if (!is.na(old_2) && !is.na(new_2)) {
                file.symlink(old_2, new_2)
            }
        }, original_df$fastq_1, original_df$fastq_2, updated_df$fastq_1, updated_df$fastq_2)
    }

    # Create symlinks
    create_symlinks(original_sample_sheet, df)

    # Return the updated data frame
    return(df)
}
