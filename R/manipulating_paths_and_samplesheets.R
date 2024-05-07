#' Process HPAP 10x Genomics Fastq File Paths Sample Sheet
#'
#' Takes a vector of 10x Genomics fastq file paths and returns a dataframe
#' with extracted sample IDs, R1 and R2 reads.
#'
#' @param fastq_10x A character vector of fastq file paths.
#'
#' @return A dataframe with columns for sample ID, R1 and R2 file paths.
#' @export
#'
#' @examples
#' fastq_10x <- c("path/to/HPAP-001_R1_file.fastq.gz", "path/to/HPAP-001_R2_file.fastq.gz")
#' process_10x_samples(fastq_10x)
hpap_fastq_to_sample_sheet <- function(fastq_10x) {
  data.frame(file_path = fastq_10x) %>%
    dplyr::mutate(
      sample_id = stringr::str_extract(file_path, "HPAP-\\d{3}"),
      file_base = stringr::str_replace(file_path, "(R[12]_).*.gz$", ""),
      read_type = dplyr::if_else(stringr::str_detect(file_path, "R1_"), "R1", "R2")
    ) %>%
    dplyr::group_by(sample_id, file_base) %>%
    dplyr::summarise(
      R1 = list(file_path[read_type == "R1"]),
      R2 = list(file_path[read_type == "R2"]),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      R1 = purrr::map(R1, ~if(length(.x) == 0) NA else .x),
      R2 = purrr::map(R2, ~if(length(.x) == 0) NA else .x)
    ) |> 
    tidyr::unnest(c(R1, R2)) |>
    dplyr::select(sample_id, fastq_1 = R1, fastq_2 = R2)
}

#' Write nf-core Sample Sheet
#'
#' This function generates a sample sheet for NFCore pipelines based on input 10x Genomics FASTQ data. It allows for specifying the chemistry version, modifies the sample sheet to include strandedness information, writes it to a CSV file within a data directory, and returns the file path.
#'
#' @param fastq_10x A data frame or list containing paths to FASTQ files or relevant metadata for 10x Genomics sequencing data.
#' @param chemistry The version of the 10x Genomics chemistry used for sequencing. Defaults to "10xv3".
#'
#' @return
write_nfcore_sample_sheet <- function(fastq_10x, chemistry = "10xv3") {

# Define the file path
file_path <- glue::glue("{analysis_cache}/data/sample_sheet_{chemistry}.csv")

# Write the dataframe to the file
readr::write_csv(hpap_fastq_to_sample_sheet(fastq_10x) |> 
  mutate(strandedness = "reverse"), file_path)

# Return the file path
file_path

}

#' Process 10x Genomics Fastq File Paths for kallisto-bustools
#'
#' Takes a vector of 10x Genomics fastq file paths and returns a dataframe
#' with extracted sample IDs, file base names, and separated R1 and R2 reads.
#'
#' @param fastq_10x A character vector of fastq file paths.
#'
#' @return A dataframe with columns for sample ID, file base, and separated lists of R1 and R2 file paths.
#' @export
#'
#' @examples
#' fastq_10x <- c("path/to/HPAP-001_R1_file.fastq.gz", "path/to/HPAP-001_R2_file.fastq.gz")
#' process_10x_samples(fastq_10x)
arrange_by_sample_10x <- function(fastq_10x) {
    fastq_10x |>
    hpap_fastq_to_sample_sheet() |>
    concat_fastq_paths()
}


#' Process SureSeq Fastq File Paths
#'
#' Takes a vector of SureSeq fastq file paths and returns a dataframe
#' with extracted sample IDs and file paths. This function is designed
#' for handling single-end libraries.
#'
#' @param fastq_paths A character vector of fastq file paths.
#'
#' @return A dataframe with columns for sample ID and file path.
#' @export
#'
#' @examples
#' fastq_paths <- c("/path/to/HPAP-051_scRNA_SSq2-073A01-R1_fastq-data.fastq.gz")
#' process_sureseq_samples(fastq_paths)
arrange_by_sample_ssq <- function(fastq_paths) {
  data.frame(file_path = fastq_paths) %>%
    dplyr::mutate(
      sample_id = stringr::str_extract(file_path, "HPAP-\\d{3}"),
      file_base = stringr::str_replace(file_path, "_R[12]_fastq-data.fastq.gz$", "")
    ) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      fastq_1 = list(file_path),
      .groups = 'drop'
    ) %>%
    tidyr::unnest(c(fastq_1)) |>
    concat_fastq_paths()
}


#' Concatenate Fastq Paths
#'
#' Takes a dataframe containing FASTQ file paths grouped by sample ID and
#' concatenates the paths into a single character string for each sample.
#' Handles datasets with or without R2 reads.
#'
#' @param df A dataframe with columns `sample_id`, `R1`, and optionally `R2`.
#' @return A dataframe with `sample_id` and `fastq_paths`, where `fastq_paths`
#'         contains concatenated FASTQ file paths, each surrounded by single
#'         quotes and separated by spaces. If `R2` paths are present, they are
#'         interleaved with `R1` paths.
#' @export
#' @examples
#' df <- data.frame(
#'   sample_id = c("Sample1", "Sample1", "Sample2"),
#'   R1 = c("path/to/Sample1_R1_file1.fastq.gz", "path/to/Sample1_R1_file2.fastq.gz", "path/to/Sample2_R1_file1.fastq.gz"),
#'   R2 = c("path/to/Sample1_R2_file1.fastq.gz", "path/to/Sample1_R2_file2.fastq.gz", NA) # Optional
#' )
#' concat_fastq_paths(df)
concat_fastq_paths <- function(df) {
  # Ensure 'R2' column exists for consistent processing
  if(!"fastq_2" %in% names(df)) {
    df$fastq_2 <- NA
  }
  
  df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      fastq_paths = purrr::map2_chr(fastq_1, fastq_2, ~{
        r1 <- .x
        r2 <- ifelse(is.na(.y), "", .y)
        paths <- c(r1, r2)
        paths <- paths[paths != ""]
        paste0(paths, collapse = " ")
      }) %>%
        purrr::reduce(paste, collapse = " "), # Ensure a single string per group
      .groups = 'drop'
    )
}
