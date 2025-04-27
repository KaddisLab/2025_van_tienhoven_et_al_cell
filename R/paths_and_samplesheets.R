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
#' hpap_fastq_to_sample_sheet(fastq_10x)
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
