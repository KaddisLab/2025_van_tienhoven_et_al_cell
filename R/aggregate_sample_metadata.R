#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param pancdb_metadata
#' @param nameme1
#' @return
#' @author Denis O'Meally
#' @export
aggregate_sample_metadata <- function(pancdb_metadata,
                                    protected_cohort,
                                    rs3842753_cohort,
                                    rs689_cohort,
                                    xbp1_psi_per_sample) {

# # Sample metadata -------------------------------------------------------------------------------
pancdb_metadata$protected <- pancdb_metadata$donor_id %in% protected_cohort$sample_id

# Function to perform the merge and add the genotype column
merge_cohort <- function(metadata, cohort, new_col_name) {
    genotype_key <- cohort %>%
        select(sample_id, ref, gt)

    metadata <- metadata %>%
        left_join(genotype_key, by = c("sample_id" = "sample_id")) %>%
        mutate(!!new_col_name := if_else(is.na(gt), paste0(genotype_key$ref[1], genotype_key$ref[1]), gt)) %>%
        select(-gt, -ref)

    return(metadata)
}

# Apply the function to each cohort
pancdb_metadata <- merge_cohort(pancdb_metadata, protected_cohort, "rs3842752_consensus")
pancdb_metadata <- merge_cohort(pancdb_metadata, rs3842753_cohort, "rs3842753_consensus")
pancdb_metadata <- merge_cohort(pancdb_metadata, rs689_cohort, "rs689_consensus")

# Add the XBP1 PSI per-sample values
pancdb_metadata <- pancdb_metadata %>%
    left_join(xbp1_psi_per_sample, by = c("sample_id" = "sample_id"))

return(pancdb_metadata)

}
