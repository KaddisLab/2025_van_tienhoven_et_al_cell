#' Collate and Merge Sample Metadata with Genotype Information
#'
#' Combines sample metadata with genotype information from multiple cohorts
#' and adds XBP1 splicing data. Also marks protected samples based on donor IDs.
#'
#' @param pancdb_metadata Data frame containing PANC-db sample metadata
#' @param protected_cohort Data frame containing protected cohort information with sample_id and genotype data
#' @param rs3842753_cohort Data frame containing rs3842753 cohort information with sample_id and genotype data
#' @param rs689_cohort Data frame containing rs689 cohort information with sample_id and genotype data
#' @param xbp1_psi_per_sample Data frame containing XBP1 PSI values per sample
#'
#' @return A data frame containing the merged metadata with genotype information and XBP1 PSI values
#'
#' @author Denis O'Meally
#' @export
collate_sample_metadata <- function(pancdb_metadata,
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
