#' Get Cell-level Genotypes from CellSNP-lite Output
#'
#' Processes CellSNP-lite output files to extract allele-specific expression data
#' for specific SNPs of interest (rs3842752, rs3842753, rs689) at the single-cell level.
#'
#' @param vcf_paths Character vector of paths to CellSNP-lite base VCF files.
#'        These paths are used to locate the corresponding AD, DP, and OTH matrices.
#'
#' @return A tibble with cell-level genotype information containing:
#'   - cell: Cell barcode with sample ID prefix
#'   - [snp]_AD: Allele depth for each SNP
#'   - [snp]_OTH: Other allele counts for each SNP
#'   - [snp]_DP: Total read depth for each SNP
#'   where [snp] is one of rs3842752, rs3842753, or rs689
#'
#' @import Matrix
#' @importFrom tidyverse map_dfr
#' @importFrom purrr map_dfr
#' @importFrom progress progress_bar
#' @export
get_cellsnp_lite_cell_genotypes <- function(vcf_paths) {
    require(Matrix)
    require(tidyverse)
    require(purrr)
    require(progress)

    snps <- c("rs3842752" = "11:2159843", "rs3842753" = "11:2159830", "rs689" = "11:2160994")

    process_sample <- function(base_vcf_path) {
        folder_path <- dirname(base_vcf_path)
        sample_name <- basename(folder_path)

        # Read and transpose matrices
        ad_mat <- t(readMM(file.path(folder_path, "cellSNP.tag.AD.mtx")))
        dp_mat <- t(readMM(file.path(folder_path, "cellSNP.tag.DP.mtx")))
        oth_mat <- t(readMM(file.path(folder_path, "cellSNP.tag.OTH.mtx")))

        # Get SNP positions
        vcf_data <- read.table(base_vcf_path, comment.char = "#", stringsAsFactors = FALSE)
        snp_positions <- setNames(seq_len(nrow(vcf_data)), paste(vcf_data$V1, vcf_data$V2, sep = ":"))

        # Get cell barcodes and create result tibble
        cell_barcodes <- readLines(file.path(folder_path, "cellSNP.samples.tsv"))
        result <- tibble(cell = paste0(sample_name, "_", cell_barcodes))

        for (snp in names(snps)) {
            position <- snps[snp]
            idx <- snp_positions[position]
            if (!is.na(idx) && idx <= ncol(ad_mat)) {
                result[[paste0(snp, "_AD")]] <- ad_mat[, idx]
                result[[paste0(snp, "_OTH")]] <- oth_mat[, idx]
                result[[paste0(snp, "_DP")]] <- dp_mat[, idx]
            } else {
                result[[paste0(snp, "_AD")]] <- 0
                result[[paste0(snp, "_OTH")]] <- 0
                result[[paste0(snp, "_DP")]] <- 0
            }
        }

        return(result)
    }

    # Create progress bar
    pb <- progress_bar$new(
        format = "Processing sample [:bar] :current/:total (:percent) :eta",
        total = length(vcf_paths),
        clear = FALSE,
        width = 60
    )

    # Process all samples with progress bar
    all_results <- map_dfr(vcf_paths, function(path) {
        result <- process_sample(path)
        pb$tick()
        result
    })

    return(all_results)
}

# # Usage example
# result_tibble <- parse_cellsnp_outputs(cellsnp_lite)
# print(head(result_tibble))
# print(paste("Total rows:", nrow(result_tibble)))
# print(paste("Total samples:", length(unique(sub("_.*", "", result_tibble$cell)))))
