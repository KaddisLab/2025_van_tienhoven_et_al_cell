#' Calculate Percent Spliced-In (PSI) for Multiple Samples
#'
#' This function calculates the Percent Spliced-In (PSI) for the XBP1 gene
#' across multiple samples using STARsolo output files. The function reads
#' the splice junction data, identifies the critical splicing junction
#' for XBP1, and computes the PSI for each sample.
#'
#' @param paths A character vector of file paths to the STARsolo output 
#'              summary files. The function expects the splice junction
#'              data to be in a subdirectory `raw/features.tsv`.
#'
#' @return A data frame with two columns: `orig.ident` (the sample identifier)
#'         and `xbp1_psi` (the calculated PSI for the XBP1 critical junction).
#'
#' @examples
#' \dontrun{
# ' paths <- c(
# '   "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-047//Solo.out/SJ/Summary.csv",
# '   "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-020//Solo.out/SJ/Summary.csv",
# '   "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-021//Solo.out/SJ/Summary.csv"
# ' )
#' result_table <- calculate_psi(paths)
#' print(result_table)
#' }
#'
#' @importFrom qs qread
#' @importFrom readr read_tsv
#' @importFrom dplyr filter
#' @importFrom GenomicRanges GRanges reduce seqnames
#' @importFrom IRanges IRanges start end
#' @importFrom S4Vectors queryHits
#' @importFrom stringr str_extract
#' @export
calculate_xbp1_psi_per_sample <- function(paths) {
require(GenomicRanges)
require(dplyr)
require(readr)
require(stringr)
require(qs)
require(glue)
require(rtracklayer)

# Load the GTF file to get exon information for XBP1
if (file.exists(glue::glue("{analysis_cache}/data/xbp1_exons.qs"))) {
    # Load the exon information from a serialized file
    xbp1_exons <- qs::qread(glue::glue("{analysis_cache}/data/xbp1_exons.qs"))
} else {
    # Load the GTF file to get exon information for XBP1
    gtf_path <- "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    all_annotations <- import(gtf_path)

    # Filter for XBP1 exons
    xbp1_exons <- subset(all_annotations, gene_name == "XBP1" & type == "exon")
    seqlevels(xbp1_exons) <- seqlevelsInUse(xbp1_exons)
    glue::glue("{analysis_cache}/data/xbp1_exons.qs")
}

    xbp1_span <- range(xbp1_exons)

    # Prepare a data frame to store results
    results <- data.frame(orig.ident = character(), xbp1_psi = numeric(), stringsAsFactors = FALSE)

    # Iterate over each path in the input vector
    for (path in paths) {
        # Extract the identifier from the path
        orig.ident <- stringr::str_extract(path, "HPAP-\\d+")

        # Adjust the path to point to the raw/features.tsv file
        sj_file <- gsub("Summary.csv", "raw/features.tsv", path, fixed = TRUE)

        # Read the SJ.out.tab file
        sj_data <- readr::read_tsv(
            sj_file,
            col_names = c("chromosome", "start", "end", "strand", "intron_motif", "annotated", "unique_reads", "multi_reads", "max_overhang"),
            show_col_types = FALSE
        )

        # Convert the junction data frame to a GRanges object
        junctions <- GenomicRanges::GRanges(
            seqnames = as.factor(sj_data$chromosome),
            ranges = IRanges::IRanges(start = sj_data$start, end = sj_data$end),
            strand = ifelse(sj_data$strand == 1, "+", "-")
        )

        # Define and find overlaps with the critical splicing junction for XBP1s
        critical_junction <- GenomicRanges::GRanges(
            seqnames = "chr22",
            ranges = IRanges::IRanges(start = 28796122, end = 28796147),
            strand = "-"
        )
        critical_overlaps <- GenomicRanges::findOverlaps(junctions, critical_junction, type = "equal")

        # Count reads for the critical splicing junction (spliced isoform XBP1s)
        spliced_counts <- sum(sj_data$unique_reads[S4Vectors::queryHits(critical_overlaps)])

        # Filter SJ data for XBP1 gene span
        xbp1_sj <- sj_data %>%
            filter(
                chromosome == as.character(GenomicRanges::seqnames(xbp1_span)),
                start >= IRanges::start(xbp1_span),
                end <= IRanges::end(xbp1_span)
            )

        # Calculate the total reads supporting inclusion and exclusion
        total_reads <- sum(xbp1_sj$unique_reads)

        # Calculate PSI for the critical junction
        psi <- if (total_reads > 0) {
            (spliced_counts / total_reads) * 100
        } else {
            NA
        }

        # Append results
        results <- rbind(results, data.frame(sample_id = orig.ident, xbp1u_psi = psi))
    }

    return(results)
}

