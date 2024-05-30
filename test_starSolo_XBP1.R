#tar_source()
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(dplyr)
library(readr)


## Get gene data
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

# Define the span of XBP1
xbp1_span <- range(xbp1_exons)


# Adjust the path to your SJ.out.tab file
sj_file <- "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/starsolo_out/HPAP-047/Solo.out/SJ/raw/features.tsv"

# Read the SJ.out.tab file
sj_data <- read_tsv(
    sj_file,
    col_names = c("chromosome", "start", "end", "strand", "intron_motif", "annotated", "unique_reads", "multi_reads", "max_overhang"),
    show_col_types = FALSE
)

# Convert the junction data frame to a GRanges object
junctions <- GRanges(
    seqnames = as.factor(sj_data$chromosome),
    ranges = IRanges(start = sj_data$start, end = sj_data$end),
    strand = ifelse(sj_data$strand == 1, "+", "-")
)

# Define the critical splicing junction for XBP1s
critical_junction <- GRanges(
    seqnames = "chr22",
    ranges = IRanges(start = 28796122, end = 28796147),
    strand = "-"
)

# Find overlaps with the critical splicing junction
critical_overlaps <- findOverlaps(junctions, critical_junction, type = "equal")

# Count reads for the critical splicing junction (spliced isoform XBP1s)
spliced_counts <- sum(sj_data$unique_reads[queryHits(critical_overlaps)])

# Filter SJ data for XBP1 gene span
xbp1_sj <- sj_data %>%
    filter(chromosome == as.character(seqnames(xbp1_span)) & start >= start(xbp1_span) & end <= end(xbp1_span))

# Calculate the total reads supporting inclusion and exclusion
total_reads <- sum(xbp1_sj$unique_reads)

# Calculate PSI for the critical junction
psi <- (spliced_counts / total_reads) * 100
print(paste("PSI for XBP1 critical junction:", psi))
