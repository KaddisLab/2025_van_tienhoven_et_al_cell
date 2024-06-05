
calculate_xbp1_psi_per_cell <- function(path) {
    
require(dplyr)
require(Matrix)
require(readr)
require(GenomicRanges)
require(IRanges)
require(S4Vectors)
require(stringr)
require(qs)
require(tibble)
require(glue)
require(rtracklayer)

    orig.ident <- stringr::str_extract(path, "HPAP-\\d+")

    message("Start calculating XBP1 PSI per cell for sample ", orig.ident)

    # Load the GTF file to get exon information for XBP1
    if (file.exists(glue::glue("{analysis_cache}/data/exons/xbp1_exons.qs"))) {
        xbp1_exons <- qs::qread(glue::glue("{analysis_cache}/data/exons/xbp1_exons.qs"))
    } else {
        gtf_path <- "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
        all_annotations <- import(gtf_path)

        xbp1_exons <- subset(all_annotations, gene_name == "XBP1" & type == "exon")
        seqlevels(xbp1_exons) <- seqlevelsInUse(xbp1_exons)
        qs::qsave(xbp1_exons, glue::glue("{analysis_cache}/data/exons/xbp1_exons.qs"))
    }

    xbp1_span <- GenomicRanges::reduce(xbp1_exons)
    message("Loaded XBP1 exons")

    # Correct the paths
    matrix_path <- gsub("Summary.csv", "raw/matrix.mtx", path, fixed = TRUE)
    barcodes_path <- gsub("Summary.csv", "raw/barcodes.tsv", path, fixed = TRUE)
    features_path <- gsub("Summary.csv", "raw/features.tsv", path, fixed = TRUE)

    spliced <- as(Matrix::readMM(matrix_path), "CsparseMatrix") # Convert explicitly to dgCMatrix

    barcodes <- readr::read_lines(barcodes_path)
    features <- readr::read_tsv(features_path, col_names = c("chromosome", "start", "end", "type"), show_col_types = FALSE, progress = FALSE)
    message("Loaded spliced junction data")

    # Convert features to a GRanges object
    features_gr <- GRanges(
        seqnames = features$chromosome,
        ranges = IRanges(start = features$start, end = features$end),
        strand = ifelse(features$type == 1, "+", "-")
    )

    # Define the critical splicing junction for XBP1s
    critical_junction <- GRanges(
        seqnames = "chr22",
        ranges = IRanges(start = 28796122, end = 28796147),
        strand = "-"
    )

    # Find overlaps with XBP1 span
    xbp1_overlaps <- findOverlaps(features_gr, xbp1_span)
    xbp1_indices <- queryHits(xbp1_overlaps)

    # Find overlap for the critical junction
    critical_overlaps <- findOverlaps(features_gr, critical_junction, type = "equal")
    critical_index <- queryHits(critical_overlaps)

    if (length(critical_index) == 0) {
        stop("Critical junction not found in the features.")
    }

    # Calculate the total reads supporting inclusion and exclusion for each cell
    xbp1_counts <- spliced[xbp1_indices, ]
    total_reads_per_cell <- colSums(xbp1_counts)

    # Calculate reads for the critical splicing junction (spliced isoform XBP1s)
    spliced_counts_per_cell <- spliced[critical_index, ]

    # Calculate PSI for the critical junction per cell
    psi_per_cell <- (spliced_counts_per_cell / total_reads_per_cell) * 100
    message("Calculated PSI per cell")

    # Create tibble with cell and PSI values
    psi_tibble <- tibble(cell = barcodes, xbp1_psi = psi_per_cell, xbp1_u = spliced_counts_per_cell, xpb1_tot = total_reads_per_cell) |>
        dplyr::filter(!is.na(xbp1_psi)) |>
        dplyr::mutate(cell = glue::glue("{orig.ident}_{cell}-1"))

    psi_tibble_path <- glue::glue("{analysis_cache}/xbp1_psi_out/{orig.ident}_xbp1_psi_per_cell.csv")
    dir.create(dirname(psi_tibble_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(psi_tibble, psi_tibble_path, row.names = FALSE, quote = FALSE)

    return(psi_tibble_path)
}

# # Example usage with the provided path
# path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-122/Solo.out/SJ/Summary.csv"
# result_tibble <- calculate_xbp1_psi_per_cell(path)
# print(result_tibble)
