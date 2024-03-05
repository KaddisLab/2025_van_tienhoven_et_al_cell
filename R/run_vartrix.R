#' Run VarTrix Pipeline on Gemini
#'
#' Executes the VarTrix pipeline for single-cell genotyping based on given inputs
#' and a specified variant file (VCF). It prepares and submits a batch job to a SLURM cluster if necessary.
#'
#' @param bam_file Path to the Cellranger BAM file.
#' @param cell_barcodes_file Path to the file with cell barcodes to be evaluated.
#' @param fasta_file Genome fasta file.
#' @param vcf_file Called variant file (VCF).
#' @param min_mapq Minimum read mapping quality to consider. Default is 0.
#' @param threads Number of parallel threads to use. Default is 1.
#'
#' @return The path to the VarTrix run directory if an existing run is found, the path to the
#' batch script if a new run is initiated, or a message indicating that a run has already been
#' submitted. The function aims to ensure that each sample is processed once and efficiently utilizes
#' available genomic data for genotyping.
#'
#' @details This function checks for the existence of previous VarTrix runs for the specified
#' sample and either resumes them or initiates a new run by creating a batch submission script.
#' It handles directory creation, script preparation, and job submission, streamlining the
#' process of single-cell genotyping in a high-throughput computing environment. Requires the
#' VarTrix binary to be installed and available in the R PATH set via .Renviron.
#'
#' @examples
#' run_vartrix("path/to/bam_file.bam", "path/to/cell_barcodes.tsv", "path/to/genome.fasta", "path/to/variants.vcf", "path/to/output_directory")
#'
#' @importFrom glue glue
#' @importFrom here here
#' @export
run_vartrix <- function(cellranger_run_folder, region_vcf, mapq = 20, mode = "consensus", additional_args = "") {

    # check that mode is one of consensus, coverage, or alt_fraction
    if (!(mode %in% c("consensus", "coverage", "alt_fraction"))) {
        stop(glue::glue("Invalid mode: {mode}. Must be one of 'consensus', 'coverage', or 'alt_fraction'."))
    }
    sample_id <- basename(cellranger_run_folder)
    bam_path <- glue::glue("{cellranger_run_folder}/outs/possorted_genome_bam.bam")
    barcode_path <- glue::glue("{cellranger_run_folder}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    run_path <- glue::glue("{analysis_cache}/vartrix_out/{sample_id}_{mode}")
    
    if (!file.exists(glue::glue("{run_path}/barcodes.tsv"))) {
        # copy barcodes to run path and unzip
        system(glue::glue("cp {barcode_path} {run_path}/barcodes.tsv.gz && gunzip {run_path}/barcodes.tsv.gz"))
    }

    #TODO  - cluster specific
    ref_genome <- ifelse(hprcc::get_cluster() == "gemini", "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/fasta/genome.fa", "")

    dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    
    # Check if a successful VarTrix run exists
    if (file.exists(glue::glue("{run_path}/out_matrix.mtx"))) {
        print(glue::glue("Found existing VarTrix run, skipping run for {sample_id}..."))
        return(glue::glue("{run_path}/out_matrix.mtx"))
    } else if (!file.exists(glue::glue("{run_path}/run_vartrix.sh"))) {
        # Create a sbatch script
        script_content <- glue::glue("#!/bin/bash
#SBATCH --job-name={sample_id}_{mode}_vartrix
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

cd {run_path}

vartrix_linux --bam '{bam_path}' \\
        --cell-barcodes '{run_path}/barcodes.tsv' \\
        --fasta  {ref_genome} \\
        --out-variants '{run_path}/var_out.tsv' \\
        --ref-matrix '{run_path}/ref_out.mtx' \\
        --mapq {mapq} \\
        --threads 8 \\
        --scoring-method {mode} \\
        --vcf '{region_vcf}' \\
        {additional_args}
")

        cat(script_content, file = glue::glue("{run_path}/run_vartrix.sh"))

        # Submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_vartrix.sh && sleep 0.1"), wait = FALSE)
        return(glue::glue("{run_path}/run_vartrix.sh"))
    } else {
        print(glue::glue("The VarTrix run for \"{sample_id}\" has already been submitted. Check the SLURM job queue or navigate to {run_path} and submit the `run_vartrix.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_vartrix.sh"))
    }
}


#' Filter and Trim Variants in a VCF File
#'
#' @param vcf_path Path to the input VCF file.
#' @param fai_path Path to the FASTA index (.fai) file for chromosome lengths.
#' @param filter A string representing a genomic range to keep in the format "chr:start-end".
#' @param output_vcf_path Path for the output filtered and compressed VCF file.
#'
#' @importFrom VariantAnnotation readVcf writeVcf header 
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges GRanges rowRanges start end
#' @import utils read.table
#' @import stats setNames
#' @import base system shQuote
filter_and_trim_vcf <- function(vcf_path, fai_path, gr_filter = NULL, output_vcf_path) {
  # Read .fai file
  fai <- read.table(fai_path, stringsAsFactors = FALSE, col.names = c("chromosome", "length", "start", "basesPerLine", "bytesPerLine"))
  chrom_lengths <- setNames(fai$length, fai$chromosome)
  
  # Read VCF
  vcf <- VariantAnnotation::readVcf(vcf_path, "hg38")
  GenomeInfoDb::seqlevelsStyle(vcf) <- "NCBI"
  
  # Modify chromosome names in header & vcf
  hdr <- VariantAnnotation::header(vcf)
  contigs <- hdr@header@listData$contig
  rownames(contigs) <- gsub("^([0-9XY]+)$", "chr\\1", rownames(contigs))
  rownames(contigs) <- gsub("^MT$", "chrM", rownames(contigs))
  hdr@header@listData$contig <- contigs
  VariantAnnotation::header(vcf) <- hdr
  GenomeInfoDb::seqlevels(vcf) <- gsub("^([0-9XY]+)$", "chr\\1", GenomeInfoDb::seqlevels(vcf))
  GenomeInfoDb::seqlevels(vcf) <- gsub("^MT$", "chrM", GenomeInfoDb::seqlevels(vcf))
  
  # Trim variants based on chromosome lengths
  variant_positions <- SummarizedExperiment::rowRanges(vcf)
  seqnames_variants <- as.character(GenomeInfoDb::seqnames(variant_positions))

  valid_positions <- seqnames_variants %in% names(chrom_lengths) &
    GenomicRanges::start(variant_positions) <= chrom_lengths[seqnames_variants] &
    GenomicRanges::end(variant_positions) <= chrom_lengths[seqnames_variants]

  trimmed_vcf <- vcf[valid_positions]
  
  # Apply genomic range filter if specified
  if (!is.null(gr_filter)) {
    range <- GenomicRanges::GRanges(gr_filter)
    within_range <- IRanges::overlapsAny(variant_positions, range)
    trimmed_vcf <- trimmed_vcf[within_range]
  }
  
  # Write filtered VCF to file
  VariantAnnotation::writeVcf(trimmed_vcf, gsub(".gz$", "", output_vcf_path))
  
  # Compress using bgzip
  system(paste("bgzip -f", base::shQuote(gsub(".gz$", "", output_vcf_path))), intern = FALSE)
  return(output_vcf_path)
}

