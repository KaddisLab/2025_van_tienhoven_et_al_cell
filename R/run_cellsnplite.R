#' Run CellSNP-lite Pipeline on Gemini
#'
#' Executes the CellSNP-lite pipeline for single-cell genotyping based on given Cell Ranger outputs
#' and a specified region VCF. It prepares and submits a batch job to a SLURM cluster if necessary.
#'
#' @param cellranger_run_folder Path to the Cellranger output directory, which contains
#' the BAM file and barcode list among other outputs.
#' @param region_vcf Path to the VCF file specifying the genomic regions of interest.
#' @param minMAF Minimum minor allele frequency (MAF) threshold for filtering variants.
#' Default is 0.05.
#' @param minCOUNT Minimum read count threshold for including a variant. Default is 10.
#'
#' @return The path to the CellSNP-lite run directory if an existing run is found, the path to the
#' batch script if a new run is initiated, or a message indicating that a run has already been
#' submitted. The function aims to ensure that each sample is processed once and efficiently utilizes
#' available genomic data for genotyping.
#'
#' @details This function checks for the existence of previous CellSNP-lite runs for the specified
#' sample and either resumes them or initiates a new run by creating a batch submission script.
#' It handles directory creation, script preparation, and job submission, streamlining the
#' process of single-cell genotyping in a high-throughput computing environment. Requires the
#' cellsnp-lite binary to be installed and available in the R PATH set via .Renviron.
#'
#' @examples
#' run_cellsnp_lite("path/to/cellranger/output", "path/to/regions.vcf")
#'
#' @importFrom glue glue
#' @importFrom here here
#' @export
run_cellsnp_lite <- function(cellranger_run_folder, region_vcf, minMAF = 0.05, minCOUNT = 10) {
    
    sample_id <- basename(cellranger_run_folder)
    bam_path <- glue::glue("{cellranger_run_folder}/outs/possorted_genome_bam.bam")
    barcode_path <- glue::glue("{cellranger_run_folder}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")

    run_path <- glue::glue("{analysis_cache}/cellsnplite_out/{sample_id}")

    if (!dir.exists(run_path)) {
        dir.create(run_path, recursive = TRUE)
    }

    # if a successful cellsnp-lite run exists, use that
    if (file.exists(glue::glue("{run_path}/cellSNP.base.vcf"))) {
        print(glue::glue("Found existing cellsnp-lite, skipping run for {sample_id}..."))
        return(glue::glue("{run_path}/cellSNP.base.vcf"))
    # if theres no run script, make one and submit it to the cluster
    } else if (!file.exists(glue::glue("{run_path}/run_cellsnp-lite.sh"))) {
    # make a sbatch script
        cat(glue::glue("#!/bin/bash
#SBATCH --job-name={sample_id}_cellsnplite
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=20G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

source /etc/sysconfig/modules/init.sh
PATH=~/bin:~/miniconda3/bin:$PATH

cd {run_path}
cellsnp-lite -s '{bam_path}' \\
    -b '{barcode_path}' \\
    -O '{run_path}' \\
    -R '{region_vcf}' \\
    -p 6 \\
    --genotype \\
    --minMAF {minMAF} \\
    --minCOUNT {minCOUNT}
", .trim = FALSE),
            file = glue::glue("{run_path}/run_cellsnp-lite.sh")
        )
        # submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_cellsnp-lite.sh && sleep 0.1"), wait = FALSE)
        return(NULL)
    } else {
        # If there's no cellsnp output but there's a run script, just return the run script path
        # with a message about what to do next
        print(glue::glue("The cellsnp-lite run for \"{sample_id}\" has already been submitted.
Check the SLURM job queue or navigate to
{run_path}
and submit the `run_cellsnp-lite.sh` script to resume the run."))
        return(NULL)
    }
}

#' Get CellSNP-lite Genotypes
#'
#' Searches through `.vcf` files within the `cellsnplite_out` directory to find genotypes
#' matching a specific pattern and locus. Filters results based on minimum minor allele frequency (MAF)
#' and read count.
#'
#' @param pattern A pattern to search within `.vcf` files, typically involving chromosome
#' and position information.
#' @param locus The specific locus identifier, usually an rsID, for which genotypes are being queried.
#' @param minMAF Minimum minor allele frequency (MAF) threshold for filtering genotypes. Default is 0.05.
#' @param minCOUNT Minimum read count threshold for filtering genotypes. Default is 10.
#'
#' @return A tibble with columns `sample_id`, `ad` (allele depth), `dp` (depth of coverage),
#' `oth` (other reads count), `maf` (minor allele frequency), and `locus`. The tibble is
#' filtered by the specified `minMAF` and `minCOUNT` thresholds and sorted by the numeric part
#' of the `sample_id`.
#'
#' @examples
#' # Assuming a pattern for chromosome 11 and position 2159830
#' get_cellsnp_lite_genotypes("11\t2159830", "rs3842753")
#'
#' @export
#'
#' @importFrom dplyr filter arrange mutate select
#' @importFrom tidyr separate
#' @importFrom stringr str_extract
#' @importFrom glue glue
#' @importFrom tibble tibble
get_cellsnp_lite_genotypes <- function(pattern, locus, minMAF = 0.05, minCOUNT = 1000) {
        system(glue::glue("find {analysis_cache}/cellsnplite_out -type f -name '*.base.vcf' -exec grep -Hn '{pattern}' {{}} +"), intern = TRUE) %>%
        tibble(vcf_output = .) %>%
        tidyr::separate(vcf_output, into = c("File", "MatchedLine"), sep = ":", extra = "merge") %>%
        mutate(
            sample_id = stringr::str_extract(File, "HPAP-\\d+"),
            ref = stringr::str_extract(MatchedLine, "(?<=\\t)[AGCT]+(?=\\t)"),
            alt = stringr::str_extract(MatchedLine, "(?<=\\t[AGCT]\\t)[AGCT]+"),
            ad = as.numeric(stringr::str_extract(MatchedLine, "(?<=AD=)\\d+")),
            dp = as.numeric(stringr::str_extract(MatchedLine, "(?<=DP=)\\d+")),
            maf = round(ad/dp, 2),
            oth = as.numeric(stringr::str_extract(MatchedLine, "(?<=OTH=)\\d+")),
            gt = if_else(maf >= 0.7, paste0(alt, alt),  if_else(maf >= 0.3, paste0(ref, alt), paste0(ref, ref))),
            locus = locus) |>
        select(sample_id, ad, dp, oth, maf, ref, alt, gt, locus) |>
        dplyr::filter(maf >= minMAF, dp >= minCOUNT) |>
        mutate(numeric_part = as.numeric(stringr::str_extract(sample_id, "\\d+"))) %>%
        arrange(numeric_part) %>%
        select(-numeric_part)
}

