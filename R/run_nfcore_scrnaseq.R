#' Run nf-core/scrnaseq pipeline on Gemini
#'
#' Makes an `sbatch` job script for the nf-core/scrnaseq pipeline and submits it to SLURM.
#' The pipeline version is specified by the `scrnaseq_release` parameter in the `_targets.R` file. The [default
#' parameters](https://nf-co.re/scrnaseq/latest/usage) are used for the pipeline.
#'
#' Execution of the pipeline is independent of the targets script, allowing asynchronous development.
#' If no `run_folder` exists, it will be created along with a `sample_sheet.csv` and `run_{ref_genome}.sh`,
#' which is submitted via `sbatch`. A successful run is detected by the presence of a a multiqc report.
#' An existing run script wont be overwritten by this function so it must be deleted manually.
#'
#' @name run_nf_core_scrnaseq
#' @param run_folder subfolder in which to run nextflow, keeps the pipeline
#' from overwriting previous runs and enables the `-resume` feature of Nextflow.
#' @param sample_sheet nf-core sample sheet with the columns `library_id, fastq_1, fastq_2, strandedness` (see [`R/import_metadata.R`](https://github.com/drejom/haemdata/blob/HEAD/R/import_metadata.R))
#' @param protocol sequencing protocol, either `10XV1`, `10XV2` or `10XV3`
#' @return a path to the run script, or the multiqc report if it exists.
#' @author Denis O'Meally
#' @export

run_nfcore_scrnaseq <- function(run_folder, sample_sheet, protocol = "10XV3") {
    ref_genome <- "refdata-gex-GRCh38-2020-A"
    run_path <- glue::glue("{analysis_cache}/{run_folder}")
    out_folder <- glue::glue("nfcore-scrnaseq-v{scrnaseq_release}-{ref_genome}")

    if (!dir.exists(run_path)) {
        dir.create(run_path, recursive = TRUE)
    }

        # if a full-run multiqc report exists, use that
    if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/scrnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))
        # if a STAR-only multiqc report exists, use that
    } else if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/scrnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))
        # if theres no run script, make one and submit it to the cluster
    } else if (!file.exists(glue::glue("{run_path}/run_{ref_genome}.sh"))) {
        sample_sheet |>
            dplyr::rename(sample = sample_id) |>
            readr::write_csv(paste0(run_path, "/sample_sheet.csv"), quote = "all")

        # make a sbatch script
        cat(glue::glue("#!/bin/bash
#SBATCH --job-name={run_folder}_{ref_genome}
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

source /etc/sysconfig/modules/init.sh
PATH=~/bin:$PATH

cd {run_path}
module load singularity
module load Java
export NXF_ANSI_LOG=FALSE

nextflow run \\
    -profile singularity \\
    nf-core/scrnaseq -r {scrnaseq_release} -resume \\
    --publish_dir_mode copy \\
    --input {run_path}/sample_sheet.csv \\
    --outdir {out_folder} \\
    --cellranger_index /ref_genomes/cellranger/human/{ref_genome} \\
    --gtf /ref_genomes/cellranger/human/{ref_genome}/genes/genes.gtf \\
    --fasta /ref_genomes/cellranger/human/{ref_genome}/fasta/genome.fa \\
    --email domeally@coh.org --aligner cellranger --protocol {protocol}
", .trim = FALSE),
            file = glue::glue("{run_path}/run_{ref_genome}.sh")
        )
        # submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_{ref_genome}.sh"), wait = FALSE)
        return(glue::glue("{run_path}/run_{ref_genome}.sh"))
    } else {
        # If there's no multiqc report but there's a run script, just return the run script path
        # with a message about what to do next
        print(glue::glue("The nextflow run for \"{run_folder}\" has already been submitted.
Check the SLURM job queue or navigate to
{run_path}
and submit the `run_{ref_genome}.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_{ref_genome}.sh"))
    }
}
