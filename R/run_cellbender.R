
#' Run CellBender analysis on a Cell Ranger run folder
#'
#' This function runs the CellBender analysis on a Cell Ranger run folder.
#'
#' @param cellranger_run_folder The path to the Cell Ranger run folder.
#' @param additional_args Additional arguments to be passed to the CellBender command.
#'
#' @return The path to the CellBender output file.
#'
#' @examples
#' if_interactive({
#' run_cellbender("/path/to/cellranger_run_folder")
#' run_cellbender("/path/to/cellranger_run_folder", additional_args = "--param1 value1 --param2 value2")
#' })
#' 
#' @export
run_cellbender <- function(cellranger_run_folder, additional_args = "") {

    sample_id <- basename(cellranger_run_folder)
    run_path <- glue::glue("{analysis_cache}/cellbender_out/{sample_id}")

    raw_feature_matrix <- glue::glue("{cellranger_run_folder}/outs/raw_feature_bc_matrix.h5")
    out_file_h5 <- glue::glue("{run_path}/{sample_id}_cellbender.h5")

    cellbender_exe <- "/packages/singularity/3.11.5/bin/singularity exec --nv $SINGULARITY_CACHEDIR/cellbender_latest.sif cellbender"

    dir.create(run_path, recursive = TRUE, showWarnings = FALSE)
    
    # Check if a successful cellbender run exists
    if (file.exists(glue::glue("{run_path}/{sample_id}_cellbender.h5"))) {
        print(glue::glue("Found existing cellbender run, skipping run for {sample_id}..."))
        return(glue::glue("{run_path}/{sample_id}_cellbender.h5"))
    } else if (!file.exists(glue::glue("{run_path}/run_cellbender.sh"))) {
        # Create a sbatch script
        script_content <- glue::glue("#!/bin/bash
#SBATCH --job-name={sample_id}_cellbender
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --chdir={run_path}
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

{cellbender_exe} remove-background --cuda --cpu-threads 8 --input {raw_feature_matrix} --output {out_file_h5} \\
        {additional_args}
")
        cat(script_content, file = glue::glue("{run_path}/run_cellbender.sh"))
        # Submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_cellbender.sh && sleep 0.1"), wait = FALSE)
        return(glue::glue("{run_path}/run_cellbender.sh"))
    } else {
        print(glue::glue("The CellBender run for \"{sample_id}\" has already been submitted. Check the SLURM job queue or navigate to {run_path} and submit the `run_cellbender.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_cellbender.sh"))
    }
}
