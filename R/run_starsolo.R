#' Run STARsolo Analysis on CellRanger BAM Files
#'
#' Processes CellRanger BAM files using STARsolo to generate splice junction and
#' RNA velocity information. Supports 10X Genomics v2 and v3 chemistry.
#'
#' @param cellranger_folder Path to CellRanger output folder containing possorted_genome_bam.bam
#' @param technology String indicating chemistry version ("10XV2" or "10XV3")
#' @param ref_dir Path to STAR genome reference directory
#' @param out Output directory (default: "{analysis_cache}/starsolo_out/{sample_id}")
#'
#' @return Path to the generated Summary.csv file, or NULL if processing failed
#'
#' @details
#' The function:
#' 1. Checks for existing output to avoid reprocessing
#' 2. Configures STAR parameters based on chemistry version
#' 3. Runs STAR with RNA velocity and splice junction analysis
#' 4. Uses available SLURM resources if running on cluster
#'
#' Output includes Gene counts, GeneFull lengths, splice junctions (SJ), and
#' RNA velocity information in the Solo.out directory.
#'
#' @importFrom glue glue
#' @importFrom stringr str_detect
#' @export
run_starsolo <- function(cellranger_folder, technology, ref_dir, out = NULL) {
    
    sample_id <- basename(cellranger_folder)
    # build the output path
    if (is.null(out)) out <- glue::glue("{analysis_cache}/starsolo_out/{sample_id}")

    # Check for an existing run
    return_path <- glue::glue("{out}/Solo.out/SJ/Summary.csv")
    if (file.exists(return_path)) {
        return(return_path)
    }

    resources <- hprcc::slurm_allocation()
    if (!is.null(resources)) {
        cpus <- resources$CPUs
        job_id <- resources$JobID
    } else {
        cpus <- 2
    }
    bam_file <- glue::glue("{cellranger_folder}/outs/possorted_genome_bam.bam")
    if (stringr::str_detect(technology, "10XV2")) whitelist <- "/ref_genomes/cellranger/human/cellranger_whitelists/10x_version2_whitelist.txt"
    if (stringr::str_detect(technology, "10XV3")) whitelist <- "/ref_genomes/cellranger/human/cellranger_whitelists/10x_version3_whitelist.txt"

    dir.create(out, showWarnings = FALSE, recursive = TRUE)
        
    # Construct the command with dynamic parameters
    command <- glue::glue("cd {out};$HOME/bin/STAR --runThreadN {cpus} \\
     --genomeDir {ref_dir} \\
     --soloType CB_UMI_Simple \\
     --readFilesIn {bam_file} \\
     --readFilesCommand 'samtools view -F 0x100' \\
     --readFilesType SAM SE \\
     --soloInputSAMattrBarcodeSeq CR UR \\
     --soloInputSAMattrBarcodeQual CY UY \\
     --soloCBwhitelist {whitelist} \\
     --soloFeatures Gene GeneFull SJ Velocyto \\
     --soloBarcodeReadLength 0 \\
     --outFileNamePrefix {out}")

    # Execute the command using system()
    result <- system(command, intern = TRUE)
    Log_message(result)
    
    if (file.exists(return_path)) return(return_path) else return(NULL)
}
