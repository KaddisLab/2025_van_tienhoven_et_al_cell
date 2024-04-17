#' Build a kallisto index and transcript-to-gene mapping with kb ref
#'
#' This function acts as an R wrapper for the `kb ref` command-line utility, facilitating the 
#' building of a kallisto index and the generation of a transcript-to-gene (T2G) mapping. It 
#' supports the full range of `kb ref` options, including handling of various workflows and 
#' species-specific index building.
#'
#' @param index Path to the kallisto index to be constructed. If `n` is also specified, this is the prefix for the n indices to construct.
#' @param t2g Path to transcript-to-gene mapping to be generated.
#' @param fasta Genomic FASTA file(s), comma-delimited.
#' @param gtf Reference GTF file(s), comma-delimited.
#' @param feature [Optional; `kite` workflow only] Path to TSV containing barcodes and feature names.
#' @param tmp Override default temporary directory.
#' @param keep_tmp Do not delete the tmp directory post-operation.
#' @param verbose Print debugging information.
#' @param n Number of files to split the index into for parallel processing.
#' @param species Download a pre-built kallisto index for `human`, `mouse`, or `linnarsson` instead of building it locally.
#' @param k Override the k-mer length of the index.
#' @param workflow Type of workflow (`standard`, `lamanno`, `nac`, `kite`) to prepare files for.
#' @param overwrite Overwrite existing kallisto index.
#' @param f1 Path to the cDNA FASTA (for `lamanno`, `nac` workflows) or mismatch FASTA (for `kite`) to be generated.
#' @param f2 [Required for `lamanno` and `nac` workflows] Path to the intron FASTA to be generated.
#' @param c1 [Required for `lamanno` and `nac` workflows] Path to generate cDNA transcripts-to-capture.
#' @param c2 [Required for `lamanno` and `nac` workflows] Path to generate intron transcripts-to-capture.
#'
#' @return Invisible list containing the stdout and stderr from the `kb ref` command execution.
#' @examples
#' \dontrun{
#' kb_ref(index = "path/to/index", t2g = "path/to/t2g.txt", fasta = "path/to/genome.fa",
#'        gtf = "path/to/genes.gtf", workflow = "standard", overwrite = TRUE)
#' }
#' @export
#'
#' @importFrom utils system2
kb_ref <- function(index, t2g, fasta = NULL, gtf = NULL, feature = NULL, tmp = NULL, keep_tmp = FALSE, 
                   verbose = FALSE, n = NULL, species = NULL, k = NULL, 
                   workflow = NULL, overwrite = FALSE, f1 = NULL, f2 = NULL, c1 = NULL, c2 = NULL) {

  #make output folders if they don't exist
    if (!dir.exists(dirname(index))) dir.create(dirname(index), recursive = TRUE)
    if (!dir.exists(dirname(t2g))) dir.create(dirname(t2g), recursive = TRUE)
  
  # Validate `species`
  valid_species <- c("human", "mouse", "linnarsson")
  if (!is.null(species) && !species %in% valid_species) {
    stop("Invalid species value. Must be one of: ", paste(valid_species, collapse = ", "), ".")
  }
  
  # Validate `workflow`
  valid_workflows <- c("standard", "lamanno", "nac", "kite")
  if (!is.null(workflow) && !workflow %in% valid_workflows) {
    stop("Invalid workflow value. Must be one of: ", paste(valid_workflows, collapse = ", "), ".")
  }
  
  # Adjust initial command construction for optional fasta and gtf
  args <- c("ref", "-i", shQuote(index), "-g", shQuote(t2g))
  if (!is.null(fasta)) args <- c(args, "-f1", shQuote(fasta))
  if (!is.null(gtf)) args <- c(args, "-gtf", shQuote(gtf))
  
  # Adjustments for other optional arguments
  if (!is.null(feature)) args <- c(args, "--feature", shQuote(feature))
  if (!is.null(tmp)) args <- c(args, "--tmp", shQuote(tmp))
  if (keep_tmp) args <- c(args, "--keep-tmp")
  if (verbose) args <- c(args, "--verbose")
  if (!is.null(n)) args <- c(args, "-n", as.character(n))
  if (!is.null(species)) args <- c(args, "-d", shQuote(species))
  if (!is.null(k)) args <- c(args, "-k", as.character(k))
  if (!is.null(workflow)) args <- c(args, "--workflow", shQuote(workflow))
  if (overwrite) args <- c(args, "--overwrite")
  if (!is.null(f1)) args <- c(args, "-f1", shQuote(f1))
  if (!is.null(f2)) args <- c(args, "-f2", shQuote(f2))
  if (!is.null(c1)) args <- c(args, "-c1", shQuote(c1))
  if (!is.null(c2)) args <- c(args, "-c2", shQuote(c2))
  
  # Execute the command
  result <- system2("kb", args, stdout = TRUE, stderr = TRUE)
  
  # update mtime so that files are not purged from /scratch
  system(sprintf("find '%s' -type f -exec touch {} +", dirname(here::here(index))))
  
  here::here(index)
}


#--------------------------------------------------------------------------------

#' Generate count matrices from single-cell FASTQ files
#'
#' This function wraps the `kb count` command to generate count matrices from a set of single-cell FASTQ files.
#' It supports various workflows including standard, lamanno, nac, kite, and kite:10xFB.
#'
#' @param index Path to kallisto index/indices, comma-delimited.
#' @param t2g Path to transcript-to-gene mapping.
#' @param technology Single-cell technology used (run `kb --list` to view).
#' @param fastqs Vector of paths to FASTQ files.
#' @param tmp Override default temporary directory.
#' @param keep_tmp Do not delete the tmp directory.
#' @param verbose Print debugging information.
#' @param out Path to output directory (default: current directory).
#' @param whitelist Path to file of whitelisted barcodes to correct to.
#' @param threads Number of threads to use (default: 8).
#' @param memory Maximum memory used (default: 4G).
#' @param workflow Type of workflow (standard, lamanno, nac, kite, kite:10xFB).
#' @param mm Include reads that pseudoalign to multiple genes.
#' @param tcc Generate a TCC matrix instead of a gene count matrix.
#' @param filter Produce a filtered gene count matrix (default: bustools).
#' @param c1 Path to cDNA transcripts-to-capture (required for lamanno and nac workflows).
#' @param c2 Path to intron transcripts-to-captured (required for lamanno and nac workflows).
#' @param overwrite Overwrite existing output.bus file.
#' @param dry_run Perform a dry run without executing.
#' @param loom Generate loom file from count matrix.
#' @param h5ad Generate h5ad file from count matrix.
#' @return Invisible NULL. The function is called for its side effects.
#' @examples
#' \dontrun{
#' kb_count(index = "path/to/index", t2g = "path/to/t2g", technology = "10xv3", 
#'          fastqs = c("fastq1.fastq.gz", "fastq2.fastq.gz"))
#' }
#' @export
kb_count <- function(index, t2g, technology, fastqs, tmp = NULL, keep_tmp = FALSE, 
                     verbose = FALSE, out = NULL, whitelist = NULL, threads = NULL, 
                     memory = NULL, workflow = NULL, mm = FALSE, tcc = FALSE, 
                     filter = "bustools", c1 = NULL, c2 = NULL, overwrite = FALSE, 
                     dry_run = FALSE, loom = FALSE, h5ad = FALSE) {

  # if kb_info.json exits, return its path
  # else zap the output directory if overwrite = TRUE
  if (file.exists(glue::glue("{out}/kb_info.json"))) {
      return(glue::glue("{out}/kb_info.json"))
  } else {
      if (dir.exists(out)) {
          if (overwrite) unlink(out, recursive = TRUE)
          if (!overwrite) stop("Output directory already exists. Set overwrite = TRUE to overwrite.")
      }
  }
  dir.create(out, recursive = TRUE)

  # Resources
  if (nzchar(Sys.getenv("SLURM_JOB_ID"))) {
    alloc <- hprcc::slurm_allocation()
    memory = glue::glue("{alloc$Memory_GB - 10}G")
    threads = alloc$CPUs}
  
  args <- c("count", "-i", shQuote(index), "-g", shQuote(t2g), "-x", shQuote(technology), fastqs)
  
  if (!is.null(tmp)) args <- c(args, "--tmp", shQuote(tmp))
  if (keep_tmp) args <- c(args, "--keep-tmp")
  if (verbose) args <- c(args, "--verbose")
  if (verbose) log <- glue::glue(" | tee {out}/kb_count.log")
  if (!is.null(out)) args <- c(args, "-o", shQuote(out))
  if (!is.null(whitelist)) args <- c(args, "-w", shQuote(whitelist))
  args <- c(args, "-t", as.character(threads), "-m", shQuote(memory))
  if (!is.null(workflow)) args <- c(args, "--workflow", shQuote(workflow))
  if (mm) args <- c(args, "--mm")
  if (tcc) args <- c(args, "--tcc")
  if (!is.null(filter)) args <- c(args, "--filter", shQuote(filter))
  if (!is.null(c1)) args <- c(args, "-c1", shQuote(c1))
  if (!is.null(c2)) args <- c(args, "-c2", shQuote(c2))
  if (overwrite) args <- c(args, "--overwrite")
  if (dry_run) args <- c(args, "--dry-run")
  if (loom) args <- c(args, "--loom")
  if (h5ad) args <- c(args, "--h5ad")
 
  # which kb
  #kb <- Sys.which("kb")
  kb<-"/home/domeally/miniconda3/bin/kb"
  # Constructing the command string for verbose output
  cmd_str <- paste(kb, paste(c(args, log), collapse = " "))

  # If verbose is TRUE, echo the command
  if (verbose) {
    message("Executing system command: ", cmd_str)
  }
  
  # Execute the command
  result <- system2(kb, args, stdout = TRUE, stderr = TRUE)
  Log_message(result)
  # Check for a result and return
  if (file.exists(glue::glue("{out}/kb_info.json"))) {
        return(glue::glue("{out}/kb_info.json"))
  } else {
        return(NULL) 
  }
  

#--------------------------------------------------------------------------------
#' Display package and citation information for kb
#'
#' Executes the `kb info` command to display version and citation information for kb, kallisto, and bustools.
#'
#' @return Prints version information and citations for kb and its dependencies.
#' @examples
#' kb_info()
#' @export
kb_info <- function() {
  result <- system2("kb", args = "info", stdout = TRUE, stderr = TRUE)
  return(result)
}

Log_message <- function(message_str, log = TRUE, log_dir = "logs") {

    # Get the SLURM job ID, or set to 'NA' if not under SLURM control
    job_id <- Sys.getenv("SLURM_JOB_ID", unset = "NA")

    log_dir_path <- here::here(log_dir)
    if (!dir.exists(log_dir_path)) {
        dir.create(log_dir_path, recursive = TRUE, showWarnings = FALSE)
    }

    log_file <- glue::glue("{log_dir_path}/job_{job_id}.log")

    message(message_str)

    if (isTRUE(log)) {
        writeLines(message_str, con = file(log_file, open = "a"))
    }
}
