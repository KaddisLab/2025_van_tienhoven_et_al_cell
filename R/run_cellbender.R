#' @title Run CellBender
#'
#' @description Runs CellBender using the 10x h5 file as input
#' @param rawFeatureMatrix The path to 10x's raw_feature_bc_matrix.h5
#' @param expectedCells Passed to CellBender --expected-cells
#' @param totalDropletsIncluded Passed to CellBender --total-droplets-included
#' @param fpr Passed to CellBender --fpr
#' @param epochs Passed to CellBender --epochs
#'
#' @export
run_cellbender <- function(run_folder, expectedCells = 5000, totalDropletsIncluded = 20000, fpr = 0.01, epochs = 150) {
  # Check if CellBender is in the PATH
  cellBenderPath <- Sys.which("cellbender")
  if (nzchar(cellBenderPath)) {
    print(paste("Running CellBender from:", cellBenderPath))
  } else {
    stop("CellBender is not found in the PATH. Please ensure CellBender is installed and accessible.")
  }
run_folder<-"/scratch/domeally/DCD.tienhoven_scRNAseq.2024/cellranger_out/10xv2/nfcore-scrnaseq-v2.5.1-refdata-gex-GRCh38-2020-A/cellranger/count/HPAP-019"
  # Define the output file path
  outH5File <- tempfile(fileext = '.h5')
  rawFeatureMatrix <- paste0(run_folder, "outs/raw_feature_bc_matrix.h5")


  # Construct CellBender command
  cmd <- sprintf("cellbender remove-background --input %s --output %s --expected-cells %d --total-droplets-included %d --fpr %f --epochs %d",
                 shQuote(rawFeatureMatrix),
                 shQuote(outH5File),
                 expectedCells,
                 totalDropletsIncluded,
                 fpr,
                 epochs)

  # Check if CUDA is available and append --cuda if truesbatch --chdir=$(pwd)
  cudaCheckCmd <- "/usr/bin/python3 -c 'import torch; print(\"--cuda\" if torch.cuda.is_available() else \"\")'"
  cudaAvailable <- system(cudaCheckCmd, intern = TRUE)
  if (cudaAvailable == "--cuda") {
    cmd <- paste(cmd, "--cuda")
  }

  # Execute the CellBender command
  system(cmd)

  # Check for the existence of the output file
  outputFiltered <- sub(pattern = '.h5$', replacement = '_filtered.h5', outH5File)
  if (!file.exists(outputFiltered)) {
    stop(paste0('Missing file: ', outputFiltered))
  }

  # Load the filtered data into a Seurat object
  seuatRawData <- Seurat::Read10X_h5(filename = outputFiltered, use.names = TRUE)
  print(paste0('Cells in CellBender filtered matrix: ', ncol(seuatRawData)))

  return(seuatRawData)
}


srun -n 6 -N 1-1 --mem=60G --gres=gpu:1 -p gpu --time=4:00:00 --pty bash

library(crew)
library(hprcc)
controller <- 

create_controller(
  "gpu",
  slurm_cpus = 6,
  slurm_mem_gigabytes = 60,
  slurm_walltime_minutes = 360,
  slurm_workers = 35,
  slurm_partition = "gpu",
  slurm_log_dir = "logs"
)

crew_controller_local(
  name = "example",
  workers = 2,
  seconds_idle = 10
)