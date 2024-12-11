# Pipeline for Processing HPAP scRNAseq Data

This repository contains R scripts and analysis pipelines for processing single-cell RNA sequencing data from [PancDB](https://hpap.pmacs.upenn.edu/about-pancdb) pancreatic tissue samples. The analysis is performed using R and various Bioconductor packages.

## Project Overview

The project focuses on analyzing single-cell RNA sequencing data with particular attention to:
- INS gene expression analysis
- Differential expression analysis
- Cell type-specific analyses
- Housekeeping gene normalization

## Requirements

### Computing Environment
The analysis is configured to run on SLURM-based HPC systems with appropriate resource allocation.

### Prerequisites
- R version 4.3 or higher
- Bioconductor 3.19
- Key packages: Seurat, muscat, tidyverse, ComplexHeatmap

### Installation

The project uses `renv` for package management with a shared cache configuration:

```
git clone https://github.com/cohmathonc/mmu_10X_preleukemic.git
cd mmu_10X_preleukemic
Rscript -e 'renv::restore()'
```

## Project Structure

- `R/`: Custom R functions and utilities
- `test_*.r`: Analysis scripts for different aspects of the project
- `.Rprofile`: Project-specific R configuration
- `_targets/`: Target pipeline outputs
- Set `analysis_cache` in `_targets.R` to a folder on you local system

## Usage

### Running with SLURM
To run the pipeline in a SLURM session, submit the jobscript:

```
sbatch _targets.job
```

### Running Interactively
To run interactively, open an R session in the folder and launch the targets pipeline:

```
targets::tar_make()
```
## Data
Data files should be placed in the appropriate directories as specified in the targets pipeline.

## Contributing
Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests.

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.

## Contact
Denis O'Meally (@drejom)
