# scRNA-seq Analysis of Insulin Gene Variants in Type 1 Diabetes

## Research Overview

This repository contains code and analysis pipelines for the single-cell RNA sequencing (scRNA-seq) component of our research on genetic protection from Type 1 diabetes (T1D). Our analysis reveals how the protective insulin gene (*INS*) variant (rs3842752) affects beta-cell stress and gene expression at the single-cell level.

### Key Analysis Findings

- Beta cells carrying the protective *INS* variant show an inverse correlation between insulin mRNA levels and cellular stress markers
- This correlation is not present in beta cells carrying only the susceptible *INS* variant
- Beta cells with the protective variant show significantly reduced expression of TGM2 (tissue transglutaminase-2), which is associated with T1D immunopathogenesis
- Differential expression analysis identified additional genes associated with beta-cell stress and function

### Citation

If you use these analysis methods, please cite our paper:

> van Tienhoven, R., O'Meally, D., Scott, T.A., Morris, K.V., Williams, J.C., Kaddis, J.S., Zaldumbide, A., Roep, B.O. (2025). Genetic protection from type 1 diabetes resulting from accelerated insulin mRNA decay. Cell, 188, 1-10. https://doi.org/10.1016/j.cell.2025.02.018

## Repository Contents

[Targets](https://docs.ropensci.org/targets/) pipeline for processing scRNA-seq data from the Human Pancreas Analysis Program ([HPAP](https://hpap.pmacs.upenn.edu/)). The workflow includes targets for

- Cellranger counts and library quality control ([nf-core/scrnaseq](https://nf-co.re/scrnaseq/2.7.1))
- SNP genotyping using CellSNP-lite
- Ambient RNA removal with CellBender
- Per cell quality control & doublet detection with scDblFinder
- Cell type annotation using SingleR with reference atlas from [Elgamal et al. 2023](https://doi.org/10.2337/db23-0130)
- XBP1 splicing analysis and stress scoring (not used, insufficient read depth)
- INS expression normalization using housekeeping genes
- Integration of donor libraries using Harmony
- Publication figures 3 and S2

### Resource Files

- [MultiQC Report v2](assets/multiqc_report_v2.html) - nf-core/scrnaseq MultiQC report for 10X V2 chemistry libraries
- [MultiQC Report v3](assets/multiqc_report_v3.html) - nf-core/scrnaseq MultiQC report for 10X V3 chemistry libraries
- [1000 Genomes SNP data](assets/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) - Reference variants for genotyping, from [CellSNP](http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz)

The processed Seurat object (~4.5Gb) is available via Globus. Make your request by raising an [issue](https://github.com/KaddisLab/2025_van_tienhoven_et_al_cell/issues).

## Data Sources

The scRNA-seq data analyzed in this repository were obtained from the Human Pancreas Analysis Program (HPAP) Database (https://hpap.pmacs.upenn.edu/). Computational infrastructure was provided by the City of Hope High Performance Research Computing Center.

## Technical Requirements

### Computing Environment
- SLURM-based HPC system
- R version 4.3 or higher
- Bioconductor 3.19+

### R Packages
- Seurat v5.1.0+
- Harmony v0.2.0+
- tidyverse & tidyomics framework
- other packages as documented in the code

### Installation

```bash
git clone https://github.com/KaddisLab/2025_van_tienhoven_et_al_cell.git
```

## Running the Pipeline

```R
targets::tar_make()
```

Note: Set `analysis_cache` in `_targets.R` to a folder on your local system.

## Project Structure

- `R/`: Custom R functions and analysis utilities
- `assets/`: Resource files including reference data and QC reports
- `figures/`: Generated figures and visualizations
- `_targets/`: Target pipeline outputs

## Contact

For inquiries regarding the single-cell analysis, please contact:
- Denis O'Meally (data analysis): [@drejom](https://github.com/drejom)
- Bart O. Roep (principal investigator): boroep@lumc.nl

## Collaborating Institutions

- Leiden University Medical Center, Netherlands
- City of Hope, Duarte, CA, USA
- Queensland University of Technology, Australia
