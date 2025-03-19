# scRNA-seq Analysis of Insulin Gene Variants in Type 1 Diabetes

## Project Status

This repository is currently under active development. Code and documentation for single-cell RNA sequencing analysis related to our recent Cell publication will be updated shortly. **'Watch' this repo to get notified when updates are posted or raise an [issue](https://github.com/KaddisLab/2025_van_tienhoven_et_al_cell/issues)**

## Research Overview

This repository contains code and analysis pipelines for the single-cell RNA sequencing (scRNA-seq) component of our research on genetic protection from Type 1 diabetes (T1D). Our analysis reveals how the protective insulin gene (INS) variant (rs3842752) affects beta-cell stress and gene expression at the single-cell level.

### Key Analysis Findings

- Beta cells carrying the protective INS variant show an inverse correlation between insulin mRNA levels and cellular stress markers
- This correlation is not present in beta cells carrying only the susceptible INS variant
- Beta cells with the protective variant show significantly reduced expression of TGM2 (tissue transglutaminase-2), which is associated with T1D immunopathogenesis
- Differential expression analysis identified additional genes associated with beta-cell stress and function

### Citation

If you use these analysis methods, please cite our paper:

> van Tienhoven, R., O'Meally, D., Scott, T.A., Morris, K.V., Williams, J.C., Kaddis, J.S., Zaldumbide, A., Roep, B.O. (2025). Genetic protection from type 1 diabetes resulting from accelerated insulin mRNA decay. Cell, 188, 1-10. https://doi.org/10.1016/j.cell.2025.02.018

## Repository Contents

This repository will include:

- Preprocessing workflow for scRNA-seq data from human pancreatic islets
- Code for INS genotype analysis from scRNA-seq data
- Scripts for identifying nascent vs mature INS transcript counts
- Differential gene expression analysis between INS variants
- Beta-cell stress score calculation methodology
- Visualization code for UMAP, feature plots, and correlation analyses

## Data Sources

The scRNA-seq data analyzed in this repository were obtained from the Human Pancreas Analysis Program (HPAP) Database (https://hpap.pmacs.upenn.edu/). Data processing utilized the Human Islet Research Network (HIRN) consortium resources.

## Requirements

The analysis pipelines require:
- R v4.3+
- Bioconductor packages
- Seurat v5.1.0+
- Harmony v0.2.0+
- tidyomics framework

Detailed environment specifications and package versions will be provided.

## Contact

For inquiries regarding the single-cell analysis, please contact:
- Denis O'Meally (data analysis): [@drejom](https://github.com/drejom)
- Bart O. Roep (principal investigator): boroep@lumc.nl

## Collaborating Institutions

- Leiden University Medical Center, Netherlands
- City of Hope, Duarte, CA, USA
- Queensland University of Technology, Australia
