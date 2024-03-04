# Pipeline for Processing HPAP scRNAseq data

This repository contains the R script `_targets.R` for processing 10X samples from [PancDB](https://hpap.pmacs.upenn.edu/about-pancdb). 

## Requirements

Set `analysis_cache` in `_targets.R` to a folder on Apollo or Gemeni.

## Description

The script sets up a pipeline for processing HPAP data. It defines the paths for the `analysis_cache` and .... The `analysis_cache` path is determined by the cluster, and the transcriptome path is always located in the `analysis_cache`.

The script also sets some options for the {targets}, such as the packages to load and the error handling behaviour.

## Usage

To run the script, clone this repository to your $HOME folder on _Apollo_ or _Gemini_: 

```{sh}
git clone https://github.com/cohmathonc/mmu_10X_preleukemic.git
cd mmu_10X_preleukemic
```

To run the pipeline in a SLURM session, submit the jobscript:

```{sh}
sbatch --chdir=$(pwd) _targets.job
```

You can submit the script from within a SSH shell on Apollo, or from the terminal in RStudio running on Apollo.

To run interactively, open an R session in the folder and launch the targets pipeline:

```{r}
targets::tar_make()
```

[This issue](https://github.com/cohmathonc/hprcc/issues/9) is currently preventing the pipeline from running interactively in RStudio on Apollo.

### Contributing

Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests.

### License

This project is licensed under the MIT License - see the LICENSE.md file for details
