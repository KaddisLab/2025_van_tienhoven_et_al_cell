library(targets)
library(tarchetypes)

# HPC Support - optional
# https://github.com/cohmathonc/hprcc
options(hprcc.slurm_logs = TRUE, hprcc.slurm_jobs = TRUE)
library(hprcc)

tar_source()

tar_option_set(
    packages = c("tidyverse", "Seurat"),
    error = "abridge",
    workspace_on_error = TRUE,
)

list(
    # check fastq md5sums
    tar_target(fastq_md5sums, read_md5sums(here::here(glue::glue("{analysis_cache}/data/hpapdata/fastq"))), resources = small),
    tar_target(md5sums_check, check_md5sum(fastq_md5sums), map(fastq_md5sums), resources = tiny),
    # import metadata --------------------------------------------------------
    tar_target(pancdb_metadata, get_pancdb_metadata(), deployment = "main"),
    ## get fastq files --------------------------------------------------------
    tar_target(
        # Exclude Fluidigm fastq files
        fastq_10x,
        grep(
            "SSq2",
            list.files(
                path = glue::glue("{analysis_cache}/data/hpapdata/fastq"),
                pattern = "R[12](_001)?(_fastq-data)?\\.fastq\\.gz$",
                full.names = TRUE,
                recursive = TRUE
            ),
            value = TRUE,
            invert = TRUE,
            perl = TRUE
        ),
        deployment = "main"
    ),
    tar_target(
        fastq_10xv2,
        fastq_10x[grep(
            pancdb_metadata$donor_id[
                pancdb_metadata$reagent_kit == "10X-Chromium-GEX-3p-v2"
            ] |>
                na.omit() |>
                paste(collapse = "|"),
            fastq_10x
        )],
        deployment = "main"
    ),
    tar_target(
        fastq_10xv3,
        fastq_10x[grep(
            pancdb_metadata$donor_id[
                pancdb_metadata$reagent_kit %in%
                    c("10X-Chromium-GEX-3p-v3", "10X-Chromium-GEX-3p-v3.1")
            ] |>
                na.omit() |>
                paste(collapse = "|"),
            fastq_10x
        )],
        deployment = "main"
    ),
    tar_target(
        fastq_ss2,
        list.files(
            path = glue::glue("{analysis_cache}/data/hpapdata/fastq"),
            pattern = "_scRNA_SSq2.*\\.fastq\\.gz$",
            full.names = TRUE,
            recursive = TRUE
        ),
        deployment = "main"
    ),
    tar_target(
        fastq_ss3,
        list.files(
            path = glue::glue("{analysis_cache}/data/hpapdata/fastq"),
            pattern = "_scRNA_\\d{5}_.*\\.fastq\\.gz$",
            full.names = TRUE,
            recursive = TRUE
        ),
        deployment = "main"
    ),
    ## Cellranger -----------------------------------------------------------------------
    # Make nf-core sample sheet
    tar_target(
        sample_sheet_10xv2,
        hpap_fastq_to_sample_sheet(fastq_10xv2),
        deployment = "main"
    ),
    tar_target(
        sample_sheet_10xv3,
        hpap_fastq_to_sample_sheet(fastq_10xv3),
        deployment = "main"
    ),
    # rename fastq files as per cellranger requirements
    tar_target(
        renamed_sample_sheet_10xv2,
        rename_sample_sheet(sample_sheet_10xv2),
        deployment = "main"
    ),
    tar_target(
        renamed_sample_sheet_10xv3,
        rename_sample_sheet(sample_sheet_10xv3),
        deployment = "main"
    ),
    # run nfcore-scrnaseq-cellranger
    tar_target(
        nfcore_scrnaseq_multiqc_10xv2,
        run_nfcore_scrnaseq(
            run_folder = "cellranger_out/10xv2",
            sample_sheet = renamed_sample_sheet_10xv2,
            protocol = "10XV2"
        ),
        deployment = "main",
    ),
    tar_target(
        nfcore_scrnaseq_multiqc_10xv3,
        run_nfcore_scrnaseq(
            run_folder = "cellranger_out/10xv3",
            sample_sheet = renamed_sample_sheet_10xv3,
            protocol = "10XV3"
        ),
        deployment = "main"
    ),
    tar_target(
        cellranger_run_folders_10xv3,
        list.files(
            gsub(
                "multiqc/multiqc_report.html",
                "cellranger/count",
                nfcore_scrnaseq_multiqc_10xv3
            ),
            pattern = "HPAP",
            full.names = TRUE
        ),
        deployment = "main"
    ),
    tar_target(
        cellranger_run_folders_10xv2,
        list.files(
            gsub(
                "multiqc/multiqc_report.html",
                "cellranger/count",
                nfcore_scrnaseq_multiqc_10xv2
            ),
            pattern = "HPAP",
            full.names = TRUE
        ),
        deployment = "main"
    ),
    tar_target(
        cellranger_run_folders,
        c(cellranger_run_folders_10xv3, cellranger_run_folders_10xv2),
        deployment = "main"
    ),
    ## STARSolo -------------------------------------------------------------------
    # Run STARsolo on cellranger BAMs 
    tar_target(
        starsolo_10xv2,
        run_starsolo(
            technology = "10XV2",
            ref_dir = "/ref_genomes/cellranger/human/star_cr_ref",
            cellranger_folder = cellranger_run_folders_10xv2
        ),
        resources = large,
        pattern = map(cellranger_run_folders_10xv2),
        format = "file_fast"
    ),
    tar_target(
        starsolo_10xv3,
        run_starsolo(
            technology = "10XV3",
            ref_dir = "/ref_genomes/cellranger/human/star_cr_ref",
            cellranger_folder = cellranger_run_folders_10xv3
        ),
        resources = large,
        pattern = map(cellranger_run_folders_10xv3),
        format = "file_fast"
    ),
    tar_target(
        starsolo_alltech,
        c(starsolo_10xv2, starsolo_10xv3),
        deployment = "main"
    ),
    ## XBP1u score -------------------------------------------------------------------
    # per sample XBP1u score
    tar_target(
        xbp1_psi_per_sample,
        calculate_xbp1_psi_per_sample(starsolo_alltech),
        resources = small,
    ),
    # per cell XBP1u score
    tar_target(
        xbp1_psi_per_cell,
        calculate_xbp1_psi_per_cell(starsolo_alltech),
        resources = tiny,
        pattern = map(starsolo_alltech)
    ),
    # percent spliced per cell
    tar_target(
        percent_spliced_per_cell_INS,
        percent_spliced_per_cell(starsolo_alltech, "INS", housekeeping_genes),
        resources = tiny,
        pattern = map(starsolo_alltech)
    ),
    tar_target(
        percent_spliced_per_cell_XBP1,
        percent_spliced_per_cell(starsolo_alltech, "XBP1"),
        resources = tiny,
        pattern = map(starsolo_alltech)
    ),
    tar_target(
        percent_spliced_per_cell_GAPDH,
        percent_spliced_per_cell(starsolo_alltech, "GAPDH"),
        resources = tiny,
        pattern = map(starsolo_alltech)
    ),
    ## Download SNPs ---------------------------------------------------------------------
    tar_target(
        snp_vcf,
        {
            destfile <- glue::glue(
                "{analysis_cache}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
            )
            download.file(
                url = "http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz",
                destfile = destfile,
                mode = "wb"
            )
            return(destfile)
        },
        deployment = "main",
        format = "file_fast"
    ),
    ## CellSNP-lite ---------------------------------------------------------------------
    tar_target(
        cellsnp_lite,
        run_cellsnp_lite(
            cellranger_run_folder = cellranger_run_folders,
            region_vcf = glue::glue(
                "{analysis_cache}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
            )
        ),
        deployment = "main",
        pattern = map(cellranger_run_folders)
    ),
    # parse sample genotypes
    tar_target(
        protected_cohort,
        get_cellsnp_lite_genotypes("11\t2159843", "rs3842752"),
        deployment = "main"
    ),
    tar_target(
        rs3842753_cohort,
        get_cellsnp_lite_genotypes("11\t2159830", "rs3842753"),
        deployment = "main"
    ),
    tar_target(
        rs689_cohort,
        get_cellsnp_lite_genotypes("11\t2160994", "rs689"),
        deployment = "main"
    ),
    tar_target(
        rs13266634_cohort,
        get_cellsnp_lite_genotypes("8\t117172544", "rs13266634"),
        deployment = "main"
    ),
    # parse cell genotypes
    tar_target(cell_genotypes, get_cellsnp_lite_cell_genotypes(cellsnp_lite)),
    ## Collate sample metadata -------------------------------------------------
    tar_target(
        pancdb_metadata_agg,
        collate_sample_metadata(
            pancdb_metadata,
            protected_cohort,
            rs3842753_cohort,
            rs689_cohort,
            xbp1_psi_per_sample
        ),
        deployment = "main"
    ),
    ## Cellbender using CellRanger counts ------------------------------------------------------------
    tar_target(
        cellbender_h5,
        run_cellbender(cellranger_run_folders),
        pattern = map(cellranger_run_folders),
        format = "file_fast",
    ),
    tar_target(
        cellbender_seurat_objects,
        make_seurat_cellbender(
            cellbender_h5,
            cellranger_run_folders,
            pancdb_metadata_agg
        ),
        pattern = map(cellbender_h5, cellranger_run_folders),
        format = "file_fast"
    ),
    tar_target(
        cellbender_qc_plots,
        seurat_plot_cellbender(cellbender_seurat_objects, "INS"),
        pattern = slice(map(cellbender_seurat_objects), 1L),
        format = "file_fast",
        resources = small
    ),
    ## INS housekeeping normalisation -----------------------------------------------------
    tar_target(
        INS_hknorm_per_cell,
        seurat_INS_hknorm(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        resources = tiny
    ),
    ## UCell Signature scores ---------------------------------------------------------------------
    tar_target(
        sig_scores_per_cell,
        seurat_sig_scores(cellbender_seurat_objects, features = signatures),
        pattern = map(cellbender_seurat_objects),
        resources = medium
    ),
    ## SingleR cell type annotation --------------------------------------------------------------
    # Cell atlas from Elgamal et al. 2023
    # https://doi.org/10.2337/db23-0130
    # https://islet-hpap.s3.us-west-2.amazonaws.com/hpap_islet_scRNAseq.rds
    tar_target(
        elgamal_etal_seurat_object,
        glue::glue(
            "{analysis_cache}/data/hpapdata/gaulton_lab/hpap_islet_scRNAseq.rds"
            # https://islet-hpap.s3.us-west-2.amazonaws.com/hpap_islet_scRNAseq.rds
        ),
        format = "file_fast",
        deployment = "main"
    ),
    tar_target(
        elgamal_cell_type_csv,
        seurat_singleR_transfer_label(
            cellbender_seurat_objects,
            elgamal_etal_seurat_object,
            cell_type_col = "Cell Type"
        ),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = medium
    ),
    ## Cell cycle annotation --------------------------------------------------------------
    tar_target(
        cell_cycle_csv,
        seurat_cell_cycle(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = small
    ),
    ## Doublet annotation --------------------------------------------------------------
    tar_target(
        scDblFinder_csv,
        seurat_scDblFinder(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = medium
    ),
    ## ddqc -------------------------------------------------------------------------------
    tar_target(
        ddqc_seurat_objects,
        seurat_ddqc(
            cellbender_seurat_objects,
            scDblFinder_csv,
            use_mad_mito = FALSE,
            percent.mito.upper.bound = 30L
        ),
        pattern = map(cellbender_seurat_objects, scDblFinder_csv),
        format = "file_fast",
        resources = small
    ),
    ## Convert seurat_objects to BPcells matrices -------------------------------------------
    tar_target(
        ddqc_bpcells_all,
        seurat_to_bpcells(ddqc_seurat_objects),
        format = "file_fast",
        pattern = map(ddqc_seurat_objects)
    ),

    ## Drop failed samples before merge/clustering/DEG ------------------------------------
    # First, create a target with a vector of valid samples
    targets::tar_target(
        ddqc_bpcells,
        {
            all_paths <- ddqc_bpcells_all
            stringr::str_subset(
                all_paths,
                stringr::str_c(failed_qc_donor_ids, collapse = "|"),
                negate = TRUE
            )
        },
        iteration = "vector"
    ),

    ## Aggregate cell annotation -------------------------------------------------------
    tar_target(
        aggregated_cell_annot_csv,
        seurat_collate_cell_annot(
            ddqc_seurat_objects,
            elgamal_cell_type_csv,
            scDblFinder_csv,
            cell_cycle_csv,
            INS_hknorm_per_cell,
            xbp1_psi_per_cell,
            percent_spliced_per_cell_INS,
            percent_spliced_per_cell_XBP1,
            percent_spliced_per_cell_GAPDH,
            sig_scores_per_cell
        ),
        resources = small
    ),
    ## Log-normalised counts --------------------------------------------------------------
    tar_target(
        seurat_object_merged,
        seurat_merge(ddqc_bpcells, "all_merged"),
        resources = large
    ),
    tar_target(
        seurat_object_lognorm,
        seurat_lognorm(
            seurat_object_merged,
            assay = "RNA",
            vars_to_regress = NULL
        ),
        resources = large
    ),
    tar_target(
        seurat_object_lognorm_ann_samples,
        seurat_annotate_samples(seurat_object_lognorm, pancdb_metadata_agg),
        resources = large
    ),
    tar_target(
        seurat_object_lognorm_harmony,
        seurat_harmony(
            seurat_object_lognorm_ann_samples,
            "RNA",
            c("tissue_source", "reagent_kit", "orig.ident")
        ),
        resources = large_mem
    ),
    tar_target(
        seurat_object_lognorm_annotated,
        seurat_annotate_cells(
            seurat_object_lognorm_harmony,
            aggregated_cell_annot_csv
        ),
        resources = large
    ),
    ## Sparse matrix (portable, no BPCells links) --------------------------------
    tar_target(
        seurat_object_lognorm_annotated_sparse,
        seurat_bpcells2sparse(seurat_object_lognorm_annotated),
        resources = large
    ),
    ## Quarto reports for figures ------------------------------------------------
    tar_quarto(figure3, "figure3.qmd"),
    tar_quarto(figureS2, "figureS2.qmd"),
)
