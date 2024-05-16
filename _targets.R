library(targets)
library(tarchetypes)
library(hprcc)

tar_source()

tar_option_set(
    packages = c("tidyverse"),
    error = "abridge",
    workspace_on_error = TRUE
)

list(
    # check fastq md5sums
    # tar_target(fastq_md5sums, read_md5sums(here::here(glue::glue("{analysis_cache}/data/hpapdata/fastq"))), resources = small),
    # TODO vectorise this function and make groups of 20 or more files to avoid overhead of calling md5sum for each file
    # TODO look at tar_rep() for this...
    # tar_target(md5sums_check, check_md5sum(fastq_md5sums), map(fastq_md5sums), resources = tiny),
    # import metadata --------------------------------------------------------
    tar_target(pancdb_metadata, get_pancdb_metadata(), deployment = "main"),
    # get fastq files --------------------------------------------------------
    tar_target(
        fastq_10x,
        grep("SSq2", list.files(path = glue::glue("{analysis_cache}/data/hpapdata/fastq"), pattern = "R[12](_001)?(_fastq-data)?\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
            value = TRUE, invert = TRUE, perl = TRUE
        ),
        deployment = "main"
    ),
    tar_target(
        fastq_10xv2,
        fastq_10x[grep(pancdb_metadata$donor_id[pancdb_metadata$reagent_kit == "10X-Chromium-GEX-3p-v2"] |> na.omit() |> paste(collapse = "|"), fastq_10x)],
        deployment = "main"
    ),
    tar_target(
        fastq_10xv3,
        fastq_10x[grep(pancdb_metadata$donor_id[pancdb_metadata$reagent_kit %in% c("10X-Chromium-GEX-3p-v3", "10X-Chromium-GEX-3p-v3.1")] |> na.omit() |> paste(collapse = "|"), fastq_10x)],
        deployment = "main"
    ),
    tar_target(
        fastq_ss2,
        list.files(path = glue::glue("{analysis_cache}/data/hpapdata/fastq"), pattern = "_scRNA_SSq2.*\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        deployment = "main"
    ),
    tar_target(
        fastq_ss3,
        list.files(path = glue::glue("{analysis_cache}/data/hpapdata/fastq"), pattern = "_scRNA_\\d{5}_.*\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        deployment = "main"
    ),

    # kallisto-bustools ---------------------------------------------------------------------
    # Using cellranger reference genome as index
    # Group fastqs by sample
    tar_target(fastq_10xv3_by_sample, arrange_by_sample_10x(fastq_10xv3), deployment = "main"),
    tar_target(fastq_10xv2_by_sample, arrange_by_sample_10x(fastq_10xv2), deployment = "main"),
    tar_target(fastq_ss2_by_sample, arrange_by_sample_ssq(fastq_ss2), deployment = "main"),
    tar_target(fastq_ss3_by_sample, arrange_by_sample_ssq(fastq_ss3), deployment = "main"),
    # run kb count - tcc mode
    tar_target(
        kb_count_10xv3_tcc,
        kb_count(
            technology = "10XV3", tcc = TRUE,
            index = "/ref_genomes/cellranger/human/kb_cr_ref/nac_transcriptome.idx",
            t2g = "/ref_genomes/cellranger/human/kb_cr_ref/nac_t2g.txt",
            c1 = "/ref_genomes/cellranger/human/kb_cr_ref/nac_cdna_t2c.txt",
            c2 = "/ref_genomes/cellranger/human/kb_cr_ref/nac_nascent_t2c.txt",
            workflow = "nac",
            fastqs = fastq_10xv3_by_sample$fastq_paths,
            out = glue::glue("{analysis_cache}/kb_out/{fastq_10xv3_by_sample$sample_id}"),
            overwrite = FALSE, verbose = TRUE
        ),
        resources = large,
        pattern = map(fastq_10xv3_by_sample),
        format = "file_fast"
    ),
    tar_target(
        kb_count_10xv2_tcc,
        kb_count(
            technology = "10XV2", tcc = TRUE,
            index = "/ref_genomes/cellranger/human/kb_cr_ref/nac_transcriptome.idx",
            t2g = "/ref_genomes/cellranger/human/kb_cr_ref/nac_t2g.txt",
            c1 = "/ref_genomes/cellranger/human/kb_cr_ref/nac_cdna_t2c.txt",
            c2 = "/ref_genomes/cellranger/human/kb_cr_ref/nac_nascent_t2c.txt",
            workflow = "nac",
            fastqs = fastq_10xv2_by_sample$fastq_paths,
            out = glue::glue("{analysis_cache}/kb_out/{fastq_10xv2_by_sample$sample_id}"),
            overwrite = FALSE, verbose = TRUE
        ),
        resources = large,
        pattern = map(fastq_10xv2_by_sample),
        format = "file_fast"
    ),
    tar_target(
        kb_count_tcc_alltech,
        c(kb_count_10xv3_tcc, kb_count_10xv2_tcc),
        deployment = "main"
    ),
    # Cellranger -----------------------------------------------------------------------
    # Make nf-core sample sheet
    tar_target(sample_sheet_10xv2, hpap_fastq_to_sample_sheet(fastq_10xv2), deployment = "main"),
    tar_target(sample_sheet_10xv3, hpap_fastq_to_sample_sheet(fastq_10xv3), deployment = "main"),
    # rename fastq files as per cellranger requirements
    tar_target(renamed_sample_sheet_10xv2, rename_sample_sheet(sample_sheet_10xv2), deployment = "main"),
    tar_target(renamed_sample_sheet_10xv3, rename_sample_sheet(sample_sheet_10xv3), deployment = "main"),
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
    tar_target(cellranger_run_folders_10xv3, list.files(gsub("multiqc/multiqc_report.html", "cellranger/count", nfcore_scrnaseq_multiqc_10xv3), pattern = "HPAP", full.names = TRUE), deployment = "main"),
    tar_target(cellranger_run_folders_10xv2, list.files(gsub("multiqc/multiqc_report.html", "cellranger/count", nfcore_scrnaseq_multiqc_10xv2), pattern = "HPAP", full.names = TRUE), deployment = "main"),
    tar_target(cellranger_run_folders, c(cellranger_run_folders_10xv3, cellranger_run_folders_10xv2), deployment = "main"),
    # Run STAR solo on cellranger BAMs -----------------------------------------------------
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
    # Download SNPs ---------------------------------------------------------------------
    tar_target(
        snp_vcf,
        {
            destfile <- glue::glue("{analysis_cache}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz")
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
    # CellSNP-lite ---------------------------------------------------------------------
    tar_target(
        cellsnp_lite,
        run_cellsnp_lite(
            cellranger_run_folder = cellranger_run_folders,
            region_vcf = "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
        ),
        deployment = "main",
        pattern = map(cellranger_run_folders)
    ),
    # parse sample genotypes
    tar_target(protected_cohort, get_cellsnp_lite_genotypes("11\t2159843", "rs3842752"), deployment = "main"),
    tar_target(rs3842753_cohort, get_cellsnp_lite_genotypes("11\t2159830", "rs3842753"), deployment = "main"),
    tar_target(rs689_cohort, get_cellsnp_lite_genotypes("11\t2160994", "rs689"), deployment = "main"),
    tar_target(rs13266634_cohort, get_cellsnp_lite_genotypes("8\t117172544", "rs13266634"), deployment = "main"),
    # VarTrix ---------------------------------------------------------------------
    # prep reference snps
    tar_target(
        vartrix_vcf,
        filter_and_trim_vcf(
            vcf_path = snp_vcf,
            fai_path = "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai",
            output_vcf_path = "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.mod.vcf.gz"
        ),
        resources = small,
        format = "file_fast"
    ),
    tar_target(
        vartrix_coverage,
        run_vartrix(cellranger_run_folders, vartrix_vcf, mode = "coverage", mapq = 30),
        pattern = map(cellranger_run_folders),
        format = "file_fast",
    ),
    tar_target(
        vartrix_consensus,
        run_vartrix(cellranger_run_folders, vartrix_vcf, mode = "consensus", mapq = 30),
        pattern = map(cellranger_run_folders),
        format = "file_fast",
    ),
    # Add sample genotypes to metadata
    tar_target(
        pancdb_metadata_gt,
        merge_cellSNPlite_genotypes_with_sample_metadata(pancdb_metadata, protected_cohort, rs3842753_cohort, rs689_cohort),
        deployment = "main"
    ),
    # Cellbender using CellRanger counts ------------------------------------------------------------
    tar_target(cellbender_h5,
        run_cellbender(cellranger_run_folders),
        pattern = map(cellranger_run_folders),
        format = "file_fast",
    ),
    tar_target(
        cellbender_seurat_objects,
        make_seurat_cellbender(cellbender_h5, cellranger_run_folders, pancdb_metadata_gt),
        pattern = map(cellbender_h5, cellranger_run_folders),
        format = "file_fast",
        description = "Makes Seurat objects with CellBender and CellRanger data; adds cell QC metrics and renames cells using sample identifiers."
    ),
    tar_target(
        cellbender_qc_plots,
        seurat_plot_cellbender(cellbender_seurat_objects, "INS"),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = small
    ),
    # SingleR cell type annotation --------------------------------------------------------------
    # Get cell atlas from Tosti et al. 2021
    # https://doi.org/10.1053/j.gastro.2020.11.010
    # http://singlecell.charite.de/cellbrowser/pancreas/
    tar_target(
        tosti_etal_seurat_object,
        make_seurat_tosti_etal(),
        format = "file_fast", resources = xlarge
    ),
    tar_target(
        tosti_cell_type_csv,
        seurat_singleR_transfer_label(cellbender_seurat_objects, tosti_etal_seurat_object, cell_type_col = "Cluster"),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = small
    ),
    # Consolidated tosti cell types
    tar_target(
        tosti_consolidated_csv,
        {
            df <- tosti_cell_type_csv %>%
                readr::read_csv(show_col_types = FALSE, progress = FALSE) %>%
                mutate(cell_type = case_when(
                    pruned.labels %in% c("Alpha", "Beta", "Delta", "Gamma") ~ pruned.labels,
                    pruned.labels %in% c("Acinar-i", "Acinar-REG+", "Acinar-s") ~ "Acinar",
                    pruned.labels %in% c("Ductal", "MUC5B+ Ductal") ~ "Ductal",
                    TRUE ~ "Other"
                )) |>
                select(cell, cell_type)
            write.csv(df, glue::glue("{analysis_cache}/cell_type_out/tosti_consolidated.csv"), quote = FALSE)
            glue::glue("{analysis_cache}/cell_type_out/tosti_consolidated.csv")
        },
        format = "file_fast",
        resources = small
    ),
    # Azimuth --------------------------------------------------------
    tar_target(
        azimuth_reference_path,
        download_zenodo_files(
            urls = c(
                "https://zenodo.org/records/4546926/files/idx.annoy?download=1",
                "https://zenodo.org/records/4546926/files/ref.Rds?download=1"
            ),
            dest_dir = glue::glue("{analysis_cache}/data/azimuth")
        ),
        format = "file_fast",
        resources = tiny
    ),
    tar_target(
        azimuth_mapped_seurat_objects,
        seurat_azimuth(cellbender_seurat_objects, azimuth_reference_path),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = medium
    ),
    # Cell cycle annotation --------------------------------------------------------------
    tar_target(
        cell_cycle_csv,
        seurat_cell_cycle(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = small
    ),
    # Doublet annotation --------------------------------------------------------------
    tar_target(
        scDblFinder_csv,
        seurat_scDblFinder(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file_fast",
        resources = medium
    ),
    # HPAP annotation ---------------------------------------------------------------------
    tar_target(hpap_annotation_csv, get_hpap_cell_metadata(), format = "file_fast", resources = medium),
    # Project into reference annotation - Tosti et al 2021
    # * performs very poorly, not used for now
    # tar_target(
    #     ref_mapped_seurat_objects,
    #     seurat_project_into_ref(cellbender_seurat_objects, tosti_etal_seurat_object, reduction_model = "umap_harmony"),
    #     pattern = map(cellbender_seurat_objects),
    #     format = "file_fast",
    #     resources = large
    # ),
    # ddqc -------------------------------------------------------------------------------
    tar_target(
        ddqc_seurat_objects,
        # TODO update this from Nadia's project
        seurat_ddqc(cellbender_seurat_objects, scDblFinder_csv),
        pattern = map(cellbender_seurat_objects, scDblFinder_csv),
        format = "file_fast",
        resources = small
    ),
    # Convert seurat_objects to BPcells matrices -------------------------------------------
    tar_target(
        ddqc_bpcells_all,
        seurat_to_bpcells(ddqc_seurat_objects),
        format = "file_fast",
        pattern = map(ddqc_seurat_objects)
    ),

    # Drop failed samples before merge/clustering/DEG ------------------------------------
    tar_target(ddqc_bpcells, grep(failed_qc_donor_ids, ddqc_bpcells_all, value = TRUE, invert = TRUE), format = "file_fast", iteration = "vector", deployment = "main"),
    # Sketch
    tar_target(
        seurat_sketch_750,
        seurat_sketch(ddqc_bpcells, 750),
        format = "file_fast",
        pattern = map(ddqc_bpcells)
    ),
    # Merge
    tar_target(
        merged_seurat_sketch_750,
        seurat_merge(seurat_sketch_750, "merged_sketch_750"),
        format = "file_fast",
        resources = large
    ),
    # Integrate
    tar_target(
        integrated_seurat_sketch_750,
        seurat_sketch_harmony(merged_seurat_sketch_750, batch = "batch", pancdb_metadata),
        format = "file_fast",
        resources = large
    ),
    # Clustering --------------------------------------------------------------------------
    tar_target(
        cluster_merged_sketch_csv,
        seurat_cluster_ari(merged_seurat_sketch_750, assay = "sketch"),
        format = "file_fast",
        resources = large_mem
    ),
    # choose cluster manually
    tar_target(
        cluster_merged_sketch_man_csv,
        {
            cluster_merged_sketch <- readr::read_csv(cluster_merged_sketch_csv, show_col_types = FALSE, progress = FALSE) |>
                mutate(manual_clusters = SCT_snn_res.0.6) |> select(cell, manual_clusters)
            write.csv(cluster_merged_sketch, stringr::str_replace(cluster_merged_sketch_csv, ".csv", "_res0.6.csv"), quote=FALSE)
            stringr::str_replace(cluster_merged_sketch_csv, ".csv", "_res0.6.csv")
        },
        format = "file_fast",
        resources = tiny
    ),
    # Clustering with CHOIR ----------------------------------------------------------------
    #* Not working - fails for both skethc and RNA assay*#
    # tar_target(
    #     CHOIR_merged_seurat_sketch_750,
    #     seurat_CHOIR(merged_seurat_sketch_750, assay = "RNA", batch = "batch",  pancdb_metadata[c("sample_name", "batch")]),
    #     format = "file_fast",
    #     resources = xlarge
    # ),
    # GPT cell type annotation --------------------------------------------------------------
    # TODO
    # tar_target(
    #     gpt_cell_type_merged_sketch_750_csv,
    #     seurat_gpt_cell_type(merged_seurat_sketch_750, "sketch", group.by = "seurat_clusters", clusters_table = cluster_merged_sketch_csv),
    #     format = "file_fast",
    #     resources = large
    # ),
    # Aggregate cell annotation -------------------------------------------------------
    tar_target(
        aggregated_cell_annot_csv,
        seurat_aggregate_cell_annot(ddqc_seurat_objects, tosti_cell_type_csv, tosti_consolidated_csv, scDblFinder_csv, cell_cycle_csv, hpap_annotation_csv, cluster_merged_sketch_man_csv),
        resources = tiny
    ),
    # Annotate sketch object --------------------------------------------------------------
        tar_target(
        seurat_object_annotated_sketch,
        make_annotated_seurat_object(integrated_seurat_sketch_750, aggregated_cell_annot_csv, pancdb_metadata_gt),
        resources = small,
        format = "file_fast"
    ),
    # Project sketch onto full dataset 
    tar_target(
        seurat_object_annotated_full,
        seurat_project_sketch(seurat_object_annotated_sketch),
        format = "file_fast",
        resources = large
    ),    
    # Reports -------------------------------------------------------------------------------
    
    tar_target(report_one,
        render_report(here::here("quarto/pancdb_metadata.qmd")),
        deployment = "main", format = "file_fast"
    ),
    # Housekeeping --------------------------------------------------------------------------
    # Update the mtime of all files in the cache
    tar_target(
        touch_cache,
        sapply(
            c(analysis_cache, list.files(analysis_cache, full.names = TRUE, recursive = TRUE)),
            function(f) Sys.setFileTime(f, Sys.time())
        ),
        deployment = "main",
        cue = tarchetypes::tar_cue_age(touch_cache, as.difftime(1, units = "weeks"))
    )
)
#! TODO 
# TODO annotate beta like cells 
## Common Considerations for Quality Control Filters for Single Cell RNA-seq Data
## https://10xgenomics.com/resources/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data
## Essentially we follow this guide with respect to pre-processing
# 1. Filtering cell barcodes by UMI counts
# 2. Filtering cells by number of features
# 3. Filtering cells by percent mitochondrial expression
#    seurat_remove_outliers()
# 4. Filtering cells by doublet detection using community tools
#    Seurat::HTOdemux(), {scDblFinder}        (potentially {DoubletFinder}, {DoubletDecon})
# 5. Identifying and removing empty droplets based on the expression profile
#    {DropletQC}
# 6. Removing ambient RNAs associated with barcodes
#    {SoupX}, cellbender
# TODO: genotype viz. canadian study (rs...53)
# TODO: genotype viz. INS promoter (rs...689)
# TODO: kallisto quantification
# TODO: maybe https://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples
# TODO: vartrix genotyping of donors
# TODO: vartrix genotyping of cells
# TODO: stress score (XPB1 un/spliced)
# TODO: https://github.com/horsedayday/eQTLsingle
# TODO: sccomp


# clustering with anti-correlation
# https://www.nature.com/articles/s41467-023-43406-9#code-availability

## Kmerator
# https://academic.oup.com/nargab/article/3/3/lqab058/6308460

## sccomp - differential cell type analysis
# https://bioconductor.org/packages/release/bioc/vignettes/sccomp/inst/doc/introduction.html

## GPTCelltype - cell type annotation
# https://github.com/Winnie09/GPTCelltype
#
# https://github.com/const-ae/lemur
# If you have collected a single-cell RNA-seq dataset with more than one condition, lemur predicts
# for each cell and gene how much the expression would change if the cell had been in the other condition

# https://github.com/zhanghao-njmu/SCP
# Companion to scCustomize

# Variable Gene Selection:
# https://github.com/RuzhangZhao/mixhvg

# Genomic Plots
# https://www.bioconductor.org/packages/release/bioc/vignettes/GenomicPlot/inst/doc/GenomicPlot_vignettes.html#Introduction


# Splicing graphs
# https://www.bioconductor.org/packages/release/bioc/vignettes/SplicingGraphs/inst/doc/SplicingGraphs.pdf

# Splicing analysis
# MARVEL
# ** https://www.bioconductor.org/packages/release/bioc/vignettes/SGSeq/inst/doc/SGSeq.html

# XBP1 transcripts:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6982609/
# Unspliced XBP1-207
# discussion of similar approach using TCGA/recount2
# https://support.bioconductor.org/p/126127/
