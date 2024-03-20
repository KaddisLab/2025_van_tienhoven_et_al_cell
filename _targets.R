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
    #tar_target(fastq_md5sums, read_md5sums(here::here(glue::glue("{analysis_cache}/data/hpapdata"))), resources = small),
    #TODO vectorise this function and make groups of 20 or more files to avoid overhead of calling md5sum for each file
    #TODO look at tar_rep() for this...
    #tar_target(md5sums_check, check_md5sum(fastq_md5sums), map(fastq_md5sums), resources = tiny),
    # import metadata --------------------------------------------------------
    tar_target(pancdb_metadata, get_pancdb_metadata(), deployment = "main"),
    # get fastq files --------------------------------------------------------
    tar_target(
        fastq_10x,
        grep("SSq2", list.files(path = glue::glue("{analysis_cache}/data/hpapdata"), pattern = "R[12](_001)?(_fastq-data)?\\.fastq\\.gz$", full.names = TRUE,  recursive = TRUE),
            value = TRUE, invert = TRUE, perl = TRUE),
        deployment = "main"),
    tar_target(
        fastq_10xv2,
        fastq_10x[grep(pancdb_metadata$donor_id[pancdb_metadata$reagent_kit == "10X-Chromium-GEX-3p-v2"] |> na.omit() |> paste(collapse="|"), fastq_10x)],
        deployment = "main"),
    tar_target(
        fastq_10xv3,
        fastq_10x[grep(pancdb_metadata$donor_id[pancdb_metadata$reagent_kit %in% c("10X-Chromium-GEX-3p-v3", "10X-Chromium-GEX-3p-v3.1")] |> na.omit() |> paste(collapse="|"), fastq_10x)],
        deployment = "main"),
    tar_target(
        fastq_ss2,
        list.files(path = glue::glue("{analysis_cache}/data/hpapdata"), pattern = "_scRNA_SSq2.*\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        deployment = "main"),
    tar_target(
        fastq_ss3,
        list.files(path = glue::glue("{analysis_cache}/data/hpapdata"), pattern = "_scRNA_\\d{5}_.*\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        deployment = "main"),

# kallisto-bustools ---------------------------------------------------------------------
    # Group fastqs by sample
    tar_target(fastq_10xv3_by_sample, arrange_by_sample_10x(fastq_10xv3), deployment = "main"),
    tar_target(fastq_10xv2_by_sample, arrange_by_sample_10x(fastq_10xv2), deployment = "main"),
    tar_target(fastq_ss2_by_sample, arrange_by_sample_ssq(fastq_ss2), deployment = "main"),
    tar_target(fastq_ss3_by_sample, arrange_by_sample_ssq(fastq_ss3), deployment = "main"),
    # fetch kallisto reference
    tar_target(
        kb_ref_hg38_std,
        kb_ref(species = "human", index = glue::glue("{analysis_cache}/data/kb_ref/hg38_std.idx"), t2g = glue::glue("{analysis_cache}/data/kb_ref/hg38_std_t2g.txt"), overwrite = TRUE),
        format = "file", resources = tiny
    ),
    # fetch kallisto reference - nac
    tar_target(
        kb_ref_hg38_nac,
        kb_ref(
            species = "human",  workflow = "nac", index = glue::glue("{analysis_cache}/data/kb_ref/hg38_nac.idx"), t2g = glue::glue("{analysis_cache}/data/kb_ref/hg38_nac_t2g.txt"),
            c1 = glue::glue("{analysis_cache}/data/kb_ref/hg38_nac_cdna.txt"), c2 = glue::glue("{analysis_cache}/data/kb_ref/hg38_nac_nascent.txt"), overwrite = TRUE),
        format = "file", resources = tiny
    ),
    # kb_count
    tar_target(
        kb_count_10xv3,
        kb_count(technology = "10XV3",
            index = kb_ref_hg38_std,
            t2g = gsub(".idx", "_t2g.txt", kb_ref_hg38_std),
            fastqs = fastq_10xv3_by_sample$fastq_paths,
            out = glue::glue("{analysis_cache}/kb_out/{fastq_10xv3_by_sample$sample_id}"),
            overwrite = FALSE),
        resources = medium,
        pattern = map(fastq_10xv3_by_sample),
        format = "file"
        ),
    tar_target(
        kb_count_10xv2,
        kb_count(technology = "10XV2",
            index = kb_ref_hg38_std,
            t2g = gsub(".idx", "_t2g.txt", kb_ref_hg38_std),
            fastqs = fastq_10xv2_by_sample$fastq_paths,
            out = glue::glue("{analysis_cache}/kb_out/{fastq_10xv2_by_sample$sample_id}"),
            overwrite = FALSE),
        resources = medium,
        pattern = map(fastq_10xv2_by_sample),
        format = "file"
        ),
    tar_target(
        kb_count_alltech,
        c(kb_count_10xv3, kb_count_10xv2),
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
            protocol = "10XV2"),
        deployment = "main",
        ),
    tar_target(
        nfcore_scrnaseq_multiqc_10xv3,
        run_nfcore_scrnaseq(
            run_folder = "cellranger_out/10xv3",
            sample_sheet = renamed_sample_sheet_10xv3,
            protocol = "10XV3"),
        deployment = "main"
        ),
    tar_target(cellranger_run_folders_10xv3, list.files(gsub("multiqc/multiqc_report.html", "cellranger/count", nfcore_scrnaseq_multiqc_10xv3), pattern = "HPAP", full.names = TRUE), deployment = "main"),
    tar_target(cellranger_run_folders_10xv2, list.files(gsub("multiqc/multiqc_report.html", "cellranger/count", nfcore_scrnaseq_multiqc_10xv2), pattern = "HPAP", full.names = TRUE), deployment = "main"),
    tar_target(cellranger_run_folders, c(cellranger_run_folders_10xv3, cellranger_run_folders_10xv2), deployment = "main"),
    # Failed QC cohort
    tar_target(
        cellranger_run_folders_failed_qc,
        grep(paste0(
            c("HPAP-021|HPAP-023|HPAP-027|"), # MultiQC v2
            c("HPAP-038|HPAP-093"), # MultiQC v3 
            collapse="|"), 
            cellranger_run_folders, value = TRUE),
        deployment = "main"),
    # Non-diabetic cohort
    tar_target(
        cellranger_run_folders_nodx,
        grep(pancdb_metadata$donor_id[pancdb_metadata$diabetes_status == "NODM" ] |> na.omit() |> paste(collapse="|"), cellranger_run_folders, value = TRUE) |>
            setdiff(cellranger_run_folders_failed_qc),
        deployment = "main"),
# Download SNPs ---------------------------------------------------------------------
    tar_target(
        snp_vcf,
        {destfile <- glue::glue("{analysis_cache}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz")
            download.file(
                url = "http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz",
                destfile = destfile,
                mode = "wb")
            return(destfile)},
        deployment = "main",
        format = "file"
    ),
    # CellSNP-lite ---------------------------------------------------------------------
    tar_target(
        cellsnp_lite,
        run_cellsnp_lite(
            cellranger_run_folder = cellranger_run_folders,
            region_vcf = "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"),
        pattern = map(cellranger_run_folders),
        deployment = "main"
    ),
    # parse sample genotypes
    tar_target(protected_cohort, get_cellsnp_lite_genotypes("11\t2159843", "rs3842752"), deployment = "main"),
    tar_target(rs3842753_cohort, get_cellsnp_lite_genotypes("11\t2159830", "rs3842753"), deployment = "main"),
    tar_target(rs13266634_cohort, get_cellsnp_lite_genotypes("8\t117172544", "rs13266634"), deployment = "main"),
    # VarTrix ---------------------------------------------------------------------
    # prep reference snps
    tar_target(
        vartrix_vcf,
        filter_and_trim_vcf(
            vcf_path = snp_vcf,
            fai_path = "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai",
            output_vcf_path = "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.mod.vcf.gz"),
        resources = small,
        format = "file"
        ),
    tar_target(
        vartrix_coverage,
        run_vartrix(cellranger_run_folders, vartrix_vcf, mode = "coverage", mapq = 30),
        deployment = "main",
        pattern = map(cellranger_run_folders),
        format = "file",
        ),
    tar_target(
        vartrix_consensus,
        run_vartrix(cellranger_run_folders, vartrix_vcf, mode = "consensus", mapq = 30),
        deployment = "main",
        pattern = map(cellranger_run_folders),
        format = "file",
        ),
# Cellbender using CellRanger counts ------------------------------------------------------------
    tar_target(cellbender_h5,
        run_cellbender(cellranger_run_folder = cellranger_run_folders),
        pattern = map(cellranger_run_folders),
        deployment = "main",
        format = "file",
    ),
    tar_target(
        cellbender_seurat_objects,
        make_seurat_cellbender(cellbender_h5, cellranger_run_folders),
        pattern = map(cellbender_h5, cellranger_run_folders),
        deployment = "main", # must run on "main" for some reason??
        format = "file"
    ),
    tar_target(
        cellbender_qc_plots,
        seurat_plot_cellbender(cellbender_seurat_objects, "INS"),
        pattern = map(cellbender_seurat_objects),
        format = "file",
        resources = small
    ),
# SingleR cell type annotation --------------------------------------------------------------
    # Get cell atlas from Tosti et al. 2021
    # https://doi.org/10.1053/j.gastro.2020.11.010
    # http://singlecell.charite.de/cellbrowser/pancreas/
    tar_target(
        tosti_etal_seurat_object,
        make_seurat_tosti_etal(),
        format = "file", resources = xlarge
    ),
    tar_target(
        tosti_cell_type_csv,
        seurat_singleR_transfer_label(cellbender_seurat_objects, tosti_etal_seurat_object, cell_type_col = "Cluster"),
        pattern = map(cellbender_seurat_objects),
        format = "file",
        resources = small
    ),
    # Azimuth - human pancreas --------------------------------------------------------
    tar_target(
        azimuth_reference_path,
        download_zenodo_files(
            urls = c(
                "https://zenodo.org/records/4546926/files/idx.annoy?download=1",
                "https://zenodo.org/records/4546926/files/ref.Rds?download=1"),
            dest_dir = glue::glue("{analysis_cache}/data/azimuth")),
        format = "file",
        resources = tiny
    ),
    tar_target(
        azimuth_mapped_seurat_objects,
        seurat_azimuth(cellbender_seurat_objects, azimuth_reference_path),
        pattern = map(cellbender_seurat_objects),
        format = "file",
        resources = medium
    ),
# Cell cycle annotation --------------------------------------------------------------
    tar_target(
        cell_cycle_csv,
        seurat_cell_cycle(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file",
        resources = small
    ),
# Doublet annotation --------------------------------------------------------------
    tar_target(
        scDblFinder_csv,
        seurat_scDblFinder(cellbender_seurat_objects),
        pattern = map(cellbender_seurat_objects),
        format = "file",
        resources = medium
    ),
# Seurat reference annotation - Tosti et al 2021
#! performs very poorly, not used for now
    # tar_target(
    #     ref_mapped_seurat_objects,
    #     seurat_project_into_ref(cellbender_seurat_objects, tosti_etal_seurat_object, reduction_model = "umap_harmony"),
    #     pattern = map(cellbender_seurat_objects),
    #     format = "file",
    #     resources = large
    # ),
# ddqc -------------------------------------------------------------------------------
    tar_target(
        ddqc_seurat_objects,
        seurat_ddqc(cellbender_seurat_objects, scDblFinder_csv),
        pattern = map(cellbender_seurat_objects, scDblFinder_csv),
        format = "file",
        resources = small
    ),
# Clustering --------------------------------------------------------------------------
#! Not needed for now
    # tar_target(
    #     seurat_cluster_ari_csv,
    #     seurat_cluster_ari(ddqc_seurat_objects),
    #     pattern = map(ddqc_seurat_objects),
    #     format = "file",
    #     resources = medium
    # ),
# Housekeeping --------------------------------------------------------------------------
    # Update the mtime of all files in the cache
    tar_target(
        touch_cache,
                sapply(c(analysis_cache, list.files(analysis_cache, full.names = TRUE, recursive = TRUE)), 
                    function(f) Sys.setFileTime(f, Sys.time())),
        deployment = "main",
        cue = tarchetypes::tar_cue_age(touch_cache, as.difftime(3, units = "weeks"))
    )

)
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

    # demultiplexing
    # get the path of the raw matrix file for sample mixes
    # tar_files_input(raw_matrix_paths, get_raw_matrix_path(readr::read_csv(glue::glue("{analysis_cache}/data/sample_metadata.csv"), show_col_types = FALSE)), format = "file"),

    # # make a seurat object for each sample mix, including CMO tags and sample mapping
    # tar_target(
    #     seurat_cellplex_objects,
    #     create_cellplex_seurat_object(raw_matrix_paths, sample_metadata),
    #     pattern = map(raw_matrix_paths),
    #     iteration = "list"),

    # # perform  cell type annotation
    # tar_target(seurat_celltype,
    #     seurat_singleR_mmu(seurat_cellplex_objects),
    #     pattern = map(seurat_cellplex_objects),
    #     iteration = "list",
    #     resources = medium),


# clustering with anti-correlation
# https://www.nature.com/articles/s41467-023-43406-9#code-availability

## Kmerator
# https://academic.oup.com/nargab/article/3/3/lqab058/6308460