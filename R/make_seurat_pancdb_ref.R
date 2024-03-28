#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @return
#' @author Denis O'Meally
#' @export
make_seurat_pancdb_ref <- function(pancdb_metadata, protected_cohort, annotated_seurat_bp) {
    require(Seurat)
    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    #--------------------------------------------------------------------------------
    # Sample metadata
    # rs3842752 annotation
    pancdb_metadata$protected <- pancdb_metadata$donor_id %in% protected_cohort$sample_id
    # 10X libraries that passed QC
    pancdb_metadata <- pancdb_metadata |>
        dplyr::filter(str_detect(reagent_kit, "10X") & !str_detect(donor_id, failed_qc_donor_ids)) |>
        dplyr::mutate(batch = as.integer(as.factor(reagent_kit)))

    meta_data <- select(pancdb_metadata, c(donor_id, diabetes_status, protected, batch, sample_sex, sample_age))

    # function to sample within each diabetes_status and protection status
    sample_protected_unprotected <- function(data) {
        data %>%
            group_by(protected) %>%
            sample_n(size = 2) %>%
            ungroup()
    }
    # Take a subset of each disease status and protected/unprotected status
    set.seed(42) # For reproducibility
    ref_samples <- meta_data %>%
        group_by(diabetes_status) %>%
        group_modify(~ sample_protected_unprotected(.x)) %>%
        ungroup()

    #--------------------------------------------------------------------------------
    # Seurat object
    seurat_object <- annotated_seurat_bp |> load_seurat()
    Idents(seurat_object) <- "orig.ident"
    seurat_object <- subset(x = seurat_object, idents = ref_samples$donor_id)
    # Run PCA
    seurat_object <- NormalizeData(seurat_object) |>
        FindVariableFeatures() |>
        ScaleData() |>
        RunPCA()
    # Integrate across donors
    seurat_object <- IntegrateLayers(
        object = seurat_object,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        verbose = FALSE
    )
    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", return.model = TRUE)
    seurat_object <- JoinLayers(seurat_object)
    #--------------------------------------------------------------------------------
    p1 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "tosti_cell_type",
        cols = custom_palette,
        label = TRUE, shuffle = TRUE, alpha = 0.4,
        repel = TRUE, label.size = 4
    ) & NoLegend() & labs(title = "Tosti cell types") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    p2 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "sample_sex",
        cols = c("#F8766D", "#00BFC4"),
        label = FALSE, shuffle = TRUE, alpha = 0.4,
    ) & labs(title = "Sex") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.8, 0.9))
    p3 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "batch",
        cols = c("#2ECC40", "#FF851B", "#5acdfa"),
        label = FALSE, shuffle = TRUE, alpha = 0.4,
    ) & labs(title = "Batch") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.8, 0.9))
    p4 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "cell_cycle_phase",
        cols = c("#0c6416", "#B10DC9", "#FF851B"),
        label = FALSE, shuffle = TRUE, alpha = 0.4,
    ) & labs(title = "Cell cycle") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.8, 0.9))
    p5 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "protected",
        cols = c("TRUE" = "red", "FALSE" = "black"),
        label = FALSE, shuffle = TRUE, alpha = 0.4,
    ) & labs(title = "Protected") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.8, 0.9))
    p6 <- DimPlot(
        seurat_object,
        reduction = "umap_harmony",
        group.by = "diabetes_status",
        cols = diabetes_palette,
        label = FALSE, shuffle = TRUE, alpha = 0.4,
    ) & labs(title = "Disease") & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.8, 0.9))
    umap_harmony <- patchwork::wrap_plots((p1 | p2) / (p3 | p4) / (p5 | p6))

    azimuth_vp <- scCustomize::Stacked_VlnPlot(seurat_object, features = cell_type_markers, x_lab_rotate = TRUE, colors_use = custom_palette, group.by = "azimuth_label")
    tosti_vp   <- scCustomize::Stacked_VlnPlot(seurat_object, features = cell_type_markers, x_lab_rotate = TRUE, colors_use = custom_palette, group.by = "tosti_cell_type")

    #--------------------------------------------------------------------------------
    # Save results
    seurat_object_path <- (glue::glue("{analysis_cache}/pancdb_ref/pancdb_harmony_ref.qs"))
    dir.create(dirname(seurat_object_path), showWarnings = FALSE, recursive = TRUE)

    ggsave(sub("\\.qs", "\\.png", seurat_object_path), umap_harmony, width = 12, height = 20, dpi = 300)
    ggsave(sub("\\.qs", "\\_tosti_vlnPlt.png", seurat_object_path), tosti_vp, width = 10, height = 10, dpi = 300)
    ggsave(sub("\\.qs", "\\_azimuth_vlnPlt.png", seurat_object_path), azimuth_vp, width = 10, height = 10, dpi = 300)

    qs::qsave(seurat_object, seurat_object_path)

    return(seurat_object_path)
}
