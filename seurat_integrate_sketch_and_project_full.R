seurat_integrate_sketch_and_project_full <- function(
    seurat_object_sketch,
    res = seq(0.1, 1.5, by = 0.05),
    dims = 30,
    regress_out = NULL,
    do.plot = TRUE, ...) {
    set.seed(42)

    #seurat_object <- load_seurat(seurat_object_sketch)
    tar_source()
    seurat_object<-"/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/data/merged_seurat_bp_sketch.rds"|>load_seurat()
    ref_seurat_object <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/pancdb_ref/pancdb_harmony_ref.qs" |> load_seurat() |> Seurat::NormalizeData(verbose = FALSE)
    format(object.size(seurat_object), units = "Gb")
    res = 2
    dims = 30
    options(parallelly.availableCores.methods = "BiocParallel")
    hprcc::init_multisession()
    future::plan("multisession", workers = 4)


    # Inspired by Sketch integration using a 1 million cell dataset from Parse Biosciences -------------------
    # https://satijalab.org/seurat/articles/mousebrain_sketch_clustering
    # https://github.com/satijalab/seurat/issues/6048
    # Perform integration on the sketched cells across samples
    DefaultAssay(seurat_object) <- "sketch"
        seurat_object <- FindVariableFeatures(seurat_object, verbose = FALSE)
        seurat_object <- ScaleData(seurat_object, verbose = FALSE)
        seurat_object <- RunPCA(seurat_object, verbose = FALSE)
        # integrate the datasets - takes about 1 hour
        seurat_object <- IntegrateLayers(seurat_object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
            dims = 1:dims, k.anchor = 20, verbose = TRUE)
        # cluster the integrated data
        seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.rpca", dims = 1:dims)
        seurat_object <- FindClusters(seurat_object, resolution = res)
        seurat_object <- RunUMAP(seurat_object, reduction = "integrated.rpca", dims = 1:dims, return.model = TRUE, verbose = FALSE)

        # rejoin the layers in the sketched assay this is required to perform differential expression
        joined_seurat_object<-seurat_object
        joined_seurat_object[["sketch"]] <- JoinLayers(seurat_object[["sketch"]])
        Idents(joined_seurat_object) <- "tosti_cell_type"
        beta_markers <- FindMarkers(object = joined_seurat_object, ident.1 = "Beta")
        beta_markers|>filter(p_val_adj < 0.05)|>head(30)
        beta_markers_protected<-FindMarkers(joined_seurat_object, ident.1 = "TRUE", group.by = 'protected', subset.ident = "Beta")
        beta_markers_protected|>filter(p_val_adj < 0.05) |> arrange(desc(abs(avg_log2FC)))|>head(30   )

        # Integrate the full datasets
        seurat_object <- ProjectIntegration(object = seurat_object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
        seurat_object <- ProjectData(
            object = seurat_object,
            sketched.assay = "sketch",
            assay = "RNA",
            sketched.reduction = "integrated.rpca.full",
            full.reduction = "integrated.rpca.full",
            dims = 1:30,
            refdata = list(
                azimuth_labels = "azimuth_label",
                tosti_cell_type = "tosti_cell_type"
            ))

        seurat_object <- RunUMAP(
            seurat_object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
            reduction.key = "UMAPfull_")
        DimPlot(seurat_object, group.by = "tosti_cell_type", reduction = "umap.full", cols = custom_palette, label = TRUE, label.size = 4, repel = TRUE) +
            NoLegend() + theme_void() + theme(legend.position = "none")
        ggsave("plot.png")
#qs::qsave(seurat_object, file = "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/merged_seurat_bp_sketch_integrated_20240327.rds")

######################################################################



# Compare healthy and diabetic samples
bulk <- AggregateExpression(seurat_object, return.seurat = TRUE, slot = "counts", assays = "RNA", group.by = c("azimuth_label",
    "orig.ident", "protected", "disease_status"))
# each sample is an individual-specific celltype-specific pseudobulk profile
tail(Cells(bulk))

library(ggrepel)
  beta.bulk <- subset(bulk, azimuth_label == "beta")
  Idents(beta.bulk) <- "protected"
  de_markers <- FindMarkers(beta.bulk, ident.1 = "TRUE", ident.2 = "FALSE", slot = "counts", test.use = "wilcox", verbose = TRUE)
  de_markers$gene <- rownames(de_markers)
  ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
      ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
      "")), colour = "red", size = 3)
ggsave("volcano.png")







        # Compare healthy and diabetic samples
        pseudobulk <- AggregateExpression(
            seuratseurat__object,
            return.seurat = TRUE, slot = "counts", assays = "RNA",
            group.by = c("azimuth_labels", "orig.ident", "disease_status"))
# each sample is an individual-specific celltype-specific pseudobulk profile
tail(Cells(pseudobulk))
          beta.bulk <- subset(bulk, celltype.full == "beta")
  Idents(beta.bulk) <- "disease_status"
  de_markers <- FindMarkers(beta.bulk, ident.1 = "TRUE", ident.2 = "FALSE", slot = "counts", test.use = "DESeq2", verbose = F)
  de_markers$gene <- rownames(de_markers)
  ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) +
    geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
    ylab("-log10(unadjusted p-value)") +
    geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene, "")), colour = "red", size = 3)

    return(seurat_object_path)
}
