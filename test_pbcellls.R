suppressPackageStartupMessages({
library(targets)
library(Seurat)
library(scCustomize)
library(BPCells)
library(progress)
library(tidyverse)
})
hprcc::init_multisession()
tar_source()
tar_load(annotated_seurat_bp)
#seurat_object<-"/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/merged_seurat_bp.rds"|>load_seurat()
#seurat_object<-annotated_seurat_bp |> load_seurat()
seurat_object<-annotated_seurat_bp |> LoadSeuratRds()
object.size(merged_seurat) |> print(units = "GB")

object <- JoinLayers(merged_seurat)
object <- NormalizeData(seurat_object)
# split assay into 77 layers
object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
object <- FindVariableFeatures(object, verbose = FALSE)

#  Sample representative cells from each dataset
object <- SketchData(object = object, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
object
# Save seurat object
saveRDS(object, file = '/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/data/merged_bpcells_sketch.rds')

# Perform integration on the sketched cells across samples
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)
# integrate the datasets
object <- IntegrateLayers(
    object,
    method = RPCAIntegration,
    orig = "pca",
    new.reduction = "integrated.rpca",
    dims = 1:30, k.anchor = 20, verbose = F)

# cluster the integrated data
object <- FindNeighbors(object, reduction = "integrated.rpca", dims = 1:30)
object <- FindClusters(object, resolution = 2)

object <- RunUMAP(object, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)

# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
object[["sketch"]] <- JoinLayers(object[["sketch"]])
c10_markers <- FindMarkers(object = object, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
head(c10_markers)

# You can now annotate clusters using marker genes.  We performed this step, and include the
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(object, group.by = "orig.ident", reduction = "umap")
plot.s2 <- DimPlot(object, group.by = "azimuth_label", reduction = "umap")
# Compare healthy and diabetic samples

# Save seurat object
saveRDS(object, file = '/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/data/study_cohort_merged_bpcells_integrated.rds')


# Project full datasets

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
object[["sketch"]] <- split(object[["sketch"]], f = object$orig.ident)

object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


object <- ProjectData(object = object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
    full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "azimuth_label"))

object <- RunUMAP(object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAPfull_")

# Save seurat object
saveRDS(object, file = '/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/data/study_cohort_merged_bpcells_integrated_full.rds')

# Compare healthy and diabetic samples
bulk <- AggregateExpression(object, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("azimuth_label",
    "orig.ident", "protected"))
# each sample is an individual-specific celltype-specific pseudobulk profile
tail(Cells(bulk))

library(ggrepel)
  beta.bulk <- subset(bulk, azimuth_label == "beta")
  Idents(beta.bulk) <- "protected"
  de_markers <- FindMarkers(beta.bulk, ident.1 = "TRUE", ident.2 = "FALSE", slot = "counts", test.use = "DESeq2",
      verbose = F)
  de_markers$gene <- rownames(de_markers)
  ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
      ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
      "")), colour = "red", size = 3)
ggsave("volcano.png")


## Load merged object
hprcc::init_multisession()
tar_source()
seurat_object<-"/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/merged_seurat_bp.rds"|>load_seurat()
object.size(seurat_object) |> print(units = "GB")
joined_obj <- JoinLayers(seurat_object)

library(Seurat)
library(BPCells)

t1_CreateSketchAssay <- system.time({
    obj <- NormalizeData(seurat_object)
    obj <- FindVariableFeatures(obj, layer = "counts")
    obj <- SketchData(
        object = obj,
        ncells = 750,
        method = "LeverageScore",
        sketched.assay = "sketch"
        )

})

joined_obj <- JoinLayers(obj)
joined_obj <- SketchData(
        object = joined_obj,
        ncells = 5000,
        method = "LeverageScore",
        sketched.assay = "sketch"
        )



##sketched clustering
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object)
object <- FindNeighbors(object, dims = 1:50)
object <- FindClusters(object, resolution = 2)
object <- RunUMAP(object, dims = 1:50, return.model = TRUE)
DimPlot(object, reduction = "umap") + NoLegend()
ggplot("plot.png")



###_----------#####
suppressPackageStartupMessages({
        require(tidyverse)
        require(Seurat)
        require(scCustomize)
        require(BPCells)
        require(progress)
    })

tar_load(c(ddqc_seurat_objects, pancdb_metadata, protected_cohort, azimuth_mapped_seurat_objects, cell_cycle_csv, tosti_cell_type_csv))
tar_source()
# Sample metadata
    # rs3842752 annotation
    pancdb_metadata$protected <- pancdb_metadata$donor_id %in% protected_cohort$sample_id
    # 10X libraries that passed QC
    pancdb_metadata <- pancdb_metadata |>
        dplyr::filter(str_detect(reagent_kit, "10X") & !str_detect(donor_id, failed_qc_donor_ids)) |>
        dplyr::mutate(batch = as.integer(as.factor(reagent_kit)))
    ddqc_seurat_object_paths <- data.frame(
        "qs_path" = ddqc_seurat_objects,
        "donor_id" = gsub(".*/(HPAP-\\d+)_ddqc\\.qs", "\\1", ddqc_seurat_objects)
    )
    pancdb_metadata <- pancdb_metadata |> left_join(ddqc_seurat_object_paths, by = "donor_id")

    seurat_paths <- pancdb_metadata |> pull(qs_path)
    meta_data <- select(pancdb_metadata, c(donor_id, protected, batch, sample_sex, sample_age))

    # Initialise progress bar
    pb <- progress_bar$new(
        format = "  Loading [:bar] :percent eta: :eta (Loading :current of :total)",
        total = length(seurat_paths), clear = FALSE, width = 60
    )

    data_list <- c()
    metadata_list <- c()

    for (path in seurat_paths) {
        seurat_object <- load_seurat(path)
        # add sample metadata
        seurat_object <- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, join_by_seurat = "orig.ident", join_by_meta = "donor_id")
        sample_id <- seurat_object$orig.ident[1]
        # add cell metadata 
        azimuth_data <- sub("\\.qs", "\\.csv", grep(sample_id, azimuth_mapped_seurat_objects, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        cell_cycle_data <- sub("\\.qs", "\\.csv", grep(sample_id, cell_cycle_csv, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        tosti_cell_type_data <- sub("\\.qs", "\\.csv", grep(sample_id, tosti_cell_type_csv, value = TRUE)) |> readr::read_csv(show_col_types = FALSE, progress = FALSE)
        cell_metadata <- azimuth_data |>
            left_join(cell_cycle_data, by = "cell") |>
            left_join(tosti_cell_type_data, by = "cell") |>
            select(c(cell, azimuth_label = "predicted.annotation.l1", cell_cycle_phase = "Phase", tosti_cell_type = "labels")) |>
            column_to_rownames(var = "cell")
        seurat_object <- Seurat::AddMetaData(seurat_object, cell_metadata)

        bp_dir <- glue::glue("{analysis_cache}/bpcells_out/{sample_id}")

        iterable_matrix <- as(seurat_object[["RNA"]]$counts, "IterableMatrix")
        # Convert counts to BPCells / load the BPCells folder
        if (!dir.exists(bp_dir)) {
            write_matrix_dir(mat = iterable_matrix, dir = bp_dir, compress = TRUE)
        }
        mat <- open_matrix_dir(dir = bp_dir)

        # Get mat and metadata for later merging
        data_list[[sample_id]] <- mat
        metadata_list[[sample_id]] <- seurat_object[[]]
        # Increment progress bar after each iteration
        pb$tick()
    }

    # Build the merged seurat object
    metadata <- Reduce(rbind, metadata_list)

    merged_seurat <- Seurat::CreateSeuratObject(counts = data_list, meta.data = metadata)
## ok now sketch
object.size(merged_seurat) |> print(units = "GB")

joined_object <- JoinLayers(merged_seurat)
object <- NormalizeData(joined_object)
# split assay into 3 layers
#split_object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
object <- FindVariableFeatures(object, verbose = TRUE)

#  Sample representative cells from each dataset
sketch_object <- SketchData(object = object, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

sketch_object<-tar_read(merged_seurat_bp_sketch) |> load_seurat()
sketch_object<-"/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/merged_seurat_bp_sketch_200cells.rds" |> load_seurat()

joined_obj <- JoinLayers(sketch_object)
t2_SketchClustering <- system.time({
    DefaultAssay(joined_obj) <- "sketch"
    obj <- FindVariableFeatures(joined_obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    obj <- FindNeighbors(obj, dims = 1:50)
    obj <- FindClusters(obj)
})
obj <- RunUMAP(obj, dims = 1:50, return.model = TRUE)
DimPlot(obj, label = TRUE, group.by = "tosti_cell_type", cols= custom_palette, reduction = "umap") + NoLegend()
ggsave("plot200.png")
