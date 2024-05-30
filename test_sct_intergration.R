tar_source()
library(Seurat)
library(dplyr)
library(stringr)

tar_load(sct_ddqc_bpcells)
tar_load(pancdb_metadata_gt)

nodm <- pancdb_metadata_gt[pancdb_metadata_gt$diabetes_status == "NODM",]
prot<- nodm[nodm$protected == TRUE,1]
susc <- nodm[nodm$protected == FALSE, 1]

prot <- sct_ddqc_bpcells[stringr::str_detect(sct_ddqc_bpcells, stringr::str_c(prot$donor_id, collapse = "|"))]
susc <- sct_ddqc_bpcells[stringr::str_detect(sct_ddqc_bpcells, stringr::str_c(susc$donor_id, collapse = "|"))]

prot_list <- lapply(prot, load_seurat)
susc_list <- lapply(susc, load_seurat)

## gpt https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette

# Step 1: Identify integration features across all objects
features <- SelectIntegrationFeatures(object.list = c(prot_list, susc_list), nfeatures = 3000)

# Step 2: Prepare for SCT integration
prot_list <- PrepSCTIntegration(object.list = prot_list, anchor.features = features)
susc_list <- PrepSCTIntegration(object.list = susc_list, anchor.features = features)

# Step 3: Find anchors within each group
# prot_anchors <- FindIntegrationAnchors(object.list = prot_list, normalization.method = "SCT", anchor.features = features)
# qs::qsave(prot_anchors, file = "prot_anchors.qs")
prot_anchors <- qs::qread("prot_anchors.qs")

# susc_anchors <- FindIntegrationAnchors(object.list = susc_list, normalization.method = "SCT", anchor.features = features)
#qs::qsave(susc_anchors, file = "susc_anchors.qs")
susc_anchors <- qs::qread("susc_anchors.qs")

# Step 4: Integrate within each group


prot_integrated <- IntegrateData(anchorset = prot_anchors, normalization.method = "SCT")
qs::qsave(prot_integrated, file = "prot_integrated.qs")
#prot_integrated <- qs::qread("prot_integrated.qs")

susc_integrated <- IntegrateData(anchorset = susc_anchors, normalization.method = "SCT")
qs::qsave(susc_integrated, file = "susc_integrated.qs")
#susc_integrated <- qs::qread("susc_integrated.qs")

# Step 5: Find anchors between the integrated groups
group_anchors <- FindIntegrationAnchors(object.list = list(prot_integrated, susc_integrated), normalization.method = "SCT", anchor.features = features)
qs::qsave(group_anchors, file = "group_anchors.qs")

# Step 6: Integrate the groups
combined_integrated <- IntegrateData(anchorset = group_anchors, normalization.method = "SCT")
qs::qsave(combined_integrated, file = "combined_integrated.qs")

# Proceed with downstream analysis
# Scale the integrated data
combined_integrated <- ScaleData(combined_integrated, verbose = FALSE)

# Run PCA
combined_integrated <- RunPCA(combined_integrated, verbose = FALSE)

# Find neighbors and clusters
combined_integrated <- FindNeighbors(combined_integrated, dims = 1:30)
combined_integrated <- FindClusters(combined_integrated, resolution = 0.5)

# Run UMAP for visualization
combined_integrated <- RunUMAP(combined_integrated, dims = 1:30)

# Plot UMAP
DimPlot(combined_integrated, reduction = "umap", group.by = "orig.ident")


# https://github.com/satijalab/seurat/issues/7542
tar_source()
library(Seurat) # make sure you are running SeuratV5
library(SeuratData)
library(tidyseurat)

tar_load(merged_seurat_sketch_750)
tar_load(pancdb_metadata_gt)
tar_load(cell_cycle_csv)
meta_data <- pancdb_metadata_gt |> select(donor_id, diabetes_status, protected, batch)

full_obj <- load_seurat(merged_seurat_sketch_750)
full_obj <- full_obj |> scCustomize::Add_Sample_Meta(meta_data, join_by_seurat =  "orig.ident", join_by_meta = "donor_id")

obj <- full_obj |> dplyr::filter(diabetes_status == "NODM")

obj <- subset(obj, nFeature_RNA > 1000)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

# run sctransform
obj <- SCTransform(obj, vst.flavor = "v2", vars.to.regress = c("batch"), verbose = TRUE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# one-liner to run Integration
obj <- IntegrateLayers(
    object = obj, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    assay = "SCT", verbose = FALSE
)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "harmony_clusters")

#
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p1 <- DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = c("batch", "protected"),
    combine = FALSE
)
patchwork::wrap_plots(p1)
ggsave("plot.png")



####
# merge NODM
tar_source()
nodm <- pancdb_metadata_gt[pancdb_metadata_gt$diabetes_status == "NODM", ]
nodm <- sct_ddqc_bpcells[stringr::str_detect(sct_ddqc_bpcells, stringr::str_c(nodm$donor_id, collapse = "|"))]

nodm_list <- lapply(nodm, load_seurat)

seurat_merge_sct(nodm_list, "NODM_sct_merged")

x <- load_seurat("/scratch/domeally/DCD.tienhoven_scRNAseq.2024/seurat_merged/NODM_sct_merged.qs")

library(tidyseurat)

seurat_object <- x |>
    dplyr::filter(diabetes_status == "NODM") |>
    dplyr::filter(cell_type == "Beta")

    Seurat::RunPCA() |>
    Seurat::FindNeighbors() |>
    Seurat::FindClusters(resolution = 0.5)