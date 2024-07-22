make_pseudo_bulk_object <- function(seurat_object) {

    suppressPackageStartupMessages({
        require(tidyseurat)
        require(Seurat)
        require(targets)
        require(Matrix)
        require(tidybulk)
        require(SummarizedExperiment)
    })

# tar_load(seurat_object_lognorm_annotated)
# seurat_object <- load_seurat(seurat_object_lognorm_annotated)

# seurat_object <- seurat_object |>
#     dplyr::filter(diabetes_status == "NODM") 

seurat_object$rs689_consensus <- relevel(factor(seurat_object$rs689_consensus), ref = "TT")
seurat_object$rs3842752_consensus <- relevel(factor(seurat_object$rs3842752_consensus), ref = "GG")
seurat_object$rs3842753_consensus <- relevel(factor(seurat_object$rs3842753_consensus), ref = "TG")
seurat_object$tissue_source <- factor(seurat_object$tissue_source)
seurat_object$technology <- factor(seurat_object$technology)

# Remove undetected genes
seurat_object <- subset(seurat_object, features = rownames(seurat_object)[Matrix::rowSums(seurat_object[["RNA"]]$counts > 0) > 0])

# Calculate QC metrics and remove cells with few or many detected genes
seurat_object$nFeature_RNA <- Matrix::colSums(seurat_object[["RNA"]]$counts > 0)
nFeature_RNA_median <- median(seurat_object$nFeature_RNA)
nFeature_RNA_mad <- mad(seurat_object$nFeature_RNA)
lower_bound <- nFeature_RNA_median - 2 * nFeature_RNA_mad
upper_bound <- nFeature_RNA_median + 2 * nFeature_RNA_mad
cells_to_keep <- seurat_object$nFeature_RNA >= lower_bound & seurat_object$nFeature_RNA <= upper_bound
seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])

# Remove lowly expressed genes
seurat_object <- subset(seurat_object, features = rownames(seurat_object)[Matrix::rowSums(seurat_object[["RNA"]]$counts > 1) >= 10])

# Pseudoulk
pseudo_bulk <-
    seurat_object |>
    aggregate_cells(c(orig.ident, cell_type), assays = "RNA") 
}
