# Load necessary libraries
library(scGate)
library(Seurat)

# Define scGate models for each cell type
models <- list(
    Acinar = gating_model(name = "Acinar", signature = c("PNLIP", "CPA1", "CPB1", "PRSS2", "KLK1", "CUZD1", "AMY2A", "CTRB2", "CTRL", "CLPS")),
    Alpha = gating_model(name = "Alpha", signature = c("GCG", "CRH", "SPOCK3", "F10", "IRX2", "POPDC3", "C5orf38", "GC", "SPINK4", "TMEM236")),
    Beta = gating_model(name = "Beta", signature = c("MAFA", "IAPP", "ADCYAP1", "INS", "AC132217.2", "LRRTM3", "HHATL", "SAMD11", "C1orf127", "WSCD2")),
    Delta = gating_model(name = "Delta", signature = c("SST", "GPC5-AS11", "LRFN5", "LY6H1", "BCHE", "CALB11", "LINC01571", "CBLN4", "CRYGD", "F5")),
    Ductal = gating_model(name = "Ductal", signature = c("FGFBP1", "CFTR", "MMP7", "KRT23", "S100A14", "KRT19", "LGALS4", "CEACAM7", "VTCN1", "TRPV6")),
    Endothelial = gating_model(name = "Endothelial", signature = c("PLVAP", "KDR", "ESM1", "VWF", "PCAT19", "PECAM1", "CLDN5", "CLEC14A", "ADGRL4", "FCN3")),
    Gamma_Epsilon = gating_model(name = "Gamma_Epsilon", signature = c("PPY", "CARTPT", "NPW", "GPC5-AS1", "CALB1", "SERTM1", "LY6H", "GHRL", "ERICH3", "PRG4")),
    Macrophage = gating_model(name = "Macrophage", signature = c("C1QA", "C1QB", "C1QC", "TYROBP", "SDS", "LAPTM5", "MS4A7", "RGS1", "ITGB2", "GPR183")),
    Mast = gating_model(name = "Mast", signature = c("TPSB2", "TPSAB1", "CD69", "SAMSN1", "CD37", "CD48", "KIT", "GCSAML", "SRGN", "ITGAX")),
    Stellate = gating_model(name = "Stellate", signature = c("COL3A1", "DCN", "LUM", "SFRP2", "PTGDS", "PRRX1", "COL1A2", "BGN", "C7", "COL6A3")),
    Cycling_alpha = gating_model(name = "Cycling_alpha", signature = c("UBE2C", "CENPF", "BIRC5", "TOP2A", "FOXM1", "NUSAP1", "TPX2", "CKAP2L", "MKI67", "CKS2")),
    Active_stellate = gating_model(name = "Active_stellate", signature = c("APOD", "SFRP2", "PTGDS", "CFD", "CCL11", "DPT", "C7", "PODN", "CHI3L1", "CYP1B1")),
    Quiescent_stellate = gating_model(name = "Quiescent_stellate", signature = c("FABP4", "RGS5", "PDK4", "ADIRF", "MCAM", "MUSTN1", "ESAM", "EFHD1", "ADAMTS9", "ADAMTS4")),
    `MUC5B+_ductal` = gating_model(name = "MUC5B_plus_ductal", signature = c("TFF2", "TFF1", "MGST1", "TFF3", "TCN1", "FXYD3", "SPINK1", "MMP1", "C19orf33", "LDHB"))
)

original <- list(
    Beta = gating_model(name = "Beta", signature = c("INS", "IAPP")),
    Alpha = gating_model(name = "Alpha", signature = c("GCG")),
    Proliferating_Cells = gating_model(name = "Proliferating_Cells", signature = c("MKI67", "CDK1")),
    Delta = gating_model(name = "Delta", signature = c("SST")),
    Gamma = gating_model(name = "Gamma", signature = c("PPY")),
    Epsilon = gating_model(name = "Epsilon", signature = c("GHRL")),
    Ductal = gating_model(name = "Ductal", signature = c("CFTR")),
    MUC5B_Ductal = gating_model(name = "MUC5B_Ductal", signature = c("MUC5B")),
    Acinar = gating_model(name = "Acinar", signature = c("REG1A", "CTRB2", "PRSS1", "PRSS2")),
    Stellate = gating_model(name = "Stellate", signature = c("PDGFRB", "COL6A1")),
    Quiescent_Stellate = gating_model(name = "Quiescent_Stellate", signature = c("RGS5")),
    Endothelial = gating_model(name = "Endothelial", signature = c("PLVAP", "ESAM", "VWF")),
    Mast = gating_model(name = "Mast", signature = c("KIT", "CD69")),
    Macrophages = gating_model(name = "Macrophages", signature = c("C1QA", "C1QB", "C1QC"))
)
# Print models to check
lapply(models, print)

tar_source()
obj <- load_seurat("/scratch/domeally/DCD.tienhoven_scRNAseq.2024/analysis_cache/ddqc_out/HPAP-020_ddqc.qs")

# Example: Run the model for Acinar cells on your Seurat object
obj999 <- scGate(data = obj, model = original, verbose = TRUE, ncores = hprcc::slurm_allocation()$CPUs, assay = "RNA", pos.thr = 0.999, slot = "counts", return.CellOntology = FALSE)
DimPlot(obj999, group.by = "scGate_multi") + theme(aspect.ratio = 1)



### UCell

library(Seurat)
library(UCell)
set.seed(123)

# Define the list of signatures
signatures <- list(
    Acinar = c("PNLIP", "CPA1", "CPB1", "PRSS2", "KLK1", "CUZD1", "AMY2A", "CTRB2", "CTRL", "CLPS"),
    Alpha = c("GCG", "CRH", "SPOCK3", "F10", "IRX2", "POPDC3", "C5orf38", "GC", "SPINK4", "TMEM236"),
    Beta = c("MAFA", "IAPP", "ADCYAP1", "INS", "AC132217.2", "LRRTM3", "HHATL", "SAMD11", "C1orf127", "WSCD2"),
    Delta = c("SST", "GPC5-AS11", "LRFN5", "LY6H1", "BCHE", "CALB11", "LINC01571", "CBLN4", "CRYGD", "F5"),
    Ductal = c("FGFBP1", "CFTR", "MMP7", "KRT23", "S100A14", "KRT19", "LGALS4", "CEACAM7", "VTCN1", "TRPV6"),
    Endothelial = c("PLVAP", "KDR", "ESM1", "VWF", "PCAT19", "PECAM1", "CLDN5", "CLEC14A", "ADGRL4", "FCN3"),
    Gamma_Epsilon = c("PPY", "CARTPT", "NPW", "GPC5-AS1", "CALB1", "SERTM1", "LY6H", "GHRL", "ERICH3", "PRG4"),
    Macrophage = c("C1QA", "C1QB", "C1QC", "TYROBP", "SDS", "LAPTM5", "MS4A7", "RGS1", "ITGB2", "GPR183"),
    Mast = c("TPSB2", "TPSAB1", "CD69", "SAMSN1", "CD37", "CD48", "KIT", "GCSAML", "SRGN", "ITGAX"),
    Stellate = c("COL3A1", "DCN", "LUM", "SFRP2", "PTGDS", "PRRX1", "COL1A2", "BGN", "C7", "COL6A3"),
    Cycling_alpha = c("UBE2C", "CENPF", "BIRC5", "TOP2A", "FOXM1", "NUSAP1", "TPX2", "CKAP2L", "MKI67", "CKS2"),
    Active_stellate = c("APOD", "SFRP2", "PTGDS", "CFD", "CCL11", "DPT", "C7", "PODN", "CHI3L1", "CYP1B1"),
    Quiescent_stellate = c("FABP4", "RGS5", "PDK4", "ADIRF", "MCAM", "MUSTN1", "ESAM", "EFHD1", "ADAMTS9", "ADAMTS4"),
    MUC5B_ductal = c("TFF2", "TFF1", "MGST1", "TFF3", "TCN1", "FXYD3", "SPINK1", "MMP1", "C19orf33", "LDHB"),
    ku_ductal = c("CFTR+", "SPP1+", "SOX9+", "KRT19+", "AMY2A-", "CPA1-", "PTF1A-", "GP2-", "CEL-"),
    chen_cancer = c("BCL2L1+", "ATP6V0D1+", "ERO1A+", "RNF139+", "BFAR+", "ARFGAP1+", "MAP3K5+", "PLA2G4B-", "TSPYL2-", "NCCRP1-"),
    chronic_er_stress = c("ATF4+", "GADD34+", "WFS1+", "VIMP+", "SEL1L+", "DERL2+", "P4HB+", "NKX2-2-", "PDX1-", "MAFA-", "INS-", "SLC2A2-"),
    active_er_stress = c("BIP+", "CHOP+", "XBP1+", "ATF6+", "IRE1α+", "PERK+")
    islet_er_stress = c("ATF4+", "GADD34+", "WFS1+", "VIMP+", "SEL1L+", "DERL2+", "P4HB+","BIP+", "CHOP+", "XBP1+", "ATF6+", "IRE1α+", "PERK+", "ERO1A+","RNF139+", "BFAR+", "ARFGAP1+", "MAP3K5+","NKX2-2-", "PDX1-", "MAFA-", "INS-", "SLC2A2-", "PLA2G4B-", "TSPYL2-", "NCCRP1-")
    er_stress = c("ERN1+", "ATF6+", "XBP1+", "EIF2AK3+", "HSPA5+", "DDIT3+", "ATF3+", "ATF4+", "TXNIP+", "INS-")
)

# chen_cancer: https://www.frontiersin.org/articles/10.3389/fmolb.2023.1298077/full

obj.score <- AddModuleScore_UCell(obj, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

VlnPlot(obj.score, features = signature.names, group.by = "seurat_clusters")


#---
library(ggplot2)
library(ggridges)
library(dplyr)
tar_load(sig_scores_per_cell)

sig_scores_long <- sig_scores_per_cell %>%
    pivot_longer(cols = -cell, names_to = "score_type", values_to = "score")
