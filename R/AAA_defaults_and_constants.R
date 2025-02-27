analysis_cache <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024_cache"

scrnaseq_release <- "2.5.1"

failed_qc_donor_ids <- paste0(
    c("HPAP-021|HPAP-023|HPAP-027|"), # v2 multiqc report
    c("HPAP-038|HPAP-093|"), # v3 multiqc report
    c("HPAP-023|HPAP-027|HPAP-070|HPAP-106|HPAP-109|"), # failed_qc: see QC_PCA report
    c("HPAP-118"), # failed nascent INS QC; also poor "NODM" as this donor has high A1C & morbid obesity
    collapse = "|"
)

metadata <- tar_read(pancdb_metadata)
nodm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "NODM"] |>
    na.omit() |>
    paste0(collapse = "|")
aabp_donor_ids <- metadata$donor_id[metadata$diabetes_status == "AABP"] |>
    na.omit() |>
    paste0(collapse = "|")
t1dm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "T1DM"] |>
    na.omit() |>
    paste0(collapse = "|")
t2dm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "T2DM"] |>
    na.omit() |>
    paste0(collapse = "|")
rm(metadata)

# BMI Matched (NODM subset)
bmim_donor_ids <- c("HPAP-047|HPAP-039|HPAP-104|HPAP-026|HPAP-099|HPAP-053|HPAP-097|HPAP-119|HPAP-117|HPAP-077|HPAP-074|HPAP-103|HPAP-037|HPAP-082|HPAP-042|HPAP-044")

# T2DM controls (NODM subset)
t2dc_donor_ids <- c("HPAP-052|HPAP-053|HPAP-054|HPAP-059|HPAP-063|HPAP-074|HPAP-075|HPAP-077|HPAP-080|HPAP-093|HPAP-097|HPAP-101|HPAP-103|HPAP-105|HPAP-117|HPAP-118|HPAP-119")

cell_type_palette <- c(
    # Endocrine
    "Cycling Alpha" = "#57ed96", # Green
    "Alpha" = "#2ECC71", # Green
    "alpha" = "#2ECC71", # Green
    "Beta" = "#3498DB", # Blue
    "Beta_like" = "#15bfd9", # Blue
    "beta" = "#3498DB", # Blue
    "Beta_like" = "#15bfd9", # Blue
    "Alpha+Beta" = "#15bfd9", # Blue
    "Delta" = "#1ABC9C", # Teal
    "delta" = "#1ABC9C", # Teal
    "Gamma" = "#16A085", # Dark Teal
    "Gamma+Epsilon" = "#16A085", # Dark Teal
    "gamma" = "#16A085", # Dark Teal
    "PP_Gamma" = "#16A085", # Dark Teal
    "epsilon" = "#27AE60", # Emerald
    "Epsilon" = "#27AE60", # Emerald

    # Exocrine
    "Acinar-s" = "#E74C3C", # Red
    "Acinar-i" = "#E67E22", # Orange
    "Acinar-REG+" = "#F39C12", # Amber
    "acinar" = "#E74C3C", # Red
    "Acinar" = "#E74C3C", # Red
    "Ductal" = "#9B59B6", # Purple
    "ductal" = "#9B59B6", # Purple
    "MUC5B+ Ductal" = "#8E44AD", # Dark Purple

    # Immune
    "Mast" = "brown", # Dark Grey
    "Macrophage" = "brown", # Dark Grey
    "macrophage" = "brown", # Dark Grey
    "immune" = "brown", # Pink
    "Immune" = "brown", # Pink

    # Other
    "Other" = "#314c4e", #
    "Endothelial" = "#4736c7", #
    "endothelial" = "#4736c7", #
    "Activated Stellate" = "#F1C40F", # Yellow
    "Active Stellate" = "#F1C40F", # Yellow
    "Stellates_Mesenchymal" = "#F1C40F", # Yellow
    "activated_stellate" = "#F1C40F", # Yellow
    "Quiescent Stellate" = "#FDFD96", # Light Yellow
    "quiescent_stellate" = "#FDFD96", # Light Yellow
    "Schwann" = "#2C3E50", # Dark blue
    "schwann" = "#2C3E50", # Dark blue
    "cycling" <- "#FF7F50",
    "Unknown" <- "cornsilk2"
)

cell_cycle_palette <- c(G1 = "#1f77b4", S = "#ff7f0e", G2M = "#2ca02c")

cell_type_markers <- c(
    "acinar" = "CPA1",
    "stellate" = "RSG10",
    "alpha" = "GCG",
    "beta" = "INS",
    "beta" = "MAFA",
    "endocrine" = "PDX1",
    "cycling" = "MKI67",
    "delta" = "SST",
    "ductal" = "KRT19",
    "endothelial" = "VWF",
    "epsilon" = "GHRL",
    "gamma" = "PPY",
    "immune" = "PTPRC",
    "schwann" = "CDH19",
    "stress" = "XBP1"
)

diabetes_palette <- c(
    "NODM" = "#3fa36b", # A cool green for 'No Diabetes'
    "AABP" = "#F4C542", # A warm yellow for 'Autoantibody Positive'
    "T1DM" = "#9353b3", # A warm purple for 'Type 1 Diabetes'
    "T2DM" = "#D8604C" # A warm red for 'Type 2 Diabetes'
)

generate_palette <- function(n) {
    palettes <- c(
        MoMAColors::moma.colors("Smith", n = 5, type = "discrete"),
        MoMAColors::moma.colors("Warhol", n = 15, type = "discrete"),
        MoMAColors::moma.colors("Klein", n = 11, type = "discrete"),
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(9, "Set1")
    )
    return(colorRampPalette(palettes)(n))
}

custom_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
        # Panel border and grid lines
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.spacing.x = ggplot2::unit(0.5, "lines"),
        # Strip plots
        # strip.background = element_blank(),
        strip.text = ggplot2::element_text(size = 16),
        strip.text.x = ggplot2::element_text(size = 20),
        strip.text.y = ggplot2::element_text(size = 20),
        # Legend
        legend.key.size = ggplot2::unit(12, "mm"),
        legend.text = ggplot2::element_text(size = 20),
        legend.title = ggplot2::element_blank(),
        legend.position = "bottom",
        # Title
        plot.title = ggplot2::element_text(size = 26, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 20),
        # Axis lines
        axis.line = ggplot2::element_line(linewidth = 0.5),
        axis.line.x = ggplot2::element_line(linewidth = 0.5),
        axis.line.y = ggplot2::element_line(linewidth = 0.5),
        # Axis ticks
        axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
        axis.ticks.y = ggplot2::element_line(linewidth = 0.5),
        # Axis text
        axis.text.x = ggplot2::element_text(size = 20),
        axis.text.y = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 20),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 20),
        # Text
        text = ggplot2::element_text(family = "Roboto")
    )

er_genes_of_interest <- c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "DDIT3", "PSMB10", "INS", "ATF3", "ATF4", "TXNIP", "GAPDH")
# XBP1spliced (XBP1s)

upr_genes <- c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "ATF4")

# ER Stress genes
er_stress_genes <- c("DDIT3", "PSMB10", "ATF3", "TXNIP")

# Housekeeping genes
housekeeping_genes <- c("ACTB", "GAPDH", "PGK1", "PPIA", "RPLP0", "SDHA", "TFRC", "GUSB", "HMBS", "HPRT1", "TBP")

# cell type signatures
signatures <- list(
    Acinar = c("PNLIP", "CPA1", "CPB1", "PRSS2", "KLK1", "CUZD1", "AMY2A", "CTRB2", "CTRL", "CLPS"),
    Alpha = c("GCG", "CRH", "SPOCK3", "F10", "IRX2", "POPDC3", "C5orf38", "GC", "SPINK4", "TMEM236"),
    Beta = c("MAFA", "IAPP", "ADCYAP1", "INS", "AC132217.2", "LRRTM3", "HHATL", "SAMD11", "C1orf127", "WSCD2"),
    Delta = c("SST", "GPC5-AS11", "LRFN5", "LY6H1", "BCHE", "CALB11", "LINC01571", "CBLN4", "CRYGD", "F5"),
    Ductal = c("FGFBP1", "CFTR", "MMP7", "KRT23", "S100A14", "KRT19", "LGALS4", "CEACAM7", "VTCN1", "TRPV6"),
    Endothelial = c("PLVAP", "KDR", "ESM1", "VWF", "PCAT19", "PECAM1", "CLDN5", "CLEC14A", "ADGRL4", "FCN3"),
    Gamma_Epsilon = c("PPY", "CARTPT", "NPW", "GPC5-AS1", "CALB1", "SERTM1", "LY6H", "GHRL", "ERICH3", "PRG4"),
    epsilon = c("ACSL1", "SPINK1", "SERPINA1", "DEFB1", "RBP4", "TTR", "TM4SF4", "TM4SF5", "HES4"),
    Macrophage = c("C1QA", "C1QB", "C1QC", "TYROBP", "SDS", "LAPTM5", "MS4A7", "RGS1", "ITGB2", "GPR183"),
    Mast = c("TPSB2", "TPSAB1", "CD69", "SAMSN1", "CD37", "CD48", "KIT", "GCSAML", "SRGN", "ITGAX"),
    Stellate = c("COL3A1", "DCN", "LUM", "SFRP2", "PTGDS", "PRRX1", "COL1A2", "BGN", "C7", "COL6A3"),
    Cycling_alpha = c("UBE2C", "CENPF", "BIRC5", "TOP2A", "FOXM1", "NUSAP1", "TPX2", "CKAP2L", "MKI67", "CKS2"),
    Active_stellate = c("APOD", "SFRP2", "PTGDS", "CFD", "CCL11", "DPT", "C7", "PODN", "CHI3L1", "CYP1B1"),
    Quiescent_stellate = c("FABP4", "RGS5", "PDK4", "ADIRF", "MCAM", "MUSTN1", "ESAM", "EFHD1", "ADAMTS9", "ADAMTS4"),
    MUC5B_ductal = c("TFF2", "TFF1", "MGST1", "TFF3", "TCN1", "FXYD3", "SPINK1", "MMP1", "C19orf33", "LDHB"),
    ku_ductal = c("CFTR+", "SPP1+", "SOX9+", "KRT19+", "AMY2A-", "CPA1-", "PTF1A-", "GP2-", "CEL-"),
    chen_cancer = c("BCL2L1+", "ATP6V0D1+", "ERO1A+", "RNF139+", "BFAR+", "ARFGAP1+", "MAP3K5+", "PLA2G4B-", "TSPYL2-", "FIGLA-"),
    chronic_er_stress = c("ATF4+", "PPP1R15A+", "WFS1+", "SELS+", "SEL1L+", "DERL2+", "P4HB+", "NKX2-2-", "PDX1-", "MAFA-", "INS-", "SLC2A2-"),
    active_er_stress = c("HSPA5+", "DDIT3+", "XBP1+", "ATF6+", "ERN1+", "EIF2AK3+"),
    islet_er_stress = c("ATF4+", "PPP1R15A+", "WFS1+", "SELS+", "SEL1L+", "DERL2+", "P4HB+", "HSPA5+", "DDIT3+", "XBP1+", "ATF6+", "ERN1+", "EIF2AK3+", "ERO1A+", "RNF139+", "BFAR+", "ARFGAP1+", "MAP3K5+", "NKX2-2-", "PDX1-", "MAFA-", "INS-", "SLC2A2-", "PLA2G4B-", "TSPYL2-", "NCCRP1-"),
    islet_stress = c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "DDIT3", "ATF3", "ATF4", "TXNIP"),
    cellular_stress = c("DDIT3", "PSMB10", "ATF3", "TXNIP"),
    core_upr_stress = c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "ATF4"),
    msigdb_upr_stress = c("ALDH18A1", "ARFGAP1", "ASNS", "ATF3", "ATF4", "ATF6", "ATP6V0D1", "BAG3", "BANF1", "CALR", "CCL2", "CEBPB", "CEBPG", "CHAC1", "CKS1B", "CNOT2", "CNOT4", "CNOT6", "CXXC1", "DCP1A", "DCP2", "DCTN1", "DDIT4", "DDX10", "DKC1", "DNAJA4", "DNAJB9", "DNAJC3", "EDC4", "EDEM1", "EEF2", "EIF2AK3", "EIF2S1", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4E", "EIF4EBP1", "EIF4G1", "ERN1", "ERO1A", "EXOC2", "EXOSC1", "EXOSC10", "EXOSC2", "EXOSC4", "EXOSC5", "EXOSC9", "FKBP14", "FUS", "GEMIN4", "GOSR2", "H2AFX", "HERPUD1", "HSP90B1", "HSPA5", "HSPA9", "HYOU1", "IARS1", "IFIT1", "IGFBP1", "IMP3", "KDELR3", "KHSRP", "KIF5B", "LSM1", "LSM4", "MTHFD2", "NFYA", "NFYB", "NHP2", "NOLC1", "NOP14", "NOP56", "NPM1", "OBFC2A", "PAIP1", "PARN", "PDIA5", "PDIA6", "POP4", "PREB", "PSAT1", "RPS14", "RRP9", "SDAD1", "SEC11A", "SEC31A", "SERP1", "SHC1", "SKIV2L2", "SLC1A4", "SLC30A5", "SLC7A5", "SPCS1", "SPCS3", "SRPR", "SRPRB", "SSR1", "STC2", "TARS1", "TATDN2", "TSPYL2", "TTC37", "TUBB2A", "VEGFA", "WFS1", "WIPI1", "XBP1", "XPOT", "YIF1A", "YWHAZ", "ZBTB17")
)
# epsilon: Dominguez Gutierrez et al https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6963699
# islet_er_stress: chatGPT
# islet_stress, cellular_stress, core_upr_stress: this study
# chronic/active stress: Chen et al https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9360004/
# chen_cancer: https://www.frontiersin.org/articles/10.3389/fmolb.2023.1298077/full
# ku_ductal: Terasa Ku, COH, personal communication
# msigdb_upr: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_UNFOLDED_PROTEIN_RESPONSE.html
