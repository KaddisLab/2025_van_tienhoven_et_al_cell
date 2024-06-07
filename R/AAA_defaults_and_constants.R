analysis_cache <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024_cache"

scrnaseq_release <- "2.5.1"

failed_qc_donor_ids <- paste0(
    c("HPAP-021|HPAP-023|HPAP-027|"), # MultiQC v2
    c("HPAP-038|HPAP-093"), # MultiQC v3
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

cell_type_palette <- c(
    # Endocrine
    "Alpha" = "#2ECC71", # Green
    "alpha" = "#2ECC71", # Green
    "Beta" = "#3498DB", # Blue
    "Beta_like" = "#15bfd9", # Blue
    "beta" = "#3498DB", # Blue
    "Beta_like" = "#15bfd9", # Blue
    "Delta" = "#1ABC9C", # Teal
    "delta" = "#1ABC9C", # Teal
    "Gamma" = "#16A085", # Dark Teal
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
    "Macrophage" = "brown", # Dark Grey
    "macrophage" = "brown", # Dark Grey
    "immune" = "brown", # Pink
    "Immune" = "brown", # Pink

    # Other
    "Other" = "#314c4e", #
    "Endothelial" = "#314c4e", #
    "endothelial" = "#314c4e", #
    "Activated Stellate" = "#F1C40F", # Yellow
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


er_genes_of_interest <- c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "DDIT3", "PSMB10", "INS", "ATF3", "ATF4", "TXNIP", "GAPDH")
# XBP1spliced (XBP1s)

upr_genes <- c("ERN1", "ATF6", "XBP1", "EIF2AK3", "HSPA5", "ATF4")

# ER Stress genes
er_stress_genes <- c("DDIT3", "PSMB10", "ATF3", "TXNIP")

# Housekeeping genes
housekeeping_genes <- c("ACTB", "GAPDH", "PGK1", "PPIA", "RPLP0", "SDHA", "TFRC", "GUSB", "HMBS", "HPRT1", "TBP")
