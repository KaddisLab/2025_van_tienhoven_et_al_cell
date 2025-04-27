analysis_cache <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024_cache"

scrnaseq_release <- "2.5.1"

failed_qc_donor_ids <- paste0(
    c("HPAP-021|HPAP-023|HPAP-027|"), # v2 multiqc report
    c("HPAP-038|HPAP-093|"), # v3 multiqc report
    c("HPAP-023|HPAP-027|HPAP-070|HPAP-106|HPAP-109|"), # failed_qc: see QC_PCA report
    c("HPAP-118"), # failed nascent INS QC; also poor "NODM" as this donor has high A1C & morbid obesity
    collapse = "|"
)

# BMI Matched (NODM subset)
bmim_donor_ids <- c(
    "HPAP-047|HPAP-039|HPAP-104|HPAP-026|HPAP-099|HPAP-053|HPAP-097|HPAP-119|HPAP-117|HPAP-077|HPAP-074|HPAP-103|HPAP-037|HPAP-082|HPAP-042|HPAP-044"
)


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

cell_cycle_palette <- c(
    G1 = "#1f77b4",
    S = "#ff7f0e",
    G2M = "#2ca02c"
)

diabetes_palette <- c(
    "NODM" = "#3fa36b", # A cool green for 'No Diabetes'
    "AABP" = "#F4C542", # A warm yellow for 'Autoantibody Positive'
    "T1DM" = "#9353b3", # A warm purple for 'Type 1 Diabetes'
    "T2DM" = "#D8604C" # A warm red for 'Type 2 Diabetes'
)

housekeeping_genes <- c(
    "ACTB",
    "GAPDH",
    "PGK1",
    "PPIA",
    "RPLP0",
    "SDHA",
    "TFRC",
    "GUSB",
    "HMBS",
    "HPRT1",
    "TBP"
)

control_genes <- c(
    "SLC11A2",
    "SLC40A1",
    "UBC",
    "RPLP0",
    "ERN1",
    "TBP",
    "ACTB",
    "GAPDH",
    "PGK1",
    "TUBB",
    "PPIA",
    "SDHA",
    "TFRC",
    "GUSB",
    "HMBS",
    "HPRT1"
)
signatures <- list(
    cellular_stress = c("DDIT3", "PSMB10", "ATF3", "TXNIP"),
    msigdb_upr_stress = c(
        "ALDH18A1",
        "ARFGAP1",
        "ASNS",
        "ATF3",
        "ATF4",
        "ATF6",
        "ATP6V0D1",
        "BAG3",
        "BANF1",
        "CALR",
        "CCL2",
        "CEBPB",
        "CEBPG",
        "CHAC1",
        "CKS1B",
        "CNOT2",
        "CNOT4",
        "CNOT6",
        "CXXC1",
        "DCP1A",
        "DCP2",
        "DCTN1",
        "DDIT4",
        "DDX10",
        "DKC1",
        "DNAJA4",
        "DNAJB9",
        "DNAJC3",
        "EDC4",
        "EDEM1",
        "EEF2",
        "EIF2AK3",
        "EIF2S1",
        "EIF4A1",
        "EIF4A2",
        "EIF4A3",
        "EIF4E",
        "EIF4EBP1",
        "EIF4G1",
        "ERN1",
        "ERO1A",
        "EXOC2",
        "EXOSC1",
        "EXOSC10",
        "EXOSC2",
        "EXOSC4",
        "EXOSC5",
        "EXOSC9",
        "FKBP14",
        "FUS",
        "GEMIN4",
        "GOSR2",
        "H2AFX",
        "HERPUD1",
        "HSP90B1",
        "HSPA5",
        "HSPA9",
        "HYOU1",
        "IARS1",
        "IFIT1",
        "IGFBP1",
        "IMP3",
        "KDELR3",
        "KHSRP",
        "KIF5B",
        "LSM1",
        "LSM4",
        "MTHFD2",
        "NFYA",
        "NFYB",
        "NHP2",
        "NOLC1",
        "NOP14",
        "NOP56",
        "NPM1",
        "OBFC2A",
        "PAIP1",
        "PARN",
        "PDIA5",
        "PDIA6",
        "POP4",
        "PREB",
        "PSAT1",
        "RPS14",
        "RRP9",
        "SDAD1",
        "SEC11A",
        "SEC31A",
        "SERP1",
        "SHC1",
        "SKIV2L2",
        "SLC1A4",
        "SLC30A5",
        "SLC7A5",
        "SPCS1",
        "SPCS3",
        "SRPR",
        "SRPRB",
        "SSR1",
        "STC2",
        "TARS1",
        "TATDN2",
        "TSPYL2",
        "TTC37",
        "TUBB2A",
        "VEGFA",
        "WFS1",
        "WIPI1",
        "XBP1",
        "XPOT",
        "YIF1A",
        "YWHAZ",
        "ZBTB17"
    )
)
