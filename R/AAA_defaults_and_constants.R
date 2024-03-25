analysis_cache <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024"

scrnaseq_release <- "2.5.1"

failed_qc_donor_ids <- paste0(
            c("HPAP-021|HPAP-023|HPAP-027|"), # MultiQC v2
            c("HPAP-038|HPAP-093"), # MultiQC v3 
            collapse="|")

metadata <- tar_read(pancdb_metadata)
nodm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "NODM" ] |> na.omit() |> paste0(collapse="|")
t1dm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "T1DM" ] |> na.omit() |> paste0(collapse="|")
t2dm_donor_ids <- metadata$donor_id[metadata$diabetes_status == "T2DM" ] |> na.omit() |> paste0(collapse="|")
rm(metadata)

custom_palette <- c(
  # Endocrine
  "Alpha" = "#2ECC71",              # Green
  "alpha" = "#2ECC71",              # Green
  "Beta" = "#3498DB",               # Blue
  "beta" = "#3498DB",               # Blue
  "Delta" = "#1ABC9C",              # Teal
  "delta" = "#1ABC9C",              # Teal
  "Gamma" = "#16A085",              # Dark Teal
  "gamma" = "#16A085",              # Dark Teal
  "epsilon" = "#27AE60",            # Emerald
  
  # Exocrine
  "Acinar-s" = "#E74C3C",           # Red
  "Acinar-i" = "#E67E22",           # Orange
  "Acinar-REG+" = "#F39C12",        # Amber
  "acinar" = "#E74C3C",             # Red
  "Ductal" = "#9B59B6",             # Purple
  "ductal" = "#9B59B6",             # Purple
  "MUC5B+ Ductal" = "#8E44AD",      # Dark Purple
  
  # Immune
  "Macrophage" = "#34495E",         # Dark Grey
  "immune" = "#34495E",             # Pink
  
  # Other
  "Endothelial" = "#95A5A6",        # Grey
  "endothelial" = "#95A5A6",        # Grey
  "Activated Stellate" = "#F1C40F", # Yellow
  "activated_stellate" = "#F1C40F", # Yellow
  "Quiescent Stellate" = "#FDFD96", # Light Yellow
  "quiescent_stellate" = "#FDFD96", # Light Yellow
  "Schwann" = "#2C3E50",            # Dark blue
  "schwann" = "#2C3E50",            # Dark blue
  "cycling" = "#FF7F50"             # Coral
)
