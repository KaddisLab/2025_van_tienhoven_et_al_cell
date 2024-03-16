

 p1<- Seurat::DimPlot(so,  group.by = "predicted.annotation.l1", label = TRUE, label.size = 3, reduction = "ref.umap") + NoLegend()


tar_source()
analysis_cache <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024"

so<-load_seurat("/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/cellbender_seurat_objects/HPAP-023.qs")
seurat_object<-so

so <- NormalizeData(so)
so <- ScaleData(so)

Idents(so) <- "predicted.celltype.l1"

p1 <- FeaturePlot(so, features = "INS", keep.scale = NULL, reduction = "ref.umap")
 ggsave("p1.png", plot = p1, width = 10, height = 10, units = "in", dpi = 300)

# Define your custom color palette
custom_palette <- c(
  "Acinar-s" = "#E74C3C",           # Red
  "Acinar-i" = "#E67E22",           # Orange
  "Acinar-REG+" = "#F39C12",        # Amber
  "Endothelial" = "#95A5A6",        # Grey
  "Ductal" = "#9B59B6",             # Purple
  "Activated Stellate" = "#F1C40F", # Yellow
  "Quiescent Stellate" = "#F7DC6F", # Light Yellow
  "Beta" = "#3498DB",               # Blue
  "Schwann" = "#ABB2B9",            # Light Grey
  "Delta" = "#1ABC9C",              # Teal
  "MUC5B+ Ductal" = "#8E44AD",      # Dark Purple
  "Macrophage" = "#34495E",         # Dark Grey
  "Alpha" = "#2ECC71",              # Green
  "Gamma" = "#16A085"               # Dark Teal
)


##ref
seurat_object<- targets::tar_read(tosti_etal_seurat_object) |> load_seurat()
pancreas.ref <- NormalizeData(seurat_object)
pancreas.ref <- FindVariableFeatures(pancreas.ref)
pancreas.ref <- ScaleData(pancreas.ref)
pancreas.ref <- RunPCA(pancreas.ref)
pancreas.ref <- FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref, res = 1.5)
pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, return.model = TRUE)
pancreas.ref$Cluster <- factor(pancreas.ref$Cluster, levels = names(custom_palette))
p1 <- DimPlot(pancreas.ref, group.by = "Cluster", cols = custom_palette, label = TRUE, label.size = 6, reduction = "umap") & NoLegend()
ggsave("p1.png", plot = p1, width = 10, height = 10, units = "in", dpi = 300)


# After preprocessing, we integrate layers with added parameters specific to Harmony:
obj <- IntegrateLayers(object = ref_seurat_object, method = HarmonyIntegration, orig.reduction = "pca",
  new.reduction = 'harmony', verbose = FALSE)



###-------


snp_vcf_matrix<-readr::read_delim("/scratch/domeally/DCD.tienhoven_scRNAseq.2024/data/snp_vcf_matrix.tsv")

rs3842752 <- snp_vcf_matrix %>% dplyr::filter(VEP_ID == "rs3842752") 

burgertools::ReadVartrix(
    genotype = rs3842752, 
    ref = "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/vartrix_out/HPAP-020_coverage/ref_matrix.mtx",
    alt = "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/vartrix_out/HPAP-020_coverage/alt_matrix.mtx",
    barcodes = "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/vartrix_out/HPAP-020_coverage/colnames.tsv")


save_filtered_vcf <- function(input_file = "", vep_ids_to_keep = c(), output_file = "filtered_vcf.vcf") {
  library(dplyr)
  library(readr)
  library(stringr)
  
  # Read the VCF file, skipping initial VCF header lines
  vcf <- read_tsv(input_file, comment = "##", show_col_types = FALSE)
  colnames(vcf)[1] <- str_replace(colnames(vcf)[1], pattern = "#", replacement = "")
  
  # Generate unique row names of the form chr-position-ref-alt, assuming the input format is consistent with your snippet
  cat("Generating unique row names of the form chr-position-ref-alt.\n")
  UNIQUE_ID <- paste0(vcf$CHROM, "-", vcf$POS, "-", vcf$REF, "-", vcf$ALT)
  vcf$UNIQUE_ID <- UNIQUE_ID
  
  # Optionally filter by VEP IDs
  if (!is.null(vep_ids_to_keep) && length(vep_ids_to_keep) > 0) {
    vcf <- vcf %>% filter(VEP_ID %in% vep_ids_to_keep)
  }
  
  # Add the FORMAT column with a placeholder value ("GT") and a GENOTYPE column ("0/0" as a placeholder)
  vcf <- vcf %>% mutate(FORMAT = "GT", GENOTYPE = "0/0")
  
  # Prepare the header lines for the VCF file
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste("#", paste(colnames(vcf), collapse = "\t"), sep = "")
  )
  
  cat("Writing filtered and formatted VCF file.\n")
  write_lines(vcf_header, output_file) # Write the header
  write_tsv(vcf, output_file, append = TRUE, col_names = FALSE) # Append the data without column names
  
  return(output_file)
}

# Usage example:
# save_filtered_vcf(input_file = "/path/to/your/input_file.vcf", 
#                   vep_ids_to_keep = c("rs123", "rs456"), 
#                   output_file = "path/to/your/output_file.vcf")

save_filtered_vcf(snp_vcf, vep_ids_to_keep = c("rs3842752"))


vcf_2_burger <- function(file = "", output_path = "filtered_vcf.vcf", vep_ids_to_keep = NULL) {
  library(readr)
  library(stringr)
  library(glue)
  library(dplyr)
  
  cat("Reading VCF file.\n")
  vcf <- read_tsv(file = file, comment = "##", show_col_types = FALSE)
  colnames(vcf)[1] <- str_replace(colnames(vcf)[1], pattern = "#", replacement = "")
  
  # Filter by VEP_ID if specified
  if (!is.null(vep_ids_to_keep) && length(vep_ids_to_keep) > 0) {
    cat("Filtering based on VEP IDs.\n")
    vcf <- vcf %>% dplyr::filter(ID %in% vep_ids_to_keep)
  }
  
  # Add the FORMAT column with "GT" as a placeholder value
  vcf <- vcf %>% mutate(FORMAT = "GT", GENOTYPE = "0/0")
  
  cat("Writing processed VCF file.\n")
  # Write the VCF header to the output file
  vcf_header <- c("##fileformat=VCFv4.2",
                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                  paste0("#", paste(names(vcf), collapse="\t")))
  write_lines(vcf_header, output_path)
  
  # Write the data to the file, appending to the header
  write_tsv(vcf, output_path, append = TRUE, col_names = FALSE)
  
  return(output_path)
}

vcf_2_burger(snp_vcf, vep_ids_to_keep = c("rs3842752"))

# Example usage:
# output_file_path <- vcf_2_burger(file = "/path/to/your/input_file.vcf", 
#                                   output_path = "/path/to/your/output_file.vcf",
#                                   vep_ids_to_keep = c("rs123", "rs456"))
# cat("Output file saved to:", output_file_path, "\n")




---

library(Matrix)

load_vartrix_mtx <- function(mtx_path) {

  # Extract the directory from the mtx_path
  run_path <- dirname(mtx_path)
  
  # Reading the matrix as a dgTMatrix
  ref_matrix <- as(readMM(mtx_path), "dgTMatrix")
  
  # Construct paths for row and column names files based on mtx_path
  rownames_path <- file.path(run_path, "rownames.tsv")
  colnames_path <- file.path(run_path, "colnames.tsv")
  
  # Reading row names and column names
  # Assuming the files have no header and the names are in the first column
  rownames_vec <- read.delim(rownames_path, header = FALSE)[,1]
  colnames_vec <- read.delim(colnames_path, header = FALSE)[,1]
  
  # Assigning row names and column names to the matrix
  dimnames(MM_matrix) <- list(rownames_vec, colnames_vec)
  
  return(MM_matrix)
}
ref_path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/vartrix_out/HPAP-037_coverage/ref_matrix.mtx"

ref_matrix <- read_and_annotate_mtx(ref_path)
alt_matrix <- read_and_annotate_mtx(sub("ref_matrix", "alt_matrix", ref_path))

ref_matching_rows <- grep("chr11_2159842", rownames(ref_matrix), value = TRUE)
alt_matching_rows <- grep("chr11_2159842", rownames(alt_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_ref <- ref_matrix[ref_matching_rows, ]
rs3842752_alt <- alt_matrix[alt_matching_rows, ]
rs3842752_prop <- round(rs3842752_alt / (rs3842752_ref + rs3842752_alt), 2)

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))


ref_path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/vartrix_out/HPAP-042_coverage/ref_matrix.mtx"
ref_matrix <- read_and_annotate_mtx(ref_path)
alt_matrix <- read_and_annotate_mtx(sub("ref_matrix", "alt_matrix", ref_path))

ref_matching_rows <- grep("chr11_2159842", rownames(ref_matrix), value = TRUE)
alt_matching_rows <- grep("chr11_2159842", rownames(alt_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_ref <- ref_matrix[ref_matching_rows, ]
rs3842752_alt <- alt_matrix[alt_matching_rows, ]
rs3842752_prop <- round(rs3842752_alt / (rs3842752_ref + rs3842752_alt), 2)

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))



ref_path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/vartrix_out/HPAP-064_coverage/ref_matrix.mtx"
ref_matrix <- read_and_annotate_mtx(ref_path)
alt_matrix <- read_and_annotate_mtx(sub("ref_matrix", "alt_matrix", ref_path))

ref_matching_rows <- grep("chr11_2159842", rownames(ref_matrix), value = TRUE)
alt_matching_rows <- grep("chr11_2159842", rownames(alt_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_ref <- ref_matrix[ref_matching_rows, ]
rs3842752_alt <- alt_matrix[alt_matching_rows, ]
rs3842752_prop <- round(rs3842752_alt / (rs3842752_ref + rs3842752_alt), 2)

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))



ref_path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/vartrix_out/HPAP-099_coverage/ref_matrix.mtx"
ref_matrix <- read_and_annotate_mtx(ref_path)
alt_matrix <- read_and_annotate_mtx(sub("ref_matrix", "alt_matrix", ref_path))

ref_matching_rows <- grep("chr11_2159842", rownames(ref_matrix), value = TRUE)
alt_matching_rows <- grep("chr11_2159842", rownames(alt_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_ref <- ref_matrix[ref_matching_rows, ]
rs3842752_alt <- alt_matrix[alt_matching_rows, ]
rs3842752_prop <- round(rs3842752_alt / (rs3842752_ref + rs3842752_alt), 2)

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))


con_path<-"/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/vartrix_out/HPAP-051_consensus/alt_matrix.mtx"
con_matrix <- read_and_annotate_mtx(ref_path)

matching_rows <- grep("chr11_2159842", rownames(con_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_con <- con_matrix[matching_rows, ]

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))


ref_path <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/vartrix_out/HPAP-051_coverage/ref_matrix.mtx"
ref_matrix <- read_and_annotate_mtx(ref_path)
alt_matrix <- read_and_annotate_mtx(sub("ref_matrix", "alt_matrix", ref_path))

ref_matching_rows <- grep("chr11_2159842", rownames(ref_matrix), value = TRUE)
alt_matching_rows <- grep("chr11_2159842", rownames(alt_matrix), value = TRUE)

# Extract the subset of the matrix with row names matching the pattern
rs3842752_ref <- ref_matrix[ref_matching_rows, ]
rs3842752_alt <- alt_matrix[alt_matching_rows, ]
rs3842752_prop <- round(rs3842752_alt / (rs3842752_ref + rs3842752_alt), 2)

View(data.frame(rs3842752_ref, rs3842752_alt, rs3842752_prop))




df <- left_join(protected_cohort, pancdb_metadata, by = c("sample_id" = "donor_id"))
df <- df %>%
  mutate(disease_status_clean = case_when(
    grepl("No Hx DIAB|No HX DIAB|No HX Diabetes", disease_status) ~ "No Hx DIAB",
    grepl("T1DM", disease_status) ~ "T1DM",
    grepl("T2DM", disease_status) ~ "T2DM",
    TRUE ~ "Other"
  ))


# Plot for 'sample_sex'
sample_sex_plot <- ggplot(df, aes(x = sample_sex, fill = sample_sex)) +
  geom_bar() +
  labs(x = "Sample Sex", y = "Count", title = "Distribution of Sample Sex") +
  scale_fill_brewer(palette = "Set3")

# Plot for 'age'
age_plot <- ggplot(df, aes(x = sample_age)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(x = "Sample Age", y = "Count", title = "Distribution of Sample Age")

# Updated plot for 'disease_status_clean'
disease_status_plot <- ggplot(df, aes(x = disease_status_clean, fill = disease_status_clean)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Disease Status", y = "Count", title = "Consolidated Distribution of Disease Status") +
  scale_fill_brewer(palette = "Set2")

# Plot for 'sample_ethnicity'
ethnicity_plot <- ggplot(df, aes(x = sample_ethnicity, fill = sample_ethnicity)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample Ethnicity", y = "Count", title = "Distribution of Sample Ethnicity") +
  scale_fill_brewer(palette = "Pastel1")

# Update combined plot to include 'ethnicity_plot'
combined_plot <- (sample_sex_plot | age_plot ) / 
                 ( disease_status_plot | ethnicity_plot)

# Save updated combined plot to PNG
ggsave("combined_plots.png", combined_plot, width = 24, height = 12, dpi = 300)





vcf <- readr::read_tsv(file = vartrix_vcf, comment = "##",show_col_types = FALSE)
vcf |>
    mutate(FORMAT = "GT", Patient_1 = "0/0") |> head()

