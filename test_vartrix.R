
tar_source()

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





vcf <- readr::read_tsv(file = vartrix_vcf, comment = "##",show_col_types = FALSE)
vcf |>
    mutate(FORMAT = "GT", Patient_1 = "0/0") |> head()

 