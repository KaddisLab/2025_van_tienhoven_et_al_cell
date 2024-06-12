

seurat_INS_hknorm <- function(seurat_object) {
    
    seurat_object <- "/scratch/domeally/DCD.tienhoven_scRNAseq.2024/analysis_cache/ddqc_out/HPAP-019_ddqc.qs"
    tar_source()
    seurat_object <- load_seurat(seurat_object)
    
    # Get the project name
    project_name <- Seurat::Project(seurat_object)

    missing_genes <- setdiff(housekeeping_genes, rownames(seurat_object[["RNA"]]$counts))
    if (length(missing_genes) > 0) {
        stop("The following housekeeping genes are not present in the Seurat object: ", paste(missing_genes, collapse = ", "))
    }

    # Extract counts for INS and housekeeping genes
    INS_counts <- seurat_object[["RNA"]]$counts["INS", ]
    housekeeping_counts <- seurat_object[["RNA"]]$counts[housekeeping_genes, ]

    # Sum of housekeeping gene counts
    sum_housekeeping_counts <- colSums(housekeeping_counts)

    # Normalize INS counts by the sum of housekeeping gene counts,
    # scaled to a total of 10,000 and natural log transformed
    INS_normalized <- log1p((INS_counts / sum_housekeeping_counts) * 10000)

    # Replace infinite values with 0, ie some housekeeping genes have 0 
    # counts so we cant normalise INS expression in those cells
    INS_normalized[is.infinite(INS_normalized)] <- 0

    INS_hk_tibble <- tibble(cell = names(INS_normalized), INS_hk = as.numeric(INS_normalized))

    INS_hk_tibble_path <- glue::glue("{analysis_cache}/INS_hknorm_out/{project_name}_INS_hknorm_per_cell.csv")
    dir.create(dirname(ratio_tibble_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(ratio_tibble, ratio_tibble_path, row.names = FALSE, quote = FALSE)
    message("Saved ratio tibble")

    return(INS_hk_tibble)
}
