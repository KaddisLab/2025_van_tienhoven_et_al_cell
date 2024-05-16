library(dplyr)
library(Matrix)

AnnotateGenes.10X <- function(MarvelObject) {
    gtf <- MarvelObject$gtf
    gene.metadata <- MarvelObject$gene.metadata

    # Read GTF data - simulate GTF reading based on your data structure
    gtf_df <- readr::read_tsv(gtf) %>%
        filter(V3 == "gene") %>%
        transmute(
            gene_short_name = str_extract(V9, "(?<=gene_name \")[^\"]+"),
            gene_type = str_extract(V9, "(?<=gene_biotype \")[^\"]+") %>% coalesce(str_extract(V9, "(?<=gene_type \")[^\"]+"))
        ) %>%
        distinct()

    # Ensure that gene_short_name in both data frames are character
    gene.metadata$gene_short_name <- as.character(gene.metadata$gene_short_name)
    gtf_df$gene_short_name <- as.character(gtf_df$gene_short_name)

    # Join with GTF annotations
    gene.metadata <- gene.metadata %>%
        left_join(gtf_df, by = "gene_short_name")

    # Error handling for non-matching entries
    if (any(is.na(gene.metadata$gene_type))) {
        warning("Some genes did not match the GTF annotations.")
    }

    # Updating the Marvel object
    MarvelObject$gene.metadata <- gene.metadata
    MarvelObject$gene.norm.matrix <- Matrix::Matrix(MarvelObject$gene.norm.matrix[gene.metadata$gene_short_name, , drop = FALSE], sparse = TRUE)

    return(MarvelObject)
}
