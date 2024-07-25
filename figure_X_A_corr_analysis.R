# Load required libraries
library(dplyr)
library(corrplot)
library(lme4)
library(performance)
library(see)
library(parameters)
library(corrplot)


# Source necessary files
source("/scratch/domeally/DCD.tienhoven_scRNAseq.2024/R/AAA_defaults_and_constants.R")

# Load or create metadata
if (file.exists("metadata.qs")) {
    metadata <- qs::qread("metadata.qs")
} else {
    tar_load(seurat_object_lognorm_annotated)
    metadata <- seurat_object_lognorm_annotated[[]]
    qs::qsave(metadata, file = "metadata.qs")
}

# Data preparation
plot_data <- metadata %>%
    dplyr::filter(cell_type == "Beta", diabetes_status == "NODM") %>%
    mutate(
        mature_counts_INS = spliced_counts_INS,
        nascent_INS = unspliced_counts_INS,
        total_counts = nCount_RNA,
        mature_cpm_INS = (mature_counts_INS / total_counts) * 1e6,
        nascent_cpm_INS = (nascent_INS / total_counts) * 1e6,
        log_mature_cpm_INS = log1p(mature_cpm_INS),
        log_nascent_cpm_INS = log1p(nascent_cpm_INS)
    ) %>%
    mutate(across(
        c(
            total_counts, percent_mt, chronic_er_stress_score,
            active_er_stress_score, islet_stress_score,
            cellular_stress_score, Beta_UCell,
            chronic_er_stress_UCell, active_er_stress_UCell,
            islet_er_stress_UCell, islet_stress_UCell,
            cellular_stress_UCell, core_upr_stress_UCell,
            msigdb_upr_stress_UCell
        ),
        ~ datawizard::standardize(.)
    ))


# Identify stress_UCell columns
stress_ucell_cols <- grep("stress_UCell", names(plot_data), value = TRUE)

# Update vars_for_correlation
vars_for_correlation <- c(
    "log_mature_cpm_INS", "log_nascent_cpm_INS", "total_counts",
    "percent_mt", "islet_stress_score", "cellular_stress_score", "Beta_UCell",
    "percent_rb", "nCount_RNA", "nFeature_RNA", "percent_hb",
    "percent_pl", "percent_xist", "percent_chrY", "sample_age",
    "sample_xbp1u_psi", "sample_ethnicity", "technology", "protected",
    "rs3842752_consensus", "rs3842753_consensus", "rs689_consensus",
    stress_ucell_cols
)

# Create a subset of data with selected variables
cor_data <- plot_data[, vars_for_correlation]

# Convert all variables to numeric
cor_data <- cor_data %>%
    mutate(across(everything(), ~ as.numeric(as.factor(.))))

# Calculate Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Create correlation plot
png("comprehensive_pearson_correlation_plot.png", width = 2000, height = 1800, res = 150)
corrplot(cor_matrix,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.srt = 45,
    addCoef.col = "black",
    number.cex = 0.5
)
dev.off()

# Create heatmap
png("correlation_heatmap.png", width = 2000, height = 1800, res = 150)
heatmap(cor_matrix,
    symm = TRUE,
    margins = c(10, 10)
)
dev.off()
