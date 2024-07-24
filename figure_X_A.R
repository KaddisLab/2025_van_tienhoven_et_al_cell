# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(targets)
library(pheatmap)
library(tidyr)
library(viridis)
library(ggpubr)
library(MASS)
library(lme4)
library(performance)
library(see)
library(parameters)
library(datawizard)

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

# Load required libraries
library(dplyr)
library(ggplot2)
library(corrplot)
library(caret)

# Assuming plot_data is already prepared as in the previous script

# Select only the stress-related variables

stress_data <- plot_data %>% dplyr::select(contains("stress"))

# Compute correlation matrix
cor_matrix <- cor(stress_data)

# Visualize correlation matrix
png("stress_measures_correlation.png", width = 800, height = 800)
corrplot(cor_matrix,
    method = "color", type = "upper", order = "hclust",
    tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7
)
dev.off()

# Function to find highly correlated pairs
find_high_correlations <- function(cor_matrix, threshold = 0.7) {
    high_cor <- which(abs(cor_matrix) > threshold & cor_matrix != 1, arr.ind = TRUE)
    high_cor_pairs <- data.frame(
        var1 = rownames(cor_matrix)[high_cor[, 1]],
        var2 = colnames(cor_matrix)[high_cor[, 2]],
        correlation = cor_matrix[high_cor]
    )
    high_cor_pairs <- high_cor_pairs[!duplicated(t(apply(high_cor_pairs[, 1:2], 1, sort))), ]
    high_cor_pairs <- high_cor_pairs[order(-abs(high_cor_pairs$correlation)), ]
    return(high_cor_pairs)
}

# Find highly correlated pairs
high_cor_pairs <- find_high_correlations(cor_matrix, threshold = 0.7)
print(high_cor_pairs)

# Focus on cellular_stress correlations
cellular_stress_cor <- cor_matrix["cellular_stress_score", ]
cellular_stress_cor <- sort(cellular_stress_cor, decreasing = TRUE)
print("Correlations with cellular_stress_score:")
print(cellular_stress_cor)

# Variance Inflation Factor (VIF) analysis
library(car)
vif_model <- lm(cellular_stress_score ~ ., data = stress_data)
vif_results <- vif(vif_model)
print("Variance Inflation Factors:")
print(vif_results[order(-vif_results)])

# Save results to a file
sink("stress_measure_collinearity_results.txt")
cat("Highly correlated pairs (|r| > 0.7):\n")
print(high_cor_pairs)
cat("\nCorrelations with cellular_stress_score:\n")
print(cellular_stress_cor)
cat("\nVariance Inflation Factors:\n")
print(vif_results[order(-vif_results)])
sink()

# Suggestions for variable removal
cat("\nSuggestions for variable removal:\n")
cat("Based on high correlations and VIF values, consider removing these variables:\n")
suggest_remove <- high_cor_pairs$var2[high_cor_pairs$var1 == "cellular_stress_score" |
    high_cor_pairs$var2 == "cellular_stress_score"]
suggest_remove <- c(suggest_remove, names(vif_results)[vif_results > 5 &
    names(vif_results) != "cellular_stress_score"])
suggest_remove <- unique(suggest_remove)
print(suggest_remove)

# Updated model without highly collinear variables
vars_to_keep <- setdiff(stress_vars, suggest_remove)
formula_string <- paste("log_mature_cpm_INS ~", paste(c(
    "protected", "rs3842752_consensus", "rs3842753_consensus",
    "rs689_consensus", "total_counts", "percent_mt", vars_to_keep, "Beta_UCell", "technology",
    "tissue_source"
), collapse = " + "), "+ (1 | orig.ident)")

updated_model <- lmer(as.formula(formula_string), data = plot_data)
print(summary(updated_model))