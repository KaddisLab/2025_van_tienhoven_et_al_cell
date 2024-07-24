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
            cellular_stress_score, Beta_UCell
        ),
        ~ datawizard::standardize(.)
    ))

# Refit the models
lmm_mature <- lmer(
    log_mature_cpm_INS ~ protected + total_counts +
        percent_mt + chronic_er_stress_score + active_er_stress_score +
        islet_stress_score + cellular_stress_score +
        Beta_UCell + (1 | orig.ident),
    data = plot_data
)

lmm_nascent <- lmer(
    log_nascent_cpm_INS ~ protected + total_counts +
        percent_mt + chronic_er_stress_score + active_er_stress_score +
        islet_stress_score + cellular_stress_score +
        Beta_UCell + (1 | orig.ident),
    data = plot_data
)

# Model diagnostics
check_model(lmm_mature)
check_model(lmm_nascent)

# Model performance
print(model_performance(lmm_mature))
print(model_performance(lmm_nascent))

# Visualize model parameters
plot_parameters <- function(model, title) {
    plot(parameters(model)) +
        ggtitle(title) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
}

mature_params_plot <- plot_parameters(lmm_mature, "Mature INS Expression Model Parameters")
nascent_params_plot <- plot_parameters(lmm_nascent, "Nascent INS Expression Model Parameters")

# Combine parameter plots
combined_params_plot <- mature_params_plot + nascent_params_plot +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "Model Parameters for INS Expression",
        theme = theme(plot.title = element_text(hjust = 0.5))
    )

print(combined_params_plot)

# Save the combined parameter plot
ggsave("figure_model_parameters.png", combined_params_plot, width = 16, height = 10, units = "in", dpi = 300)

# Model comparison
lmm_mature_no_protected <- update(lmm_mature, . ~ . - protected)
lmm_nascent_no_protected <- update(lmm_nascent, . ~ . - protected)

print(compare_performance(lmm_mature, lmm_mature_no_protected))
print(compare_performance(lmm_nascent, lmm_nascent_no_protected))

# Print model summaries
print(summary(lmm_mature))
print(summary(lmm_nascent))
