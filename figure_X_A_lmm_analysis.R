# Load required libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(performance)
library(see)
library(parameters)
library(datawizard)
library(lmerTest)
library(emmeans)

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

# Define initial variables for analysis
vars_for_analysis <- c(
    "protected", "rs3842752_consensus", "rs3842753_consensus", "rs689_consensus",
    "total_counts", "percent_mt", "islet_stress_score", "cellular_stress_score",
    "islet_stress_UCell", "cellular_stress_UCell",
    "core_upr_stress_UCell", "msigdb_upr_stress_UCell", "Beta_UCell",
    "technology", "tissue_source"
)

# Ensure all categorical variables are factors
plot_data <- plot_data %>%
    mutate(across(all_of(c(vars_for_analysis, "orig.ident")), as.factor))

# Fit models using lmerTest
mature_model <- lmerTest::lmer(
    log_mature_cpm_INS ~ protected + rs3842752_consensus + rs3842753_consensus +
        rs689_consensus + total_counts + percent_mt + islet_stress_score +
        cellular_stress_score + islet_stress_UCell + cellular_stress_UCell +
        core_upr_stress_UCell + msigdb_upr_stress_UCell + Beta_UCell +
        technology + tissue_source + (1 | orig.ident),
    data = plot_data
)

nascent_model <- lmerTest::lmer(
    log_nascent_cpm_INS ~ protected + rs3842752_consensus + rs3842753_consensus +
        rs689_consensus + total_counts + percent_mt + islet_stress_score +
        cellular_stress_score + islet_stress_UCell + cellular_stress_UCell +
        core_upr_stress_UCell + msigdb_upr_stress_UCell + Beta_UCell +
        technology + tissue_source + (1 | orig.ident),
    data = plot_data
)

# Model summaries
print(summary(mature_model))
print(summary(nascent_model))

# Model diagnostics
check_model(mature_model)
check_model(nascent_model)

# Model performance
print(model_performance(mature_model))
print(model_performance(nascent_model))

# Estimated Marginal Means
# Let's look at the effect of 'protected' while averaging over other predictors
emm_mature <- emmeans(mature_model, ~protected)
emm_nascent <- emmeans(nascent_model, ~protected)

print(emm_mature)
print(emm_nascent)

# Pairwise comparisons
pairs(emm_mature)
pairs(emm_nascent)

# Visualize emmeans results
plot(emm_mature) + ggtitle("Estimated Marginal Means for Mature INS Expression")
plot(emm_nascent) + ggtitle("Estimated Marginal Means for Nascent INS Expression")

# If you want to look at interactions, you can do something like:
emm_interaction_mature <- emmeans(mature_model, ~ protected * technology)
plot(emm_interaction_mature) + ggtitle("Interaction of Protected and Technology for Mature INS Expression")

# You can also use emmeans for more complex contrasts or custom comparisons
# For example, to compare the effect of 'protected' across different technologies:
emm_nested <- emmeans(mature_model, ~ protected | technology)
pairs(emm_nested)