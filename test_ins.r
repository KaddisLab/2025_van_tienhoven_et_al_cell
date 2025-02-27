# Load necessary libraries
library(Seurat)
library(tidyseurat)
library(ggplot2)
library(ggridges)
library(dplyr)
library(forcats)
library(ggdist)
library(ggExtra)



# Define the function to create the ridge plot
create_ridge_plot <- function(seurat_object, gene, group_var1, group_var2, palette = NULL) {
    # Fetch the data for the specified gene and grouping variables
    data <- FetchData(seurat_object, vars = c(gene, group_var1, group_var2))

    # Define endocrine and exocrine cell types
    endocrine_types <- c("Alpha", "Alpha+Beta", "Beta", "Cycling Alpha", "Delta", "Gamma+Epsilon")
    exocrine_types <- c("Acinar", "Ductal", "MUC5B+ Ductal")

    # Create a new variable to classify cell types into endocrine, exocrine, and other
    data <- data %>%
        mutate(group = case_when(
            cell_type %in% endocrine_types ~ "Endocrine",
            cell_type %in% exocrine_types ~ "Exocrine",
            cell_type == "Other" ~ "Other"
        ))

    # Ensure that the grouping variable is a factor for proper ordering
    data[[group_var2]] <- as.factor(data[[group_var2]])

    # Create a new variable to ensure T/F alternate within each cell type
    data <- data %>%
        mutate(grouped_cell_type = paste(cell_type, .data[[group_var2]], sep = "_"))

    # Create the base plot
    p <- ggplot(data, aes(x = .data[[gene]], y = fct_rev(cell_type), fill = .data[[group_var2]])) +
        geom_density_ridges(aes(group = grouped_cell_type), alpha = 0.5, scale = 0.9, color = NA) +
        facet_grid(rows = vars(group), scales = "free_y", space = "free_y") +
        labs(title = paste("Ridge Plot of", gene, "Grouped by Cell Type and", group_var2), x = gene, y = "Cell Type") +
        theme_classic() +
        theme(legend.position = "top", legend.title = element_blank())

    # Apply custom palette if supplied
    if (!is.null(palette)) {
        p <- p + scale_fill_manual(values = palette)
    }

    # Return the plot
    return(p)
}

# Example usage
# create_ridge_plot(seurat_object, "INS", "cell_type", "protected")
# create_ridge_plot(seurat_object, "INS", "cell_type", "protected", palette = c("TRUE" = "blue", "FALSE" = "red"))


# Define the function to create the dot and interval plot
create_dot_interval_plot <- function(seurat_object, gene, group_var1, group_var2, palette = NULL) {
    # Fetch the data for the specified gene and grouping variables
    data <- FetchData(seurat_object, vars = c(gene, group_var1, group_var2))

    # Define endocrine and exocrine cell types
    endocrine_types <- c("Alpha", "Alpha+Beta", "Beta", "Cycling Alpha", "Delta", "Gamma+Epsilon")
    exocrine_types <- c("Acinar", "Ductal", "MUC5B+ Ductal")

    # Create a new variable to classify cell types into endocrine, exocrine, and other
    data <- data %>%
        mutate(group = case_when(
            cell_type %in% endocrine_types ~ "Endocrine",
            cell_type %in% exocrine_types ~ "Exocrine",
            cell_type == "Other" ~ "Other"
        ))

    # Ensure that the grouping variable is a factor for proper ordering
    data[[group_var2]] <- as.factor(data[[group_var2]])

    # Create a new variable to ensure T/F alternate within each cell type
    data <- data %>%
        mutate(grouped_cell_type = paste(cell_type, .data[[group_var2]], sep = "_"))

    # Create the base plot
    p <- ggplot(data, aes(x = fct_rev(cell_type), y = .data[[gene]], fill = .data[[group_var2]])) +
        stat_dots(position = "dodgejust") +
        stat_pointinterval(color = "black", alpha = 0.6) +
        facet_grid(rows = vars(group), scales = "free_y", space = "free_y") +
        labs(title = paste("Dot and Interval Plot of", gene, "Grouped by Cell Type and", group_var2), x = "Cell Type", y = gene) +
        theme_classic() +
        theme(legend.position = "top", legend.title = element_blank())

    # Apply custom palette if supplied
    if (!is.null(palette)) {
        p <- p + scale_fill_manual(
            values = palette
        )
    }

    # Return the plot
    return(p)
}

# Example usage
# create_dot_interval_plot(seurat_object, "INS", "cell_type", "protected")
# create_dot_interval_plot(seurat_object, "INS", "cell_type", "protected", palette = c("TRUE" = "blue", "FALSE" = "red"))

tar_load(seurat_object_lognorm_annotated)
seurat_object <- seurat_object_lognorm_annotated |>
    dplyr::filter(diabetes_status == "NODM")

data <- FetchData(seurat_object, vars = c("INS_hk", "INS", "protected", "cell_type", "cell_type_extra",  "upr_score"))

ins_data <- data %>%
    mutate(cell_type_plus = paste0(cell_type, ifelse(cell_type_extra != "", paste0("_", cell_type_extra), ""))) |>
    dplyr::mutate(cell_type_plus = gsub("_Epsilon", "", cell_type_plus)) |>
    dplyr::filter(INS > 0)

# geom_point with INK_hk and uprscore by cell type
ins_data |>
    dplyr::filter(cell_type %in% c("Beta", "Alpha+Beta")) |>
    ggplot(aes(x = upr_score, y = INS_hk, color = protected, alpha = 0.2)) +
    geom_point() +
    labs(title = "INS_hk vs upr_score by cell type", x = "upr_score", y = "INS_hk") +
    theme_classic() +
    Seurat::NoLegend() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) 

# geom_point with INK_hk and uprscore by cell type
ins_data |>
    dplyr::filter(cell_type %in% c("Beta", "Alpha+Beta")) |>
    ggplot(aes(x = upr_score, y = INS, color = protected, alpha = 0.2)) +
    geom_point() +
    labs(title = "INS vs upr_score by cell type", x = "upr_score", y = "INS_hk") +
    theme_classic() +
    Seurat::NoLegend() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"))


# geom_point with INK_hk and uprscore by cell type
ggplot(ins_data, aes(x = upr_score, y = INS_hk, color = protected, alpha = 0.2)) +
    geom_point() +
    labs(title = "INS_hk vs upr_score by cell type", x = "upr_score", y = "INS_hk") +
    theme_classic() +
    Seurat::NoLegend() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    facet_wrap(~cell_type_plus, scales = "fixed", ncol = 4)

# geom_point with INK_hk and uprscore by cell type
ggplot(ins_data, aes(x = upr_score, y = INS, color = protected, alpha = 0.2)) +
    geom_point() +
    labs(title = "INS vs upr_score by cell type", x = "upr_score", y = "INS_hk") +
    theme_classic() +
    Seurat::NoLegend() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    facet_wrap(~cell_type_plus, scales = "fixed", ncol = 4)
