library(dplyr)
library(tidyr)
library(ComplexUpset)

targets::tar_load(pancdb_metadata_gt)

genotypes <- pancdb_metadata_gt |>
    dplyr::select(sample_id, rs3842752_consensus, rs3842753_consensus, rs689_consensus)

library(dplyr)
library(ComplexUpset)

# Step 1: Transform the genotype data into binary indicators
genotypes_transformed <- genotypes %>%
    pivot_longer(cols = -sample_id, names_to = "locus", values_to = "genotype") %>%
    mutate(genotype_presence = 1) %>%
    pivot_wider(names_from = c(locus, genotype), values_from = genotype_presence, values_fill = list(genotype_presence = 0))

# Step 2: Create the UpSet plot
ComplexUpset::upset(
    genotypes_transformed,
    intersect = names(genotypes_transformed)[grep("consensus", names(genotypes_transformed))],
    name = "locus_genotype",
    width_ratio = 0.1
)
