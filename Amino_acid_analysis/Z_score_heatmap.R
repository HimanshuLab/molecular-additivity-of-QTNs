library(tidyverse)
library(ggplot2)

# 1. Data Import and Preprocessing
data <- read.csv("amino_acid_data.csv")

data <- data %>%
  rename_with(~gsub(" Results", "", .), ends_with("Results")) %>%
  mutate(Timepoint = ifelse(Timepoint == 230, 2.5, as.numeric(Timepoint)))

# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = c(GLYCINE:TRYPTOPHAN),
    names_to = "AminoAcid",
    values_to = "RelativeQuantification"
  ) %>%
  mutate(
    Strains = factor(Strains),
    Timepoint = factor(Timepoint),
    Replicates = factor(Replicates)
  )

# 2. Z-Score Normalization
wide_data <- data_long %>%
  select(Strains, Timepoint, Replicates, AminoAcid, RelativeQuantification) %>%
  pivot_wider(names_from = AminoAcid, values_from = RelativeQuantification)

meta_cols <- wide_data %>% select(Strains, Timepoint, Replicates)
scaled_mat <- scale(wide_data %>% select(-Strains, -Timepoint, -Replicates))
zscore_df <- cbind(meta_cols, as.data.frame(scaled_mat))

data_long_z <- zscore_df %>%
  pivot_longer(
    cols = c(GLYCINE:TRYPTOPHAN),
    names_to = "AminoAcid",
    values_to = "RelativeQuantification"
  ) %>%
  mutate(Timepoint = as.numeric(as.character(Timepoint)))

# 3. Compute mean Z-scores per amino acid, strain, and timepoint
heatmap_data <- data_long_z %>%
  group_by(Strains, Timepoint, AminoAcid) %>%
  summarise(MeanZ = mean(RelativeQuantification, na.rm = TRUE), .groups = "drop") %>%
  mutate(Strains = factor(Strains, levels = c("MT", "S", "M", "T")))

# 4. Cluster amino acids using MT strain
mt_matrix <- heatmap_data %>%
  filter(Strains == "MT") %>%
  pivot_wider(names_from = Timepoint, values_from = MeanZ) %>%
  column_to_rownames("AminoAcid")

amino_order <- rownames(mt_matrix)[hclust(dist(mt_matrix))$order]

# 5. Plot publication-quality heatmap
heatmap_data %>%
  mutate(AminoAcid = factor(AminoAcid, levels = amino_order)) %>%
  ggplot(aes(x = factor(Timepoint), y = AminoAcid, fill = MeanZ)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~Strains, nrow = 1) +
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#A50026", midpoint = 0,
    name = "Z-Score"
  ) +
  labs(
    x = "Timepoint (hours)",
    y = "Amino Acids"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")
  )
