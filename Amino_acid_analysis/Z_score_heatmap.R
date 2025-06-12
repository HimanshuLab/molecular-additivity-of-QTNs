library(tidyverse)
library(ggplot2)

# 1. Data Import and Preprocessing
data <- read.csv("amino_acid_data.csv")

data <- data %>%
  rename_with(~gsub(" Results", "", .), ends_with("Results")) %>%
  mutate(Timepoint = ifelse(Timepoint == 230, 2.5, as.numeric(Timepoint)))


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


heatmap_data <- data_long_z %>%
  group_by(Strains, Timepoint, AminoAcid) %>%
  summarise(MeanZ = mean(RelativeQuantification, na.rm = TRUE), .groups = "drop") %>%
  mutate(Strains = factor(Strains, levels = c("MT", "S", "M", "T")))

# 3. Cluster amino acids using MT strain
mt_matrix <- heatmap_data %>%
  filter(Strains == "MT") %>%
  pivot_wider(names_from = Timepoint, values_from = MeanZ) %>%
  column_to_rownames("AminoAcid")

amino_order <- rownames(mt_matrix)[hclust(dist(mt_matrix))$order]


p <- heatmap_data %>%
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


hc <- hclust(dist(mt_matrix))

k <- 8
clusters <- cutree(hc, k = k)

cluster_df <- data.frame(AminoAcid = names(clusters), Cluster = clusters)
cluster_summary <- cluster_df %>%
  group_by(Cluster) %>%
  summarise(Count = n())

print(cluster_summary)


cluster_df <- data.frame(AminoAcid = names(clusters), Cluster = factor(clusters))


heatmap_clustered <- heatmap_data %>%
  left_join(cluster_df, by = "AminoAcid") %>%
  mutate(
    AminoAcid = factor(AminoAcid, levels = rownames(mt_matrix)[hc$order]),
    Cluster = factor(Cluster)
  )

# Plot heatmap with cluster facet wrap
hp <- ggplot(heatmap_clustered, aes(x = factor(Timepoint), y = AminoAcid, fill = MeanZ)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(Cluster ~ Strains, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#A50026", midpoint = 0,
    name = "Z-Score"
  ) +
  labs(
    x = "Timepoint (hours)",
    y = "Amino Acids",
    title = paste("Heatmap with Amino Acid Clusters (k =", k, ")")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = element_blank()
  )


ggsave('D:/mauscript2/revision/R_codes_data/figures/z_score_heatmap.png',hp, dpi = 600, height = 6, width = 6)
ggsave('D:/mauscript2/revision/R_codes_data/figures/z_score_heatmap.svg',hp, dpi = 600, height = 6, width = 6)




