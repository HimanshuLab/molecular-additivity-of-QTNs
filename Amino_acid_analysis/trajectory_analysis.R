# Install necessary packages if not already installed
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(factoextra)) install.packages("factoextra")

# Load libraries
library(tidyverse)
library(ggplot2)
library(factoextra)

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
# Reshape to wide format
wide_data <- data_long %>%
  select(Strains, Timepoint, Replicates, AminoAcid, RelativeQuantification) %>%
  pivot_wider(names_from = AminoAcid, values_from = RelativeQuantification)

meta_cols <- wide_data %>% select(Strains, Timepoint, Replicates)

# Z-score normalize across all samples for each amino acid
scaled_mat <- scale(wide_data %>% select(-Strains, -Timepoint, -Replicates))

# Combine metadata and normalized matrix
zscore_df <- cbind(meta_cols, as.data.frame(scaled_mat))
data_long_z <- zscore_df %>%
  pivot_longer(
    cols = everything()[-(1:3)],
    names_to = "AminoAcid",
    values_to = "RelativeQuantification"
  )

# Define strain colors
strain_colors <- c("S" = "#D55E00", "M" = "#0072B2", "T" = "#009E73", "MT" = "#CC79A7")
data_long_z$Timepoint <- as.numeric(as.character(data_long_z$Timepoint))

# Summary statistics
all_aa <- data_long_z %>%
  group_by(Strains, Timepoint, AminoAcid) %>%
  summarize(
    Mean = mean(RelativeQuantification), 
    SE = sd(RelativeQuantification) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Timepoint, y = Mean, color = Strains, group = Strains)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  scale_color_manual(values = strain_colors) +
  facet_wrap(~AminoAcid, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Z-Score Normalized Amino Acid Time-Course",
    x = "Timepoint (hours)",
    y = "Mean Z-Score"
  ) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.grid.major = element_line(color = "gray90")
  )

print(all_aa)

# Selected amino acids
selected_aa_list <- c("ALANINE", "HISTIDINE", "GLUTAMINE", "LEUCINE", "ARGININE", "LYSINE")
not_utilized_list <- c("TYROSINE", "TRYPTOPHAN", "PHENYLALANINE", "ASPARAGINE", "VALINE", "GLUTAMIC.ACID", "THREONINE", "HOMOSERINE")

selected_aa_plot_data <- data_long_z %>%
  filter(AminoAcid %in% selected_aa_list) %>%
  group_by(Strains, Timepoint, AminoAcid) %>%
  summarize(
    Mean = mean(RelativeQuantification),
    SD = sd(RelativeQuantification),
    .groups = "drop"
  )
selected_aa_plot_data_utilized <- data_long_z %>%
  filter(AminoAcid %in% not_utilized_list) %>%
  group_by(Strains, Timepoint, AminoAcid) %>%
  summarize(
    Mean = mean(RelativeQuantification),
    SD = sd(RelativeQuantification),
    .groups = "drop"
  )

# Rename and reorder strain labels
recode_strains <- function(df) {
  df %>%
    mutate(
      Strains = recode(Strains, "MT" = "MMTT", "S" = "SS", "M" = "MM", "T" = "TT"),
      Strains = factor(Strains, levels = c("MMTT", "SS", "MM", "TT"))
    )
}

selected_aa_plot_data <- recode_strains(selected_aa_plot_data)
selected_aa_plot_data_utilized <- recode_strains(selected_aa_plot_data_utilized)

# Color map for selected amino acids
aa_colors <- c(
  "HISTIDINE" = "#E69F00",
  "ARGININE"  = "#56B4E9",
  "GLUTAMINE" = "#009E73",
  "LYSINE"    = "#D55E00",
  "ALANINE"   = "#0072B2",
  "LEUCINE"   = "#CC79A7"
)

# Plot for selected amino acids
p <- ggplot(selected_aa_plot_data, aes(x = Timepoint, y = Mean, color = AminoAcid, group = AminoAcid)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
  facet_wrap(~Strains, nrow = 2) +
  scale_color_manual(values = aa_colors) +
  labs(
    x = "Timepoint (hours)",
    y = "Mean Z-Score",
    color = "Amino Acid"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# Plot for not-utilized amino acids
q <- ggplot(selected_aa_plot_data_utilized, aes(x = Timepoint, y = Mean, color = AminoAcid, group = AminoAcid)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
  facet_wrap(~Strains, nrow = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Timepoint (hours)",
    y = "Mean Z-Score",
    color = "Amino Acid"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )
p
# Save plots
ggsave("D:/mauscript2/revision/R_codes_data/figures/amino_acid_not_utilized_sd.png", p, dpi = 600, width = 5, height = 4)
ggsave("D:/mauscript2/revision/R_codes_data/figures/amino_acid_utilized_sd.png", q, dpi = 600, width = 5, height = 4)
