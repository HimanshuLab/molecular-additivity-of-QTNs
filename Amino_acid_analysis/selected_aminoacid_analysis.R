if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(factoextra)) install.packages("factoextra")

library(tidyverse)
library(ggplot2)
library(factoextra)

# 1. Data Import and Preprocessing
data <- read.csv("amino_acid_data.csv")

# Fix column names and timepoint
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
# Reshape to wide format
wide_data <- data_long %>%
  select(Strains, Timepoint, Replicates, AminoAcid, RelativeQuantification) %>%
  pivot_wider(names_from = AminoAcid, values_from = RelativeQuantification)

# Keep metadata
meta_cols <- wide_data %>% select(Strains, Timepoint, Replicates)

# Z-score normalize across all samples for each amino acid
scaled_mat <- scale(wide_data %>% select(-Strains, -Timepoint, -Replicates))

# Combine metadata and normalized matrix
zscore_df <- cbind(meta_cols, as.data.frame(scaled_mat))

# Convert back to long format
data_long_z <- zscore_df %>%
  pivot_longer(
    cols = c(GLYCINE:TRYPTOPHAN),
    names_to = "AminoAcid",
    values_to = "RelativeQuantification"
  )

# Define strain colors
strain_colors <- c("S" = "#D55E00", "M" = "#0072B2", "T" = "#009E73", "MT" = "#CC79A7")

# Ensure Timepoint is numeric
data_long_z$Timepoint <- as.numeric(as.character(data_long_z$Timepoint))

# Summarize and plot
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

# Display the plot
all_aa



# Rename strain labels
selected_aa_plot_data$Strains <- recode(selected_aa_plot_data$Strains,
                                        "MT" = "MMTT", "S" = "SS", "M" = "MM", "T" = "TT")
selected_aa_plot_data_utilized$Strains <- recode(selected_aa_plot_data_utilized$Strains,
                                                 "MT" = "MMTT", "S" = "SS", "M" = "MM", "T" = "TT")

# Set factor levels in the new order
selected_aa_plot_data$Strains <- factor(selected_aa_plot_data$Strains, levels = c("MMTT", "SS", "MM", "TT"))
selected_aa_plot_data_utilized$Strains <- factor(selected_aa_plot_data_utilized$Strains, levels = c("MMTT", "SS", "MM", "TT"))

# Plot 1
p <- ggplot(selected_aa_plot_data, aes(x = Timepoint, y = Mean, color = AminoAcid, group = AminoAcid)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.15) +
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

# Define color map for selected amino acids
aa_colors_selected <- c(
  "HISTIDINE" = "#E69F00",
  "ARGININE"  = "#56B4E9",
  "GLUTAMINE" = "#009E73",
  "LYSINE"    = "#D55E00",
  "ALANINE"   = "#0072B2"
)

# Plot 2
q <- ggplot(selected_aa_plot_data_utilized, aes(x = Timepoint, y = Mean, color = AminoAcid, group = AminoAcid)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.15) +
  facet_wrap(~Strains, nrow = 2) +
  scale_color_manual(values = aa_colors_selected) +
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

# Display plot
p

ggsave("amino_acid_not_utilized.png", p , dpi = 600, width = 5, height = 4 )
ggsave("amino_acid_utilized.png", q , dpi = 600, width = 5, height = 4 )



###########3identify differntailly regulated aminoacids in MT






library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
# Ensure time is numeric
data_long_z$Timepoint <- as.numeric(as.character(data_long_z$Timepoint))

# Subset for MT strain
mt_data <- data_long_z %>% filter(Strains == "MT")

# Function to perform two t-tests and return p-values and means
test_pattern <- function(df) {
  time0 <- df %>% filter(Timepoint == 0) %>% pull(RelativeQuantification)
  time2.5 <- df %>% filter(Timepoint == 2.5) %>% pull(RelativeQuantification)
  time8 <- df %>% filter(Timepoint == 8) %>% pull(RelativeQuantification)
  
  if (length(time0) >= 2 & length(time2.5) >= 2 & length(time8) >= 2) {
    p_up <- t.test(time0, time2.5, paired = FALSE)$p.value
    p_down <- t.test(time2.5, time8, paired = FALSE)$p.value
    
    mean0 <- mean(time0)
    mean2.5 <- mean(time2.5)
    mean8 <- mean(time8)
    
    return(tibble(
      p_up = p_up,
      p_down = p_down,
      mean_0 = mean0,
      mean_2.5 = mean2.5,
      mean_8 = mean8
    ))
  } else {
    return(tibble(p_up = NA, p_down = NA, mean_0 = NA, mean_2.5 = NA, mean_8 = NA))
  }
}

# Apply test per amino acid
mt_results <- mt_data %>%
  group_by(AminoAcid) %>%
  group_modify(~test_pattern(.x)) %>%
  ungroup()

# Filter for pattern: significant increase (p < 0.05), significant decrease (p < 0.05)
# and also directionality (increase then decrease)
filtered_mt <- mt_results %>%
  filter(p_up < 0.05, p_down < 0.05, mean_0 < mean_2.5, mean_2.5 > mean_8)

# View filtered results
filtered_mt
