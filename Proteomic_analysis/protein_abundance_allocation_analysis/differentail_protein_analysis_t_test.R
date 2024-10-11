#function to perform differential expression analysis and plot the volcano plots
##Srijith Sasikumar - (7-01-2024)

perform_differential_analysis_and_plot <- function(df1, df2, output_prefix) {
  # Create a folder using the output prefix
  output_folder <- paste0(getwd(), "/", output_prefix)
  dir.create(output_folder, showWarnings = FALSE)
  
  # Merge data
  merged_data <- merge(df1[, c(2, 6, 7, 3:5)], df2[, c(2, 3:5)], by = 'Entry', all = TRUE)
  
  # Make sure that all proteins occur in both sets
  merged_data_all <- na.omit(merged_data)
  
  # Perform p-value and log2 fold change calculation
  merged_data_DE <- merged_data_all
  merged_data_DE$p.value <- NA
  merged_data_DE$log2.FC <- NA
  
  for (i in 1:nrow(merged_data_DE)) {
    merged_data_DE[i, "p.value"] <- t.test(log2(merged_data_DE[i, 4:6]), log2(merged_data_DE[i, 7:9]))$p.value
    merged_data_DE[i, "log2.FC"] <- log2(rowMeans(merged_data_DE[i, 7:9]) / rowMeans(merged_data_DE[i, 4:6]))
  }
  
  # Adjust p-values
  merged_data_DE$FDR <- p.adjust(merged_data_DE$p.value, method = "fdr")
  
  # Save merged data to CSV within the created folder
  write.csv(merged_data_DE, file.path(output_folder, paste0(output_prefix, "_merged.csv")), row.names = FALSE)
  
  # Identify upregulated and downregulated proteins
  upreg <- merged_data_DE[merged_data_DE$log2.FC > 1 & merged_data_DE$p.value < 0.05,]
  downreg <- merged_data_DE[merged_data_DE$log2.FC < -1 & merged_data_DE$p.value < 0.05,]
  
  # Save upregulated and downregulated proteins separately within the created folder
  write.csv(upreg, file.path(output_folder, paste0(output_prefix, "_upregulated.csv")), row.names = FALSE)
  write.csv(downreg, file.path(output_folder, paste0(output_prefix, "_downregulated.csv")), row.names = FALSE)
  
  # Print the number of upregulated and downregulated proteins
  cat("Number of Upregulated Proteins:", nrow(upreg), "\n")
  cat("Number of Downregulated Proteins:", nrow(downreg), "\n")
  
  EnhancedVolcano(
    merged_data_DE,
    lab = merged_data_DE$Entry,
    x = 'log2.FC',
    y = 'p.value',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 2
  )
  ggsave(file.path(output_folder, paste0(output_prefix, "_volcano_plot.png")), width = 8, height = 6)
}

perform_differential_analysis_and_plot_outier_removed_df1 <- function(df1, df2, output_prefix) {
  # Create a folder using the output prefix
  output_folder <- paste0(getwd(), "/", output_prefix)
  dir.create(output_folder, showWarnings = FALSE)
  
  # Merge data
  merged_data <- merge(df1[, c(2, 6, 7, 3:4)], df2[, c(2, 3:5)], by = 'Entry', all = TRUE)
  
  # Make sure that all proteins occur in both sets
  merged_data_all <- na.omit(merged_data)
  
  # Perform p-value and log2 fold change calculation
  merged_data_DE <- merged_data_all
  merged_data_DE$p.value <- NA
  merged_data_DE$log2.FC <- NA
  
  for (i in 1:nrow(merged_data_DE)) {
    merged_data_DE[i, "p.value"] <- t.test(log2(merged_data_DE[i, 4:5]), log2(merged_data_DE[i, 6:8]))$p.value
    merged_data_DE[i, "log2.FC"] <- log2(rowMeans(merged_data_DE[i, 6:8]) / rowMeans(merged_data_DE[i, 4:5]))
  }
  
  # Adjust p-values
  merged_data_DE$FDR <- p.adjust(merged_data_DE$p.value, method = "fdr")
  
  # Save merged data to CSV within the created folder
  write.csv(merged_data_DE, file.path(output_folder, paste0(output_prefix, "_merged.csv")), row.names = FALSE)
  
  # Identify upregulated and downregulated proteins
  upreg <- merged_data_DE[merged_data_DE$log2.FC > 1 & merged_data_DE$p.value < 0.05,]
  downreg <- merged_data_DE[merged_data_DE$log2.FC < -1 & merged_data_DE$p.value < 0.05,]
  

  write.csv(upreg, file.path(output_folder, paste0(output_prefix, "_upregulated.csv")), row.names = FALSE)
  write.csv(downreg, file.path(output_folder, paste0(output_prefix, "_downregulated.csv")), row.names = FALSE)
  
  # Print the number of upregulated and downregulated proteins
  cat("Number of Upregulated Proteins:", nrow(upreg), "\n")
  cat("Number of Downregulated Proteins:", nrow(downreg), "\n")
  library(EnhancedVolcano)
  # Create an enhanced volcano plot and save it within the created folder
  EnhancedVolcano(
    merged_data_DE,
    lab = merged_data_DE$Entry,
    x = 'log2.FC',
    y = 'p.value',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 2
  )
  ggsave(file.path(output_folder, paste0(output_prefix, "_volcano_plot.png")), width = 8, height = 6)
}

#outlier_removed_df2
perform_differential_analysis_and_plot_outier_removed_df2 <- function(df1, df2, output_prefix) {
  # Create a folder using the output prefix
  output_folder <- paste0(getwd(), "/", output_prefix)
  dir.create(output_folder, showWarnings = FALSE)
  
  # Merge data
  merged_data <- merge(df1[, c(2, 6, 7, 3:5)], df2[, c(2, 3:4)], by = 'Entry', all = TRUE)
  
  # Make sure that all proteins occur in both sets
  merged_data_all <- na.omit(merged_data)
  
  # Perform p-value and log2 fold change calculation
  merged_data_DE <- merged_data_all
  merged_data_DE$p.value <- NA
  merged_data_DE$log2.FC <- NA
  
  for (i in 1:nrow(merged_data_DE)) {
    merged_data_DE[i, "p.value"] <- t.test(log2(merged_data_DE[i, 4:6]), log2(merged_data_DE[i, 7:8]))$p.value
    merged_data_DE[i, "log2.FC"] <- log2(rowMeans(merged_data_DE[i, 7:8]) / rowMeans(merged_data_DE[i, 4:6]))
  }
  
  # Adjust p-values
  merged_data_DE$FDR <- p.adjust(merged_data_DE$p.value, method = "fdr")
  
  # Save merged data to CSV within the created folder
  write.csv(merged_data_DE, file.path(output_folder, paste0(output_prefix, "_merged.csv")), row.names = FALSE)
  
  # Identify upregulated and downregulated proteins
  upreg <- merged_data_DE[merged_data_DE$log2.FC > 1 & merged_data_DE$p.value < 0.05,]
  downreg <- merged_data_DE[merged_data_DE$log2.FC < -1 & merged_data_DE$p.value < 0.05,]
  
  # Save upregulated and downregulated proteins separately within the created folder
  write.csv(upreg, file.path(output_folder, paste0(output_prefix, "_upregulated.csv")), row.names = FALSE)
  write.csv(downreg, file.path(output_folder, paste0(output_prefix, "_downregulated.csv")), row.names = FALSE)
  
  # Print the number of upregulated and downregulated proteins
  cat("Number of Upregulated Proteins:", nrow(upreg), "\n")
  cat("Number of Downregulated Proteins:", nrow(downreg), "\n")
  library(EnhancedVolcano)
  # Create an enhanced volcano plot and save it within the created folder
  EnhancedVolcano(
    merged_data_DE,
    lab = merged_data_DE$Entry,
    x = 'log2.FC',
    y = 'p.value',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 2
  )
  ggsave(file.path(output_folder, paste0(output_prefix, "_volcano_plot.png")), width = 8, height = 6)
}


setwd("D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis")
df_s_0 <- read.csv("proteomics_S_0.csv")
df_m_0 <- read.csv("proteomics_M0.csv")
df_t_0 <- read.csv("proteomics_T0.csv")
df_mt_0 <- read.csv("proteomics_MT0.csv")
df_S_230 <- read.csv("proteomics_S2h30.csv")
df_m_230 <- read.csv("proteomics_M2h30.csv")
df_t_230 <- read.csv("proteomics_T2h30.csv")
df_mt_230 <- read.csv("proteomics_MT2h30.csv")


perform_differential_analysis_and_plot (df_s_0, df_m_0, "S_M_0")
perform_differential_analysis_and_plot (df_s_0, df_t_0, "S_T_0")
perform_differential_analysis_and_plot (df_s_0, df_mt_0, "S_MT_0")


#using S_230 outlier removed function df1
perform_differential_analysis_and_plot_outier_removed_df1 (df_S_230, df_m_230, "S_M_2h30m")
perform_differential_analysis_and_plot_outier_removed_df1 (df_S_230, df_t_230, "S_T_2h30m")
perform_differential_analysis_and_plot_outier_removed_df1 (df_S_230, df_mt_230, "S_MT_2h30m")


#between timepoints outlier removed for S_230 df2
perform_differential_analysis_and_plot_outier_removed_df2 (df_s_0, df_S_230, "S_s_0_230")
#between timepoints
perform_differential_analysis_and_plot (df_m_0, df_m_230, "M_M_0_230")
perform_differential_analysis_and_plot (df_t_0, df_t_230, "T_T_0_230")
perform_differential_analysis_and_plot (df_mt_0, df_mt_230, "MT_MT_0_230")


