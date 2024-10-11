#function to perform summed mass allocation analysis on GO_simmer list
#srijith sasikumar 01-02-2024
library(ggplot2)
library(data.table)
library(ggplot2)
library(gridExtra)
setwd("D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis")
go_list <- read.csv("GO_simmer_list.csv")
unique(go_list$GO.term)
gene_GO_list <- go_list[, c(2,4)]
# Rename columns to unify names used
names(gene_GO_list)[1:2]<-c("Gene","GO")
orf_gene_conversion <- read.csv("orf_gene_conversion.csv")
gene_uniprot <- read.csv("protein_ID_kDA_maps.csv")

library(ggplot2)
library(dplyr)

perform_proteomic_analysis_without_outlier <- function(file1, file2, type1, type2, output_file, csv_file) {
  # Load data
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  
  # Merge data frames with gene information
  df1_with_gene <- merge(df1, gene_uniprot[, c(2, 8)], by = "Entry")
  df2_with_gene <- merge(df2, gene_uniprot[, c(2, 8)], by = "Entry")
  
  # Combine strains and remove NAs
  combined_strains <- merge(df1_with_gene[, c(3, 4, 5, 8)], df2_with_gene[, c(3, 4, 5, 8)], by = 'Gene')
  combined_strains <- na.omit(combined_strains)
  combined_strains_GO <- merge(combined_strains, gene_GO_list, by = 'Gene')
  
  # Perform analysis for type1
  result <- combined_strains_GO %>%
    group_by(GO) %>%
    summarise(across(c(2, 3, 4), sum))
  
  result$Mean <- apply(result[, 2:4], 1, mean)
  result$Stdev <- apply(result[, 2:4], 1, sd)
  
  # Perform analysis for type2
  result_1 <- combined_strains_GO %>%
    group_by(GO) %>%
    summarise(across(c(5, 6), sum))
  
  result_1$Mean <- apply(result_1[, 2:3], 1, mean)
  result_1$Stdev <- apply(result_1[, 2:3], 1, sd)
  
  # Merge results
  result_combined <- merge(result[, c(1:4)], result_1[, c(1:3)])
  
  # Perform t-test and add p-value to result_combined
  result_combined$p.value <- apply(result_combined[, c(2:6)], 1, function(row) {
    if (all(!is.na(row)) && all(sapply(row, is.numeric))) {
      t.test(log2(row[1:3]), log2(row[4:5]),alternative='t',paired=FALSE)$p.value
    } else {
      NA
    }
  })
  write.csv(result_combined, csv_file)
  # Add types
  result$type <- type1
  result_1$type <- type2
  
  # Set row names
  rownames(result) <- result$GO
  rownames(result_1) <- result_1$GO
  rownames(result_combined) <- result_combined$GO
  
  # Filter significant results
  results_sig <- result_combined[result_combined$p.value < 0.05, ]
  result <- result[results_sig$GO, ]
  result_1 <- result_1[results_sig$GO, ]
  
  # Combine results for plot
  result_for_plot <- rbind(result[, c(1, 5, 6, 7)], result_1[, c(1, 4, 5, 6)])
  result_for_plot$Mean <- result_for_plot$Mean * 100
  result_for_plot$Stdev <- result_for_plot$Stdev * 100
  
  # Arrange data frame in descending order of Mean
  filtered_result_for_plot <- result_for_plot[order(result_for_plot$Mean), ]
  
  # Create the ggplot
  plot <- ggplot(filtered_result_for_plot, aes(reorder(GO, -Mean), Mean, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), colour = "black") +
    geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev),
                  position = position_dodge(0.9), width = 0.2) +
    theme_bw(base_family = "serif") +
    theme(text = element_text(colour = "black", face = "bold"),
          axis.text = element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = "Proteomic Analysis Plot",
         x = "",
         y = "Mass fraction of proteome (%)",
         fill = "") +
    scale_y_continuous(limits = c(0, 30)) +
    scale_fill_manual(values = c("blue", "orange")) +
    theme(legend.position = "bottom") +
    coord_flip()
  ggsave(output_file, plot, width = 7, height = 8, dpi = 300)
  
  print(plot)
}


perform_proteomic_analysis_without_outlier("proteomics_MT2h30.csv", "proteomics_S2h30.csv", "MT", "S", "S2h30_MT2h30_go_simmer.png", "S2h30_MT2h30_go_simmer.csv")
perform_proteomic_analysis_without_outlier("proteomics_M2h30.csv", "proteomics_S2h30.csv", "M", "S", "S2h30_M2h30_go_simmer.png","S2h30_M2h30_go_simmer.csv")
perform_proteomic_analysis_without_outlier("proteomics_T2h30.csv", "proteomics_S2h30.csv", "T", "S", "S2h30_T2h30_go_simmer.png","S2h30_T2h30_go_simmer.csv")

#
perform_proteomic_analysis<- function(file1, file2, type1, type2, output_file, csv_file) {
  # Load data
  df1 <- read.csv(file1)
  df2 <- read.csv(file2)
  
  # Merge data frames with gene information
  df1_with_gene <- merge(df1, gene_uniprot[, c(2, 8)], by = "Entry")
  df2_with_gene <- merge(df2, gene_uniprot[, c(2, 8)], by = "Entry")
  
  # Combine strains and remove NAs
  combined_strains <- merge(df1_with_gene[, c(3, 4, 5, 8)], df2_with_gene[, c(3, 4, 5, 8)], by = 'Gene')
  combined_strains <- na.omit(combined_strains)
  combined_strains_GO <- merge(combined_strains, gene_GO_list, by = 'Gene')
  
  # Perform analysis for type1
  result <- combined_strains_GO %>%
    group_by(GO) %>%
    summarise(across(c(2, 3, 4), sum))
  
  result$Mean <- apply(result[, 2:4], 1, mean)
  result$Stdev <- apply(result[, 2:4], 1, sd)
  
  # Perform analysis for type2
  result_1 <- combined_strains_GO %>%
    group_by(GO) %>%
    summarise(across(c(5, 6,7), sum))
  
  result_1$Mean <- apply(result_1[, 2:4], 1, mean)
  result_1$Stdev <- apply(result_1[, 2:4], 1, sd)
  
  # Merge results
  result_combined <- merge(result[, c(1:4)], result_1[, c(1:4)])
  
  # Perform t-test and add p-value to result_combined
  result_combined$p.value <- apply(result_combined[, c(2:7)], 1, function(row) {
    if (all(!is.na(row)) && all(sapply(row, is.numeric))) {
      t.test(log2(row[1:3]), log2(row[4:6]),alternative='t',paired=FALSE)$p.value
    } else {
      NA
    }
  })
  write.csv(result_combined, csv_file)
  # Add types
  result$type <- type1
  result_1$type <- type2
  
  # Set row names
  rownames(result) <- result$GO
  rownames(result_1) <- result_1$GO
  rownames(result_combined) <- result_combined$GO
  
  # Filter significant results
  results_sig <- result_combined[result_combined$p.value < 0.05, ]
  result <- result[results_sig$GO, ]
  result_1 <- result_1[results_sig$GO, ]
  
  # Combine results for plot
  result_for_plot <- rbind(result[, c(1, 5, 6, 7)], result_1[, c(1, 5, 6, 7)])
  result_for_plot$Mean <- result_for_plot$Mean * 100
  result_for_plot$Stdev <- result_for_plot$Stdev * 100
  
  # Arrange data frame in descending order of Mean
  filtered_result_for_plot <- result_for_plot[order(result_for_plot$Mean), ]
  
  # Create the ggplot
  plot <- ggplot(filtered_result_for_plot, aes(reorder(GO, -Mean), Mean, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), colour = "black") +
    geom_errorbar(aes(ymin = Mean - Stdev, ymax = Mean + Stdev),
                  position = position_dodge(0.9), width = 0.2) +
    theme_bw(base_family = "serif") +
    theme(text = element_text(colour = "black", face = "bold"),
          axis.text = element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = "Proteomic Analysis Plot",
         x = "",
         y = "Mass fraction of proteome (%)",
         fill = "") +
    scale_y_continuous(limits = c(0, 30)) +
    scale_fill_manual(values = c("blue", "orange")) +
    theme(legend.position = "bottom") +
    coord_flip()
  ggsave(output_file, plot, width = 7, height = 8, dpi = 300)
  # Display the plot
  print(plot)
  
}

# Example usage
perform_proteomic_analysis("proteomics_MT0.csv", "proteomics_S_0.csv", "MT", "S", "S0_MT0_go_simmer.png", "S0_MT0_go_simmer.csv")
perform_proteomic_analysis("proteomics_M0.csv", "proteomics_S_0.csv", "M", "S", "S0_M0_go_simmer.png","S0_M0_go_simmer.csv")
perform_proteomic_analysis("proteomics_T0.csv", "proteomics_S_0.csv", "T", "S", "S0_T0_go_simmer.png","S0_T0_go_simmer.csv")

