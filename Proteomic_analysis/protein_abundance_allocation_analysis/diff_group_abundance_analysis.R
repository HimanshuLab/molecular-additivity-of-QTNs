#function to perform differentail protein abundance based on groups

#load libraries
library(ggplot2)
library(data.table)
library(gridExtra)
setwd("D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis")
df_s_0 <- read.csv("proteomics_S_0.csv")
df_m_0 <- read.csv("proteomics_M0.csv")
df_t_0 <- read.csv("proteomics_T0.csv")
df_mt_0 <- read.csv("proteomics_MT0.csv")
df_S_230 <- read.csv("proteomics_S2h30.csv")
df_m_230 <- read.csv("proteomics_M2h30.csv")
df_t_230 <- read.csv("proteomics_T2h30.csv")
df_mt_230 <- read.csv("proteomics_MT2h30.csv")

#load the group data
Translation<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Translation_gene_names_from_elife.csv",
                      sep=";")
Transcription<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Transcription_gene_names_from_elife.csv",
                        sep=";")
Stress<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Stress_gene_names_from_elife.csv",
                 sep=";")
Ribosomes<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Ribosomes_gene_names_from_elife.csv",
                    sep=";")
Nucleotides<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Nucleotides_gene_names_from_elife.csv",
                      sep=";")
Mitochondria<-read.csv("mitochondria_genes.csv")
Lipids<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Lipids_gene_names_from_elife.csv",
                 sep=";")
Glycolysis<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Glycolysis_gene_names_from_elife.csv",
                     sep=";")
Energy<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Energy_gene_names_from_elife.csv",
                 sep=";")
ER<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/ER_gene_names_from_elife.csv",
             sep=";")
Chaperones<-read.csv(file = "D:/S_MMTT_analysis/proteomic_analysis/protein_allocation_analysis/protein_group_ananalysis/complementaryData/Chaperons_gene_names_from_elife.csv",
                     sep=";")
AA_metabolism<-read.csv("amino_acid_genes.csv")
ppp_metabolism <- read.csv("ppp_genes.csv")
TCA_glyoxylate_cycle <- read.csv("TCA_glyxo_genes.csv")
respiration <- read.csv("respiration_genes.csv")
Mitochondrial_translation <- read.csv("mitochondial_translation.csv")
# Merge processes
Translation2<-Translation
Transcription2<-Transcription
Stress2<-Stress
Ribosomes2<-Ribosomes
Nucleotides2<-Nucleotides
Mitochondria2<-Mitochondria
Lipids2<-Lipids
Glycolysis2<-Glycolysis
Energy2<-Energy
ER2<-ER
Chaperones2<-Chaperones
AA_metabolism2<-AA_metabolism
ppp_metabolism2 <- ppp_metabolism
TCA_glyoxylate_cycle2 <- TCA_glyoxylate_cycle
respiration2<- respiration
Mitochondrial_translation2 <- Mitochondrial_translation

Translation2$Process<-"Translation"
Transcription2$Process<-"Transcription"
Stress2$Process<-"Stress"
Ribosomes2$Process<-"Ribosomes"
Nucleotides2$Process<-"Nucleotides"
Mitochondria2$Process<-"Mitochondria"
Lipids2$Process<-"Lipids"
Glycolysis2$Process<-"Glycolysis"
Energy2$Process<-"Energy"
ER2$Process<-"ER"
Chaperones2$Process<-"Chaperones"
AA_metabolism2$Process<-"AA_metabolism"
ppp_metabolism2$Process <- "PPP_metabolism"
TCA_glyoxylate_cycle2$Process <- "TCA_glyoxylate_cycle"
Chaperones2<-Chaperones2[,c(1,2,4)]
AA_metabolism2<-AA_metabolism2[,c(1,2,4)]
respiration2$Process <- "Respiration"
Mitochondrial_translation2$Process <- "Mito_translation"
names(Ribosomes2)[1]<-"ORF"
names(Glycolysis2)[1]<-"ORF"
names(ppp_metabolism2)[1:2] <- c("ORF","uniprot")

eLife_all<-rbind(Translation2, Transcription2, Stress2, Ribosomes2,
                 Nucleotides2, Mitochondria2, Lipids2, Glycolysis2, Energy2,
                 ER2, Chaperones2, AA_metabolism2, ppp_metabolism2, TCA_glyoxylate_cycle2, respiration2,Mitochondrial_translation2)
unique(eLife_all$Process)

# Rename columns to unify names used
names(eLife_all)[1:2]<-c("Gene","Entry")

#add a column type 
df_s_0$Type<-"S-0"
df_S_230$Type<-"S-2h30m" 
df_m_0$Type<-"M-0"
df_m_230$Type<-"M-2h30m"
df_t_0$Type<-"T-0"
df_t_230$Type<-"T-2h30m"
df_mt_0$Type<-"MT-0"
df_mt_230$Type<-"MT-2h30m"


library(ggplot2)
library(data.table)

diff_group_abundance <- function(df1, df2,title) {
  
  # Rename rep_columns
  names(df1)[3:5] <-  c("rep1", "rep2", "rep3")
  names(df2)[3:5] <-  c("rep1", "rep2", "rep3")
  
  # Combine dataframes
  df_combined <- rbind(df1, df2)
  df_combined <- na.omit(df_combined)
  
  # Determine mean for reps
  df_combined$Mean <- apply(df_combined[, 3:5], 1, mean)
  
  # Determine stdev for reps
  df_combined$Stdev <- apply(df_combined[, 3:5], 1, sd)
  
  # Summarize
  df_combined_subset <- merge(df_combined, eLife_all, by = "Entry")
  df_sum_means <- as.data.frame(as.data.table(df_combined_subset)[, sum(Mean, na.rm = TRUE), by = .(Type, Process)])
  df_sum_stdev <- as.data.frame(as.data.table(df_combined_subset)[, sd(Mean, na.rm = TRUE), by = .(Type, Process)])
  
  names(df_sum_means)[3] <- "Sum"
  names(df_sum_stdev)[3] <- "Stdev"
  
  df_combined_final <- merge(df_sum_means, df_sum_stdev)
  
  # Convert to %
  df_combined_final$Sum <- df_combined_final$Sum * 100
  df_combined_final$Stdev <- df_combined_final$Stdev * 100
  
  # Plotting
 
  plot <- ggplot(df_combined_final, aes(Process, Sum, fill=Type)) +
    geom_bar(stat="identity", position=position_dodge(0.9), colour="black") +
    geom_errorbar(aes(ymin=Sum-Stdev, ymax=Sum+Stdev),
                  position=position_dodge(0.9), width=0.2)+
    theme_bw(base_family= "serif") +
    theme(text = element_text(colour = "black",face="bold"),
          axis.text= element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = title,
         x="",
         y="Mass fraction of proteome (%)",
         fill="")  +
    scale_y_continuous(limits=c(0,50)) +
    scale_fill_manual(values=c("orange", "tan4")) +
    theme(legend.position = "bottom") +
    coord_flip() 
    ggsave(paste0(title, "_diff_group_abundance.png"),dpi = 600, width = 4, height = 5)
    return(plot)
}


