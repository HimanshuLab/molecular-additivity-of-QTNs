orf_gene_conversion <- read.csv("orf_gene_conversion.csv")
gene_uniprot <- read.csv("protein_ID_kDA_maps.csv")
proteomics_fmol <- read.csv("protein_fmol_microgram.csv")
rownames(proteomics_fmol) <- proteomics_fmol$X
proteomics_fmol <- proteomics_fmol[-1]
#names(proteomics_fmol)[1]<-paste("Entry")
#proteomics_fmol <- na.omit(proteomics_fmol)
proteomics_fmol_with_gene <-  proteomics_fmol[rowSums(proteomics_fmol) >= 5, ]
proteomics_fmol_with_gene <- na.omit(proteomics_fmol_with_gene)
library(tibble)
proteomics_fmol_with_gene  <- rownames_to_column(proteomics_fmol_with_gene, var = "Entry")
proteomics_fmol_with_gene <- merge(proteomics_fmol_with_gene,gene_uniprot[,c(2,8)], by="Entry")
#write.csv(proteomics_fmol_with_gene, "fmol_with_gene.csv")
#load transcritomics data
tpm_s <- read.csv("tpm_mean_s.csv")
tpm_s <- tpm_s[,-1]
tpm_s <- as.matrix(tpm_s)
tpm_s<- t(tpm_s)
colnames(tpm_s) <- tpm_s[1,]
tpm_s_final <- tpm_s[-1,]
tpm_s_final <- as.data.frame(tpm_s_final)
tpm_s_final[] <- lapply(tpm_s_final, as.numeric)
tpm_s_final <- tpm_s_final[rowSums(tpm_s_final) >= 10, ]
colnames(tpm_s_final) <- c("T0","T1","T2","T3","T4","T5","T6","T7","T8")
library(ggplot2)
#tpm_s_final$T0 <- round(as.numeric(tpm_s_final$T0))
library(tibble)
tpm_s_final <- rownames_to_column(tpm_s_final, var = "ORF")
tpm_s_final_gene <- merge(tpm_s_final,orf_gene_conversion[,1:2] , by="ORF")
tpm_fmol_merged <- merge(tpm_s_final_gene,proteomics_fmol_with_gene, by="Gene")
#calculate correlation for s_0
ID <- tpm_fmol_merged$Gene
fmol_average_s_0 <- log2(rowMeans(tpm_fmol_merged[,13:15],na.rm=TRUE))
tpm_average_s_0  <- log2(tpm_fmol_merged[,3])
tpm_fmol_averge_s_0 <- data.frame(ID = ID, fmol_average =fmol_average_s_0,tpm_averge =tpm_average_s_0 )
tpm_fmol_averge_s_0 <- na.omit(tpm_fmol_averge_s_0)




#test corelation 
cor(tpm_fmol_averge_s_0$fmol_average,tpm_fmol_averge_s_0$tpm_averge,use = "pairwise.complete.obs", method = "pearson")


library(ggplot2)


plot <- ggplot(tpm_fmol_averge_s_0, aes(x = tpm_averge, y = fmol_average)) +
  geom_point(color = "grey", size = 2, alpha = 1, shape = 16) + 
  geom_smooth(method = "lm", se = FALSE, color = "black")   

# Calculate correlation coefficient and slope
cor_coef <- cor(tpm_fmol_averge_s_0$tpm_averge, tpm_fmol_averge_s_0$fmol_average, use = "pairwise.complete.obs", method = "pearson")
slope <- coef(lm(fmol_average ~ tpm_averge, data = tpm_fmol_averge_s_0))[2]


plot + annotate("text", x = min(tpm_fmol_averge_s_0$tpm_averge), y = max(tpm_fmol_averge_s_0$fmol_average),
                label = paste("Correlation =", round(cor_coef, 3), "\nSlope =", round(slope, 3)),
                hjust = 0, vjust = 1, col = "black", size = 5) +
  
  
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  
        axis.line = element_line(color = "black"), axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
  
  xlab("log2(mRNA abundance (TPM))") + ylab("log2(protein abundance (fmol/ug))") + ggtitle("Pearson correlation s_0h")
ggsave("pearson_corelation_s_0.png", dpi = 600, height = 4, width = 4)
ggsave("D:/mauscript2/figures/figures_svg/pearson_corelation_s_0.png", dpi = 600, height = 4, width = 4)



#s-2h30m

ID <- tpm_fmol_merged$Gene
fmol_average_s_2h30 <- log2(rowMeans(tpm_fmol_merged[,25:26],na.rm=TRUE))
tpm_average_s_2h30   <- log2(tpm_fmol_merged[,8])
tpm_fmol_averge_s_2h30 <- data.frame(ID = ID, fmol_average =fmol_average_s_2h30,tpm_averge =tpm_average_s_2h30)
tpm_fmol_averge_s_2h30<- na.omit(tpm_fmol_averge_s_2h30)

#test corelation 
cor(tpm_fmol_averge_s_2h30$fmol_average,tpm_fmol_averge_s_2h30$tpm_averge,use = "pairwise.complete.obs", method = "pearson")
#0.3948377
# Assuming tpm_fmol_averge_s_0 is your data frame
library(ggplot2)
plot <- ggplot(tpm_fmol_averge_s_2h30, aes(x = tpm_averge, y = fmol_average)) +
  geom_point(color = "grey", size = 2, alpha = 1, shape = 16) +  
  geom_smooth(method = "lm", se = FALSE, color = "black")   

# Calculate correlation coefficient and slope
cor_result <- cor.test(tpm_fmol_averge_s_2h30$tpm_averge, tpm_fmol_averge_s_2h30$fmol_average, use = "pairwise.complete.obs", method = "pearson")
slope <- coef(lm(fmol_average ~ tpm_averge, data = tpm_fmol_averge_s_2h30))[2]
lm_model <- lm(fmol_average ~ tpm_averge, data = tpm_fmol_averge_s_2h30)
slope_pvalue <- summary(lm_model)$coefficients[2, 4]

plot +   annotate("text", x = min(tpm_fmol_averge_s_2h30$tpm_averge), y = max(tpm_fmol_averge_s_2h30$fmol_average),
                  label = paste("Correlation =", round(cor_result$estimate, 3),
                                #"\nP-value =", format(cor_result$p.value, scientific = TRUE, digits = 2),
                                "\nSlope =", round(lm_model$coefficients[2], 3)),
                                #"\nSlope P-value =", format(slope_pvalue, scientific = TRUE, digits = 2)),
                  hjust = 0, vjust = 1, col = "black", size = 5) +
  

  theme_minimal() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.line = element_line(color = "black"), axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
  xlab("log2(mRNA abundance (TPM))") + ylab("log2(protein abundance (fmol/ug))") + ggtitle("Pearson correlation s_2h30m") + coord_cartesian(ylim = c(-10, 10))
ggsave("D:/mauscript2/figures/figures_svg/pearson_corelation_s_2h30m.png", dpi = 600, height = 4, width = 4)


#look for whether different gene classes displayed specific expression patterns
pathways <- read.csv("GO_simmer_list.csv")


result_df <- data.frame(GO.term = character(), Number_of_genes = integer(), Correlation = numeric(), Slope = numeric(), Slope_PValue = numeric(), Corr_pvalue = character(), stringsAsFactors = FALSE)
fmol_tpm_merged_pathways <- merge(tpm_fmol_merged, pathways ,by.x = "Gene",by.y = "GENE")
unique_terms <- unique(fmol_tpm_merged_pathways$GO.term)

#s_230

for (term in unique_terms) {
  
 
  fmol_tpm_pathways <- subset(fmol_tpm_merged_pathways, GO.term == term)
  

  num_rows <- nrow(fmol_tpm_pathways)
  
 
  if (num_rows > 0) {
    
    
    ID <- fmol_tpm_pathways$Gene
    fmol_average_s_2h30 <- log2(rowMeans(fmol_tpm_pathways[, 34:36], na.rm = TRUE))
    tpm_average_s_2h30 <- log2(fmol_tpm_pathways[, 8])
    tpm_fmol_averge_s_2h30 <- data.frame(ID = ID, fmol_average = fmol_average_s_2h30, tpm_average = tpm_average_s_2h30)
    tpm_fmol_averge_s_2h30 <- na.omit(tpm_fmol_averge_s_2h30)
    
    # Correlation analysis
    cor_result <- cor.test(tpm_fmol_averge_s_2h30$tpm_average, tpm_fmol_averge_s_2h30$fmol_average, use = "pairwise.complete.obs", method = "pearson")
    
    # Linear regression analysis
    lm_model <- lm(fmol_average ~ tpm_average, data = tpm_fmol_averge_s_2h30)
    slope <- coef(lm_model)[2]
    slope_pvalue <- summary(lm_model)$coefficients[2, 4]
    P_value <- format(cor_result$p.value, scientific = TRUE, digits = 2)
    
  
    result_df <- rbind(result_df, data.frame(GO.term = term, Number_of_genes = num_rows, Correlation = cor_result$estimate, Slope = slope, Slope_PValue = slope_pvalue, Corr_pvalue = P_value))
  }
}


write.csv(result_df, file = "correlation_pathway_s_230.csv", row.names = FALSE)


#s_0h

unique_terms <- unique(fmol_tpm_merged_pathways$GO.term)
for (term in unique_terms) {
  

  fmol_tpm_pathways <- subset(fmol_tpm_merged_pathways, GO.term == term)
  

  num_rows <- nrow(fmol_tpm_pathways)

  if (num_rows > 0) {
    
    # Data transformation
    ID <- fmol_tpm_pathways$Gene
    fmol_average_s_0<- log2(rowMeans(fmol_tpm_pathways[, c("S1_0","S2_0","S3_0")], na.rm = TRUE))
    tpm_average_s_0 <- log2(fmol_tpm_pathways[, 3])
    tpm_fmol_averge_s_0 <- data.frame(ID = ID, fmol_average = fmol_average_s_0, tpm_average = tpm_average_s_0)
    tpm_fmol_averge_s_0 <- na.omit(tpm_fmol_averge_s_0)
    
    # Correlation analysis
    cor_result <- cor.test(tpm_fmol_averge_s_0$tpm_average, tpm_fmol_averge_s_0$fmol_average, use = "pairwise.complete.obs", method = "pearson")
    
    # Linear regression analysis
    lm_model <- lm(fmol_average ~ tpm_average, data = tpm_fmol_averge_s_0)
    slope <- coef(lm_model)[2]
    slope_pvalue <- summary(lm_model)$coefficients[2, 4]
    P_value <- format(cor_result$p.value, scientific = TRUE, digits = 2)
    
    result_df <- rbind(result_df, data.frame(GO.term = term, Number_of_genes = num_rows, Correlation = cor_result$estimate, Slope = slope, Slope_PValue = slope_pvalue, Corr_pvalue = P_value))
  }
}

write.csv(result_df, file = "correlation_pathway_s_0.csv", row.names = FALSE)



