library(readr)
setwd("D:/S_MMTT_analysis/DEG_analysis")
data <- read.csv("feature_counts_S_MT_wo.csv")
row.names(data) <- data$Run
colnames(data)
count_data <- data[,-1]
count_matrix <- as.matrix(count_data)
meta_data <- read.csv("sample_info_deseq2.csv")
row.names(meta_data) = meta_data$Run
meta_data
all(rownames(meta_data) == colnames(count_data))
#perform differential gene experession analysis
library(DESeq2)
meta_data$Genotype = factor(meta_data$Genotype)
meta_data$Time = factor(meta_data$Time)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_data,
                              design = ~ Genotype + Time + Genotype:Time)

#for pca 
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Time", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Time)) +
  geom_point(size=3) +  scale_shape_manual(values = c(4,8,15,16,17,18,21,22,3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + stat_ellipse()


dds <- estimateSizeFactors(dds)
c <- counts(dds,normalized=TRUE) 
counts(dds,normalized=TRUE) >= 10
keep <- rowSums(counts(dds,normalized=TRUE) >= 10) >= 2
dds <- dds[keep,]
 s <- counts(dds)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Genotype + Time)

res_LRT <- results(object = dds_lrt_time, )
library(tibble)
library(dplyr)


# Subset results for faster cluster finding
sig_res_LRT <- res_LRT%>%
  
  data.frame() %>%
  
  rownames_to_column(var="gene") %>% 
  
  as_tibble() 
  
sigLRT_genes = (sig_res_LRT[sig_res_LRT$padj <  0.001 & (abs(sig_res_LRT$log2FoldChange >0.5)|abs(sig_res_LRT$log2FoldChange < -0.5)),]) 

sigLRT_genes <- sig_res_LRT %>% 
  
  pull(gene)

setwd("D:/S_MMTT_analysis/DEG_analysis/DEGs_MMTT_LRT_spline")
write.csv(sigLRT_genes, "sig_LRT_genes_0.001_lfc0.5_nf.csv")




