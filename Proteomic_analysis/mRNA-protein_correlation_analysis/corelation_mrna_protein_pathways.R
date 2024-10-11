correlation_pathway_MT_0 <- read.csv("correlation_pathway_MT_0.csv")
correlation_pathway_MT_230 <- read.csv("correlation_pathway_MT_230.csv")
correlation_pathway_MT <- merge(correlation_pathway_MT_0[, c(1,3)],correlation_pathway_MT_230[,c(1,3)], by= "GO.term" )
rownames(correlation_pathway_MT) <- correlation_pathway_MT$GO.term
correlation_pathway_S_0 <- read.csv("correlation_pathway_s_0.csv")
correlation_pathway_S_230 <- read.csv("correlation_pathway_s_230.csv")
correlation_pathway_S <- merge(correlation_pathway_S_0[, c(1,3)],correlation_pathway_S_230[,c(1,3)], by= "GO.term" )
rownames(correlation_pathway_S) <- correlation_pathway_S$GO.term

correlation_pathway <- merge(correlation_pathway_S, correlation_pathway_MT, by = "GO.term")

correlation_pathway_1 <- correlation_pathway[c(1,3,4,8,10,33,34,35,40,46,81,88),]

rownames(correlation_pathway_1 ) <- correlation_pathway_1$GO.term
correlation_pathway_1  <- correlation_pathway_1 [-1]
correlation_pathway_1 <- as.matrix(correlation_pathway_1 )
library(RColorBrewer)

my_palette <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
library(pheatmap)

# Create the heatmap with custom options
p <- pheatmap(correlation_pathway_1, 
         display_numbers = TRUE,
         fontsize_numbers = 14,  
         cellnote_fontface = 'bold',  
         color = my_palette, legend_title= "correlation"
         )  
p

library(ggplot2)
ggsave("protein_RNA_corr.svg", p,dpi = 600, height = 5, width = 6)

slope_pathway_MT <- merge(correlation_pathway_MT_0[, c(1,4)],correlation_pathway_MT_230[,c(1,4)], by= "GO.term" )
rownames(slope_pathway_MT) <- slope_pathway_MT$GO.term
slope_pathway_S <- merge(correlation_pathway_S_0[, c(1,4)],correlation_pathway_S_230[,c(1,4)], by= "GO.term" )
rownames(slope_pathway_S) <- slope_pathway_S$GO.term
slope_pathway <- merge(slope_pathway_MT , slope_pathway_S , by = "GO.term")
slope_pathway_1 <- slope_pathway[c(1,3,4,8,10,33,34,35,40,46,81,88),]
rownames(slope_pathway_1 ) <- slope_pathway_1$GO.term
slope_pathway_1  <- slope_pathway_1 [-1]
slope_pathway_1 <- as.matrix(slope_pathway_1 )
library(RColorBrewer)

my_palette <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

p <- pheatmap(slope_pathway_1, 
              display_numbers = TRUE,
              fontsize_numbers = 15,  
              cellnote_fontface = 'bold', 
              color = my_palette,
              border_color = "black", 
              border = c(1, 1)) 
library(ggplot2)
ggsave("D:/mauscript2/figures/figures_svg/protein_RNA_slope.svg", p,dpi = 600, height = 5, width = 6)
