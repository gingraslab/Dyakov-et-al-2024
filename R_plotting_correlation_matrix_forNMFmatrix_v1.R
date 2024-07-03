# Plotting correlation matrices for NMF rank basis matrix

## Part1: Read in annotated NMF basis matrix and calculate correlation between NMF ranks
## Part2: Plot heatmaps showing the correlations

#------
## Part1: Read in annotated basis matrix and count number of overlaps and create matrix
NMF_matrix <- read.delim("0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", header=TRUE)

NMF_matrix_filtered <- NMF_matrix[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38)]
row.names(NMF_matrix_filtered) <- NMF_matrix_filtered[,1]
NMF_matrix_filtered <- NMF_matrix_filtered[,-1]

# Compute correlation matrix
correlation_matrix <- cor(NMF_matrix_filtered, method = "pearson")
correlation_matrix <- round(correlation_matrix, digits=2)
write.table(correlation_matrix, "1_correlation_matrix_NMFbasis_matrix_19ranks_Boris.txt", sep="\t", row.names=TRUE)

#-----------
## Part2: Plot heatmaps showing the correlations
library("gplots")
library("RColorBrewer")

cols <- colorRampPalette(c("blue", "white", "red"))(1000)
#hmcols <- rev(redblue(100))

pdf("Rplot_Heatmap_correlation_matrix_19ranks_Boris.pdf", width = 9, height = 8)
heatmap_plot <- heatmap.2(correlation_matrix, cellnote=ifelse(correlation_matrix <= 0, NA, correlation_matrix), notecol = "black", dendrogram='none', Colv = FALSE, Rowv = FALSE, scale=NULL, trace='none', col=cols, density.info="none", main = "Correlation between ranks (Pearson)", cexRow = 1, cexCol =  1)
dev.off()

#------
#library(corrplot)

#pdf("Rplot_Heatmap_correlation_matrix_18ranks_Boris.pdf", width = 9, height = 8)
#corrplot(correlation_matrix, method="number")
#dev.off()

