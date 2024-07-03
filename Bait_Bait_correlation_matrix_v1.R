# Plotting Bait-Bait correlation matrices for 7013 data

## Part1: Read in SAINT file and calculate normalized SPCs for each bait (row normalized)
## Part2: Calculate pearson correlations and plot heatmap showing the correlations

#----------
## Part1: Create matrix of MSPIT_SPCs for each bait across all preys
SPC_input <- read.delim("SAINT_7013_cleaned_v2.txt", sep="\t", header=TRUE)

# Create matrix with the control counts and calculate average, add to SAINT file
library(stringr)
control_counts <- data.frame(str_split_fixed(SPC_input$ctrlCounts, "\\|", 40))
control_counts <- as.data.frame(sapply(control_counts, as.numeric))
num_controls <- ncol(control_counts)

# Create new control counts data frame with only the top 10 values in each row
control_counts$ctl.1 <- apply(control_counts[, c(1:num_controls)], 1, FUN=max, na.rm=TRUE)
control_counts$ctl.2 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[2])
control_counts$ctl.3 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[3])
control_counts$ctl.4 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[4])
control_counts$ctl.5 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[5])
control_counts$ctl.6 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[6])
control_counts$ctl.7 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[7])
control_counts$ctl.8 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[8])
control_counts$ctl.9 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[9])
control_counts$ctl.10 <- apply(control_counts[, c(1:num_controls)], 1, FUN=function(x) sort(x, TRUE)[10])

# Take the average of the top 10 controls 
control_counts_top10 <- control_counts[,c((ncol(control_counts)-9):ncol(control_counts))]
control_counts_top10$AvgSpec_top10controls <- apply(control_counts_top10[, c(1:ncol(control_counts_top10))], 1, FUN=mean, na.rm=TRUE)
control_counts_top10$AvgSpec_top10controls <- as.numeric(sprintf(control_counts_top10$AvgSpec_top10controls, fmt = '%#.1f'))

# add to SAINT file (if negative set AvgSpec at 0)
SPC_input$AvgSpec_top10controls <- control_counts_top10$AvgSpec_top10controls
SPC_input$AvgSpec_minusCtls <- SPC_input$AvgSpec - SPC_input$AvgSpec_top10controls
SPC_input$AvgSpec_minusCtls <- ifelse(SPC_input$AvgSpec_minusCtls >= 0, SPC_input$AvgSpec_minusCtls, 0)
SPC_input$Prey_name_ID <- paste0(SPC_input$PreyGene, "_", SPC_input$Prey)
  
# Expand to create matrix with SPC for each prey across all baits 
SPC_input_SPC <- SPC_input[, c(1,25,24)]
SPC_input_BFDR <- SPC_input[, c(1,25,16)]

SPC_input_sign <- SPC_input[SPC_input$BFDR <= 0.01,]
SPC_input_sign <- SPC_input_sign[, c(1,25,24)]

# Expand (make horizontal) and annotate all "NA" in SPC as 0 (no values)
library(tidyr)
SPC_input_SPC_wide <- spread(SPC_input_SPC, Bait, AvgSpec_minusCtls)
SPC_input_SPC_wide[is.na(SPC_input_SPC_wide)] <- 0
SPC_input_BFDR_wide <- spread(SPC_input_BFDR, Bait, BFDR)
SPC_input_BFDR_wide[is.na(SPC_input_BFDR_wide)] <- 1
SPC_input_sign_wide <- spread(SPC_input_sign, Bait, AvgSpec_minusCtls)
SPC_input_sign_wide[is.na(SPC_input_sign_wide)] <- 0

write.table(SPC_input_SPC_wide, "1_Boris_SAINT_7013_avgSPC_subtop10ctls_SPC_table.txt", sep="\t", row.names=FALSE)
write.table(SPC_input_BFDR_wide, "1_Boris_SAINT_7013_avgSPC_subtop10ctls_BFDR_table.txt", sep="\t", row.names=FALSE)
write.table(SPC_input_sign_wide, "1_Boris_SAINT_7013_avgSPC_subtop10ctls_SPC_table_BFDRsignif.txt", sep="\t", row.names=FALSE)

#------
## Part2: Calculate pearson correlations and plot heatmap showing the correlations
# Create matrices for input
SPC_matrix <- SPC_input_SPC_wide
row.names(SPC_matrix) <- SPC_matrix[,1]
SPC_matrix <- SPC_matrix[,-1]

SPC_sign_matrix <- SPC_input_sign_wide
row.names(SPC_sign_matrix) <- SPC_sign_matrix[,1]
SPC_sign_matrix <- SPC_sign_matrix[,-1]

# Compute correlation matrix
correlation_matrix_SPC <- cor(SPC_matrix, method = "pearson")
correlation_matrix_SPC <- round(correlation_matrix_SPC, digits=2)
write.table(correlation_matrix_SPC, "2_correlation_matrix_Bait_Bait_Pearson_SPC.txt", sep="\t", row.names=TRUE)

correlation_matrix_SPC_sign <- cor(SPC_sign_matrix, method = "pearson")
correlation_matrix_SPC_sign <- round(correlation_matrix_SPC_sign, digits=2)
write.table(correlation_matrix_SPC_sign, "2_correlation_matrix_Bait_Bait_Pearson_SPC_BFDRsignif.txt", sep="\t", row.names=TRUE)

# Create long form file with primary annotations
bait_info <- read.delim("Supplementary Table 2 - BirA bait table - v2_cleaned.csv", sep=",", header=TRUE)
bait_info_tomerge <- bait_info[,c(2,6)]
colnames(bait_info_tomerge)[1] <- "BaitName"

# Create dataframe with the pairwise correlation values and add the primary annotations
long_corr_matrix_SPC_ALL <- data.frame(row=rownames(correlation_matrix_SPC)[row(correlation_matrix_SPC)], col=colnames(correlation_matrix_SPC)[col(correlation_matrix_SPC)], corr=c(correlation_matrix_SPC))
long_corr_matrix_SPC_sign_ALL <- data.frame(row=rownames(correlation_matrix_SPC_sign)[row(correlation_matrix_SPC_sign)], col=colnames(correlation_matrix_SPC_sign)[col(correlation_matrix_SPC_sign)], corr=c(correlation_matrix_SPC_sign))

long_corr_matrix_SPC_noDup <- data.frame(row=rownames(correlation_matrix_SPC)[row(correlation_matrix_SPC)[upper.tri(correlation_matrix_SPC)]], 
                                                       col=colnames(correlation_matrix_SPC)[col(correlation_matrix_SPC)[upper.tri(correlation_matrix_SPC)]], 
                                                       corr=correlation_matrix_SPC[upper.tri(correlation_matrix_SPC)])
long_corrn_matrix_SPC_sign_noDup <- data.frame(row=rownames(correlation_matrix_SPC_sign)[row(correlation_matrix_SPC_sign)[upper.tri(correlation_matrix_SPC_sign)]], 
                                                       col=colnames(correlation_matrix_SPC_sign)[col(correlation_matrix_SPC_sign)[upper.tri(correlation_matrix_SPC_sign)]], 
                                                       corr=correlation_matrix_SPC_sign[upper.tri(correlation_matrix_SPC_sign)])

# Add in the primary annotations
long_corr_matrix_SPC_ALL_2 <- merge(long_corr_matrix_SPC_ALL,bait_info_tomerge, by.x="col", by.y="BaitName")
colnames(long_corr_matrix_SPC_ALL_2)[4] <- "Primary.annotation.col"
long_corr_matrix_SPC_ALL_2 <- merge(long_corr_matrix_SPC_ALL_2,bait_info_tomerge, by.x="row", by.y="BaitName")
colnames(long_corr_matrix_SPC_ALL_2)[5] <- "Primary.annotation.row"
long_corr_matrix_SPC_ALL_2 <- long_corr_matrix_SPC_ALL_2[,c(2,1,3,4,5)]

long_corr_matrix_SPC_sign_ALL_2 <- merge(long_corr_matrix_SPC_sign_ALL,bait_info_tomerge, by.x="col", by.y="BaitName")
colnames(long_corr_matrix_SPC_sign_ALL_2)[4] <- "Primary.annotation.col"
long_corr_matrix_SPC_sign_ALL_2 <- merge(long_corr_matrix_SPC_sign_ALL_2,bait_info_tomerge, by.x="row", by.y="BaitName")
colnames(long_corr_matrix_SPC_sign_ALL_2)[5] <- "Primary.annotation.row"
long_corr_matrix_SPC_sign_ALL_2 <- long_corr_matrix_SPC_sign_ALL_2[,c(2,1,3,4,5)]

long_corr_matrix_SPC_noDup_2 <- merge(long_corr_matrix_SPC_noDup,bait_info_tomerge, by.x="col", by.y="BaitName")
colnames(long_corr_matrix_SPC_noDup_2)[4] <- "Primary.annotation.col"
long_corr_matrix_SPC_noDup_2 <- merge(long_corr_matrix_SPC_noDup_2,bait_info_tomerge, by.x="row", by.y="BaitName")
colnames(long_corr_matrix_SPC_noDup_2)[5] <- "Primary.annotation.row"
long_corr_matrix_SPC_noDup_2 <- long_corr_matrix_SPC_noDup_2[,c(2,1,3,4,5)]

long_corrn_matrix_SPC_sign_noDup_2 <- merge(long_corrn_matrix_SPC_sign_noDup,bait_info_tomerge, by.x="col", by.y="BaitName")
colnames(long_corrn_matrix_SPC_sign_noDup_2)[4] <- "Primary.annotation.col"
long_corrn_matrix_SPC_sign_noDup_2 <- merge(long_corrn_matrix_SPC_sign_noDup_2,bait_info_tomerge, by.x="row", by.y="BaitName")
colnames(long_corrn_matrix_SPC_sign_noDup_2)[5] <- "Primary.annotation.row"
long_corrn_matrix_SPC_sign_noDup_2 <- long_corrn_matrix_SPC_sign_noDup_2[,c(2,1,3,4,5)]

write.table(long_corr_matrix_SPC_ALL_2, "3_Pairwise_correlation_matrix_Bait_Bait_Pearson_SPC_ALL.txt", sep="\t", row.names=FALSE)
write.table(long_corr_matrix_SPC_sign_ALL_2, "3_Pairwise_correlation_matrix_Bait_Bait_Pearson_SPC_BFDRsignif_ALL.txt", sep="\t", row.names=FALSE)
write.table(long_corr_matrix_SPC_noDup_2, "3_Pairwise_correlation_matrix_Bait_Bait_Pearson_SPC_noDuplicates.txt", sep="\t", row.names=FALSE)
write.table(long_corrn_matrix_SPC_sign_noDup_2, "3_Pairwise_correlation_matrix_Bait_Bait_Pearson_SPC_BFDRsignif_noDuplicates.txt", sep="\t", row.names=FALSE)

#-----------
## Part2: Plot heatmaps showing the correlations
library("gplots")
library("RColorBrewer")

# SPC 
cols <- colorRampPalette(c("blue", "white", "red"))(100)
#hmcols <- rev(redblue(100))

pdf("Rplot_Heatmap_pearson_correlation_matrix_Bait_Bait.pdf", width = 9, height = 8)
heatmap_plot <- heatmap.2(correlation_matrix_SPC, dendrogram='both', Colv = TRUE, Rowv = TRUE, scale=NULL, trace='none', col=cols, density.info="none", main = "Correlation between baits (Pearson)", cexRow = 0.2, cexCol =  0.2)
dev.off()


library("gplots")
library("RColorBrewer")

# SPC (BFDR <= 0.01)
pdf("Rplot_Heatmap_pearson_correlation_matrix_Bait_Bait_BFDRsignif.pdf", width = 9, height = 8)
heatmap_plot <- heatmap.2(correlation_matrix_SPC_sign, dendrogram='both', Colv = TRUE, Rowv = TRUE, scale=NULL, trace='none', col=cols, density.info="none", main = "Correlation between baits (Pearson):Significant (BFDR <= 0.01) only", cexRow = 0.2, cexCol =  0.2)
dev.off()
