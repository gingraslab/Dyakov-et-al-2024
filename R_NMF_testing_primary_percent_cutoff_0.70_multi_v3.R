# Outputting Precision/Recall and Enrichment plots to assess different NMF score cutoffs (percent cutoffs)

##Part 1: Create columns with primary, secondary and tertiary ranks (use primary annotation only) reannotated with different filters
##Part 2: Run gProfiler for each of the samples using Boris' custom terms
##Part 3: Calculate precision/recall and enrichments and plot across clusters

#---------
##Part 1: Create columns with primary, secondary and tertiary ranks 
#Read in NMF file and re-annotate
NMF_matrix <- read.delim("7013_cleaned_v2_NMF_19_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix <- NMF_matrix[,-(ncol(NMF_matrix))]
Boris_curated <- read.delim("Boris_curated_lists_wIF_2.csv", sep=",",header = TRUE)

# Annotate the top, second and third ranked scores
func_second <- function(x){
  as.numeric(sort(x, TRUE)[2])
}  

func_third <- function(x){
  as.numeric(sort(x, TRUE)[3])
}  

NMF_matrix$Primary_score <- apply(NMF_matrix[, c(2:ncol(NMF_matrix))], 1, FUN=max, na.rm=TRUE)
NMF_matrix$Secondary_score <- apply(NMF_matrix[, c(2:(ncol(NMF_matrix)-1))], 1, FUN=func_second)
NMF_matrix$Tertiary_score <- apply(NMF_matrix[, c(2:(ncol(NMF_matrix)-2))], 1, FUN=func_third)

# Calculate "contiguous" score (divide secondary/tertiary by primary)
NMF_matrix$Multi_loc_score_secondary <- NMF_matrix$Secondary_score / NMF_matrix$Primary_score
NMF_matrix$Multi_loc_score_tertiary <- NMF_matrix$Tertiary_score / NMF_matrix$Primary_score

NMF_matrix[,-1] <- round(NMF_matrix[,-1], digits=4)

# Annotate column with primary, secondary and tertiary annotations
# Create function that outputs the matching column name 
maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]

library(dplyr)
NMF_matrix <- NMF_matrix %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix)-5))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix <- NMF_matrix %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix)-6))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix <- NMF_matrix %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix)-7))], 1, function(x) names(x)[maxn(3)(x)])) 

contig_matrix <- NMF_matrix[,c(1,(ncol(NMF_matrix)-7):ncol(NMF_matrix))]
write.table(NMF_matrix, "0.5_nmf_prey_matrix_input_19ranks.csv", sep=",", row.names=FALSE)

# Create a data frame showing the percentage of max value for each column and append
NMF_matrix_col <- NMF_matrix
rownames(NMF_matrix_col) <- NMF_matrix_col[,1]
NMF_matrix_col <- NMF_matrix_col[,c(2:(ncol(NMF_matrix)-8))]
NMF_matrix_col <- apply(NMF_matrix_col, 2, function(x) x/max(x))

# Re-format and add in the contiguous information
NMF_matrix_col <- as.data.frame(NMF_matrix_col)
NMF_matrix_col$PreyGene <- rownames(NMF_matrix_col)
NMF_matrix_col <- NMF_matrix_col[,c(ncol(NMF_matrix_col),1:(ncol(NMF_matrix_col)-1))]
rownames(NMF_matrix_col) <- NULL 

NMF_matrix_col <- merge(NMF_matrix_col, contig_matrix, by.x = "PreyGene", by.y = "PreyGene")

# Output percentage values for each secondary/tertiary rank
NMF_matrix_col$Primary_ratio <-as.numeric(apply(NMF_matrix_col, 1, function(x) { x[names(x)==x[names(x)=="Primary_cluster_annotation"]] }))
NMF_matrix_col$Secondary_ratio <-as.numeric(apply(NMF_matrix_col, 1, function(x) { x[names(x)==x[names(x)=="Secondary_cluster_annotation"]] }))
NMF_matrix_col$Tertiary_ratio <-as.numeric(apply(NMF_matrix_col, 1, function(x) { x[names(x)==x[names(x)=="Tertiary_cluster_annotation"]] }))

NMF_matrix_col <- NMF_matrix_col[,-c(21:28)]

# Add column names, merge with previous column and reorder
colnames(NMF_matrix_col)[2:20] <- c("X0_col_ratio","X1_col_ratio","X2_col_ratio","X3_col_ratio","X4_col_ratio","X5_col_ratio","X6_col_ratio","X7_col_ratio","X8_col_ratio","X9_col_ratio","X10_col_ratio","X11_col_ratio","X12_col_ratio","X13_col_ratio","X14_col_ratio","X15_col_ratio","X16_col_ratio","X17_col_ratio","X18_col_ratio")
NMF_matrix_col <- merge(NMF_matrix,NMF_matrix_col,by.x = "PreyGene",by.y = "PreyGene")

NMF_matrix_col <- NMF_matrix_col[,c(1,
                                    2,29,
                                    3,30,
                                    4,31,
                                    5,32,
                                    6,33,
                                    7,34,
                                    8,35,
                                    9,36,
                                    10,37,
                                    11,38,
                                    12,39,
                                    13,40,
                                    14,41,
                                    15,42,
                                    16,43,
                                    17,44,
                                    18,45,
                                    19,46,
                                    20,47,
                                    21:28,
                                    48:50)]

NMF_matrix_col[,c(2:44,48:50)] <- round(NMF_matrix_col[,c(2:44,48:50)], digits=4)

write.table(NMF_matrix_col, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks.csv", sep=",", row.names=FALSE)

#-----------
# Create NMF score matrices filtered by different percent of max in column filters
# Matrix_0_col_ratio
NMF_matrix_col_0_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0_multi_0.70$Primary_ratio > 0, NMF_matrix_col_0_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0_multi_0.70$Primary_ratio > 0 & NMF_matrix_col_0_multi_0.70$Secondary_ratio > 0 & NMF_matrix_col_0_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0_multi_0.70$Primary_ratio > 0 & NMF_matrix_col_0_multi_0.70$Tertiary_ratio > 0 & NMF_matrix_col_0_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.05_col_ratio
NMF_matrix_col_0.05_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05, NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Secondary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Tertiary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.05_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.1_col_ratio
NMF_matrix_col_0.10_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.10_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.10_multi_0.70$Primary_ratio >= 0.10, NMF_matrix_col_0.10_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.10_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.10_multi_0.70$Primary_ratio >= 0.10 & NMF_matrix_col_0.10_multi_0.70$Secondary_ratio >= 0.10 & NMF_matrix_col_0.10_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.10_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.10_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.10_multi_0.70$Primary_ratio >= 0.10 & NMF_matrix_col_0.10_multi_0.70$Tertiary_ratio >= 0.10 & NMF_matrix_col_0.10_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.10_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.10_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.10_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.10_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.10_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.10_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.15_col_ratio
NMF_matrix_col_0.15_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.15_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.15_multi_0.70$Primary_ratio >= 0.15, NMF_matrix_col_0.15_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.15_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.15_multi_0.70$Primary_ratio >= 0.15 & NMF_matrix_col_0.15_multi_0.70$Secondary_ratio >= 0.15 & NMF_matrix_col_0.15_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.15_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.15_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.15_multi_0.70$Primary_ratio >= 0.15 & NMF_matrix_col_0.15_multi_0.70$Tertiary_ratio >= 0.15 & NMF_matrix_col_0.15_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.15_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.15_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.15_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.15_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.15_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.15_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.20_col_ratio
NMF_matrix_col_0.20_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.20_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.20_multi_0.70$Primary_ratio >= 0.20, NMF_matrix_col_0.20_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.20_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.20_multi_0.70$Primary_ratio >= 0.20 & NMF_matrix_col_0.20_multi_0.70$Secondary_ratio >= 0.20 & NMF_matrix_col_0.20_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.20_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.20_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.20_multi_0.70$Primary_ratio >= 0.20 & NMF_matrix_col_0.20_multi_0.70$Tertiary_ratio >= 0.20 & NMF_matrix_col_0.20_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.20_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.20_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.20_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.20_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.20_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.20_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.25_col_ratio
NMF_matrix_col_0.25_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.25_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.25_multi_0.70$Primary_ratio >= 0.25, NMF_matrix_col_0.25_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.25_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.25_multi_0.70$Primary_ratio >= 0.25 & NMF_matrix_col_0.25_multi_0.70$Secondary_ratio >= 0.25 & NMF_matrix_col_0.25_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.25_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.25_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.25_multi_0.70$Primary_ratio >= 0.25 & NMF_matrix_col_0.25_multi_0.70$Tertiary_ratio >= 0.25 & NMF_matrix_col_0.25_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.25_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.25_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.25_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.25_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.25_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.25_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.30_col_ratio
NMF_matrix_col_0.30_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.30_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.30_multi_0.70$Primary_ratio >= 0.30, NMF_matrix_col_0.30_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.30_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.30_multi_0.70$Primary_ratio >= 0.30 & NMF_matrix_col_0.30_multi_0.70$Secondary_ratio >= 0.30 & NMF_matrix_col_0.30_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.30_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.30_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.30_multi_0.70$Primary_ratio >= 0.30 & NMF_matrix_col_0.30_multi_0.70$Tertiary_ratio >= 0.30 & NMF_matrix_col_0.30_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.30_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.30_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.30_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.30_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.30_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.30_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.35_col_ratio
NMF_matrix_col_0.35_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.35_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.35_multi_0.70$Primary_ratio >= 0.35, NMF_matrix_col_0.35_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.35_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.35_multi_0.70$Primary_ratio >= 0.35 & NMF_matrix_col_0.35_multi_0.70$Secondary_ratio >= 0.35 & NMF_matrix_col_0.35_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.35_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.35_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.35_multi_0.70$Primary_ratio >= 0.35 & NMF_matrix_col_0.35_multi_0.70$Tertiary_ratio >= 0.35 & NMF_matrix_col_0.35_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.35_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.35_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.35_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.35_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.35_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.35_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.40_col_ratio
NMF_matrix_col_0.40_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.40_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.40_multi_0.70$Primary_ratio >= 0.40, NMF_matrix_col_0.40_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.40_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.40_multi_0.70$Primary_ratio >= 0.40 & NMF_matrix_col_0.40_multi_0.70$Secondary_ratio >= 0.40 & NMF_matrix_col_0.40_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.40_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.40_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.40_multi_0.70$Primary_ratio >= 0.40 & NMF_matrix_col_0.40_multi_0.70$Tertiary_ratio >= 0.40 & NMF_matrix_col_0.40_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.40_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.40_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.40_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.40_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.40_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.40_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.45_col_ratio
NMF_matrix_col_0.45_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.45_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.45_multi_0.70$Primary_ratio >= 0.45, NMF_matrix_col_0.45_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.45_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.45_multi_0.70$Primary_ratio >= 0.45 & NMF_matrix_col_0.45_multi_0.70$Secondary_ratio >= 0.45 & NMF_matrix_col_0.45_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.45_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.45_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.45_multi_0.70$Primary_ratio >= 0.45 & NMF_matrix_col_0.45_multi_0.70$Tertiary_ratio >= 0.45 & NMF_matrix_col_0.45_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.45_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.45_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.45_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.45_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.45_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.45_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.50_col_ratio
NMF_matrix_col_0.50_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.50_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.50_multi_0.70$Primary_ratio >= 0.50, NMF_matrix_col_0.50_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.50_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.50_multi_0.70$Primary_ratio >= 0.50 & NMF_matrix_col_0.50_multi_0.70$Secondary_ratio >= 0.50 & NMF_matrix_col_0.50_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.50_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.50_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.50_multi_0.70$Primary_ratio >= 0.50 & NMF_matrix_col_0.50_multi_0.70$Tertiary_ratio >= 0.50 & NMF_matrix_col_0.50_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.50_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.50_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.50_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.50_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.50_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.50_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)


##Part 2: Run gProfiler for each of the samples using Boris' custom terms
library("gprofiler2")
library(openxlsx)

curated_terms <- c("Nuclear speckle curated",
                   "Paraspeckle curated", 
                   "Comprehensive nucleolus",
                   "Cajal body curated",
                   "PML body curated"
)

# Run gProfiler on each of the matrices (without sign filter) to get stats + enrichments to plot

# Filtering for only top and close secondary ranks for curated terms (from previous analysis)
ranks <- c("X0","X4","X10","X12","X14")
#ranks <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16")
#ranks <- c("X0","X1","X2")

# 0 Filter
df_0 <- data.frame()
wb_0 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0_multi_0.70[which(NMF_matrix_col_0_multi_0.70$Primary_cluster_annotation == i | 
                                               NMF_matrix_col_0_multi_0.70$Secondary_cluster_annotation == i |
                                               NMF_matrix_col_0_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0, sheet_name)
  writeData(wb_0, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0 <- rbind(df2, df_0)
}

saveWorkbook(wb_0, 'gprofiler_Results_NMF19_customToken_top10ctls_0_percentfilter.xlsx', overwrite = TRUE)
df_0$parents <- as.character(df_0$parents)
df_0$Primary_cutoff <- 0
df_0 <- df_0[,c(17,1:16)]
df_0_curated <- df_0[df_0$term_name %in% curated_terms,]
df_0_curated$F1score <- (2*(df_0_curated$precision*df_0_curated$recall)) / (df_0_curated$precision + df_0_curated$recall)
df_0_curated$F0.5score <- ((1 + 0.5^2) * df_0_curated$precision * df_0_curated$recall) / (0.5^2 * df_0_curated$precision + df_0_curated$recall)
df_0_curated$F2score <- ((1 + 2^2) * df_0_curated$precision * df_0_curated$recall) / (2^2 * df_0_curated$precision + df_0_curated$recall)

write.table(df_0, "1_gProfile_Results_NMF19_customToken_top10ctls_0_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.05 Filter
df_0.05 <- data.frame()
wb_0.05 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.05_multi_0.70[which(NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation == i | 
                                                  NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation == i |
                                                  NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.05, sheet_name)
  writeData(wb_0.05, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.05 <- rbind(df2, df_0.05)
}

saveWorkbook(wb_0.05, 'gprofiler_Results_NMF19_customToken_top10ctls_0.05_percentfilter.xlsx', overwrite = TRUE)
df_0.05$parents <- as.character(df_0.05$parents)
df_0.05$Primary_cutoff <- 0.05
df_0.05 <- df_0.05[,c(17,1:16)]
df_0.05_curated <- df_0.05[df_0.05$term_name %in% curated_terms,]
df_0.05_curated$F1score <- (2*(df_0.05_curated$precision*df_0.05_curated$recall)) / (df_0.05_curated$precision + df_0.05_curated$recall)
df_0.05_curated$F0.5score <- ((1 + 0.5^2) * df_0.05_curated$precision * df_0.05_curated$recall) / (0.5^2 * df_0.05_curated$precision + df_0.05_curated$recall)
df_0.05_curated$F2score <- ((1 + 2^2) * df_0.05_curated$precision * df_0.05_curated$recall) / (2^2 * df_0.05_curated$precision + df_0.05_curated$recall)

write.table(df_0.05, "1_gProfile_Results_NMF19_customToken_top10ctls_0.05_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.10 Filter
df_0.10 <- data.frame()
wb_0.10 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.10_multi_0.70[which(NMF_matrix_col_0.10_multi_0.70$Primary_cluster_annotation == i | 
                                                  NMF_matrix_col_0.10_multi_0.70$Secondary_cluster_annotation == i |
                                                  NMF_matrix_col_0.10_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.10, sheet_name)
  writeData(wb_0.10, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.10 <- rbind(df2, df_0.10)
}

saveWorkbook(wb_0.10, 'gprofiler_Results_NMF19_customToken_top10ctls_0.10_percentfilter.xlsx', overwrite = TRUE)
df_0.10$parents <- as.character(df_0.10$parents)
df_0.10$Primary_cutoff <- 0.10
df_0.10 <- df_0.10[,c(17,1:16)]
df_0.10_curated <- df_0.10[df_0.10$term_name %in% curated_terms,]
df_0.10_curated$F1score <- (2*(df_0.10_curated$precision*df_0.10_curated$recall)) / (df_0.10_curated$precision + df_0.10_curated$recall)
df_0.10_curated$F0.5score <- ((1 + 0.5^2) * df_0.10_curated$precision * df_0.10_curated$recall) / (0.5^2 * df_0.10_curated$precision + df_0.10_curated$recall)
df_0.10_curated$F2score <- ((1 + 2^2) * df_0.10_curated$precision * df_0.10_curated$recall) / (2^2 * df_0.10_curated$precision + df_0.10_curated$recall)

write.table(df_0.10, "1_gProfile_Results_NMF19_customToken_top10ctls_0.10_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.15 Filter
df_0.15 <- data.frame()
wb_0.15 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.15_multi_0.70[which(NMF_matrix_col_0.15_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.15_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.15_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.15, sheet_name)
  writeData(wb_0.15, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.15 <- rbind(df2, df_0.15)
}

saveWorkbook(wb_0.15, 'gprofiler_Results_NMF19_customToken_top10ctls_0.15_percentfilter.xlsx', overwrite = TRUE)
df_0.15$parents <- as.character(df_0.15$parents)
df_0.15$Primary_cutoff <- 0.15
df_0.15 <- df_0.15[,c(17,1:16)]
df_0.15_curated <- df_0.15[df_0.15$term_name %in% curated_terms,]
df_0.15_curated$F1score <- (2*(df_0.15_curated$precision*df_0.15_curated$recall)) / (df_0.15_curated$precision + df_0.15_curated$recall)
df_0.15_curated$F0.5score <- ((1 + 0.5^2) * df_0.15_curated$precision * df_0.15_curated$recall) / (0.5^2 * df_0.15_curated$precision + df_0.15_curated$recall)
df_0.15_curated$F2score <- ((1 + 2^2) * df_0.15_curated$precision * df_0.15_curated$recall) / (2^2 * df_0.15_curated$precision + df_0.15_curated$recall)

write.table(df_0.15, "1_gProfile_Results_NMF19_customToken_top10ctls_0.15_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.20 Filter
df_0.20 <- data.frame()
wb_0.20 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.20_multi_0.70[which(NMF_matrix_col_0.20_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.20_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.20_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.20, sheet_name)
  writeData(wb_0.20, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.20 <- rbind(df2, df_0.20)
}

saveWorkbook(wb_0.20, 'gprofiler_Results_NMF19_customToken_top10ctls_0.20_percentfilter.xlsx', overwrite = TRUE)
df_0.20$parents <- as.character(df_0.20$parents)
df_0.20$Primary_cutoff <- 0.20
df_0.20 <- df_0.20[,c(17,1:16)]
df_0.20_curated <- df_0.20[df_0.20$term_name %in% curated_terms,]
df_0.20_curated$F1score <- (2*(df_0.20_curated$precision*df_0.20_curated$recall)) / (df_0.20_curated$precision + df_0.20_curated$recall)
df_0.20_curated$F0.5score <- ((1 + 0.5^2) * df_0.20_curated$precision * df_0.20_curated$recall) / (0.5^2 * df_0.20_curated$precision + df_0.20_curated$recall)
df_0.20_curated$F2score <- ((1 + 2^2) * df_0.20_curated$precision * df_0.20_curated$recall) / (2^2 * df_0.20_curated$precision + df_0.20_curated$recall)

write.table(df_0.20, "1_gProfile_Results_NMF19_customToken_top10ctls_0.20_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.25 Filter
df_0.25 <- data.frame()
wb_0.25 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.25_multi_0.70[which(NMF_matrix_col_0.25_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.25_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.25_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.25, sheet_name)
  writeData(wb_0.25, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.25 <- rbind(df2, df_0.25)
}

saveWorkbook(wb_0.25, 'gprofiler_Results_NMF19_customToken_top10ctls_0.25_percentfilter.xlsx', overwrite = TRUE)
df_0.25$parents <- as.character(df_0.25$parents)
df_0.25$Primary_cutoff <- 0.25
df_0.25 <- df_0.25[,c(17,1:16)]
df_0.25_curated <- df_0.25[df_0.25$term_name %in% curated_terms,]
df_0.25_curated$F1score <- (2*(df_0.25_curated$precision*df_0.25_curated$recall)) / (df_0.25_curated$precision + df_0.25_curated$recall)
df_0.25_curated$F0.5score <- ((1 + 0.5^2) * df_0.25_curated$precision * df_0.25_curated$recall) / (0.5^2 * df_0.25_curated$precision + df_0.25_curated$recall)
df_0.25_curated$F2score <- ((1 + 2^2) * df_0.25_curated$precision * df_0.25_curated$recall) / (2^2 * df_0.25_curated$precision + df_0.25_curated$recall)

write.table(df_0.25, "1_gProfile_Results_NMF19_customToken_top10ctls_0.25_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.30 Filter
df_0.30 <- data.frame()
wb_0.30 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.30_multi_0.70[which(NMF_matrix_col_0.30_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.30_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.30_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.30, sheet_name)
  writeData(wb_0.30, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.30 <- rbind(df2, df_0.30)
}

saveWorkbook(wb_0.30, 'gprofiler_Results_NMF19_customToken_top10ctls_0.30_percentfilter.xlsx', overwrite = TRUE)
df_0.30$parents <- as.character(df_0.30$parents)
df_0.30$Primary_cutoff <- 0.30
df_0.30 <- df_0.30[,c(17,1:16)]
df_0.30_curated <- df_0.30[df_0.30$term_name %in% curated_terms,]
df_0.30_curated$F1score <- (2*(df_0.30_curated$precision*df_0.30_curated$recall)) / (df_0.30_curated$precision + df_0.30_curated$recall)
df_0.30_curated$F0.5score <- ((1 + 0.5^2) * df_0.30_curated$precision * df_0.30_curated$recall) / (0.5^2 * df_0.30_curated$precision + df_0.30_curated$recall)
df_0.30_curated$F2score <- ((1 + 2^2) * df_0.30_curated$precision * df_0.30_curated$recall) / (2^2 * df_0.30_curated$precision + df_0.30_curated$recall)

write.table(df_0.30, "1_gProfile_Results_NMF19_customToken_top10ctls_0.30_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.35 Filter
df_0.35 <- data.frame()
wb_0.35 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.35_multi_0.70[which(NMF_matrix_col_0.35_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.35_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.35_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.35, sheet_name)
  writeData(wb_0.35, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.35 <- rbind(df2, df_0.35)
}

saveWorkbook(wb_0.35, 'gprofiler_Results_NMF19_customToken_top10ctls_0.35_percentfilter.xlsx', overwrite = TRUE)
df_0.35$parents <- as.character(df_0.35$parents)
df_0.35$Primary_cutoff <- 0.35
df_0.35 <- df_0.35[,c(17,1:16)]
df_0.35_curated <- df_0.35[df_0.35$term_name %in% curated_terms,]
df_0.35_curated$F1score <- (2*(df_0.35_curated$precision*df_0.35_curated$recall)) / (df_0.35_curated$precision + df_0.35_curated$recall)
df_0.35_curated$F0.5score <- ((1 + 0.5^2) * df_0.35_curated$precision * df_0.35_curated$recall) / (0.5^2 * df_0.35_curated$precision + df_0.35_curated$recall)
df_0.35_curated$F2score <- ((1 + 2^2) * df_0.35_curated$precision * df_0.35_curated$recall) / (2^2 * df_0.35_curated$precision + df_0.35_curated$recall)

write.table(df_0.35, "1_gProfile_Results_NMF19_customToken_top10ctls_0.35_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.40 Filter
df_0.40 <- data.frame()
wb_0.40 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.40_multi_0.70[which(NMF_matrix_col_0.40_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.40_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.40_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.40, sheet_name)
  writeData(wb_0.40, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.40 <- rbind(df2, df_0.40)
}

saveWorkbook(wb_0.40, 'gprofiler_Results_NMF19_customToken_top10ctls_0.40_percentfilter.xlsx', overwrite = TRUE)
df_0.40$parents <- as.character(df_0.40$parents)
df_0.40$Primary_cutoff <- 0.40
df_0.40 <- df_0.40[,c(17,1:16)]
df_0.40_curated <- df_0.40[df_0.40$term_name %in% curated_terms,]
df_0.40_curated$F1score <- (2*(df_0.40_curated$precision*df_0.40_curated$recall)) / (df_0.40_curated$precision + df_0.40_curated$recall)
df_0.40_curated$F0.5score <- ((1 + 0.5^2) * df_0.40_curated$precision * df_0.40_curated$recall) / (0.5^2 * df_0.40_curated$precision + df_0.40_curated$recall)
df_0.40_curated$F2score <- ((1 + 2^2) * df_0.40_curated$precision * df_0.40_curated$recall) / (2^2 * df_0.40_curated$precision + df_0.40_curated$recall)

write.table(df_0.40, "1_gProfile_Results_NMF19_customToken_top10ctls_0.40_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.45 Filter
df_0.45 <- data.frame()
wb_0.45 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.45_multi_0.70[which(NMF_matrix_col_0.45_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.45_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.45_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.45, sheet_name)
  writeData(wb_0.45, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.45 <- rbind(df2, df_0.45)
}

saveWorkbook(wb_0.45, 'gprofiler_Results_NMF19_customToken_top10ctls_0.45_percentfilter.xlsx', overwrite = TRUE)
df_0.45$parents <- as.character(df_0.45$parents)
df_0.45$Primary_cutoff <- 0.45
df_0.45 <- df_0.45[,c(17,1:16)]
df_0.45_curated <- df_0.45[df_0.45$term_name %in% curated_terms,]
df_0.45_curated$F1score <- (2*(df_0.45_curated$precision*df_0.45_curated$recall)) / (df_0.45_curated$precision + df_0.45_curated$recall)
df_0.45_curated$F0.5score <- ((1 + 0.5^2) * df_0.45_curated$precision * df_0.45_curated$recall) / (0.5^2 * df_0.45_curated$precision + df_0.45_curated$recall)
df_0.45_curated$F2score <- ((1 + 2^2) * df_0.45_curated$precision * df_0.45_curated$recall) / (2^2 * df_0.45_curated$precision + df_0.45_curated$recall)

write.table(df_0.45, "1_gProfile_Results_NMF19_customToken_top10ctls_0.45_percentfilter.txt", sep="\t", row.names=FALSE)


# 0.50 Filter
df_0.50 <- data.frame()
wb_0.50 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.50_multi_0.70[which(NMF_matrix_col_0.50_multi_0.70$Primary_cluster_annotation == i | 
                                              NMF_matrix_col_0.50_multi_0.70$Secondary_cluster_annotation == i |
                                              NMF_matrix_col_0.50_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_0.50, sheet_name)
  writeData(wb_0.50, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.50 <- rbind(df2, df_0.50)
}

saveWorkbook(wb_0.50, 'gprofiler_Results_NMF19_customToken_top10ctls_0.50_percentfilter.xlsx', overwrite = TRUE)
df_0.50$parents <- as.character(df_0.50$parents)
df_0.50$Primary_cutoff <- 0.50
df_0.50 <- df_0.50[,c(17,1:16)]
df_0.50_curated <- df_0.50[df_0.50$term_name %in% curated_terms,]
df_0.50_curated$F1score <- (2*(df_0.50_curated$precision*df_0.50_curated$recall)) / (df_0.50_curated$precision + df_0.50_curated$recall)
df_0.50_curated$F0.5score <- ((1 + 0.5^2) * df_0.50_curated$precision * df_0.50_curated$recall) / (0.5^2 * df_0.50_curated$precision + df_0.50_curated$recall)
df_0.50_curated$F2score <- ((1 + 2^2) * df_0.50_curated$precision * df_0.50_curated$recall) / (2^2 * df_0.50_curated$precision + df_0.50_curated$recall)

write.table(df_0.50, "1_gProfile_Results_NMF19_customToken_top10ctls_0.50_percentfilter.txt", sep="\t", row.names=FALSE)

# Combine all the gProfiler curated results into one
df_combined_curated <- rbind(df_0_curated,
                             df_0.05_curated,
                             df_0.10_curated,
                             df_0.15_curated,
                             df_0.20_curated,
                             df_0.25_curated,
                             df_0.30_curated,
                             df_0.35_curated,
                             df_0.40_curated,
                             df_0.45_curated,
                             df_0.50_curated
)

write.table(df_combined_curated, "2_gProfile_Results_NMF19_customToken_top10ctls_0.70_multi_curated_ALL_percentfilter.txt", sep="\t", row.names=FALSE)

#--------
##Part 3: Calculate precision/recall and enrichments and plot across clusters
# PS stats (for top cluster)
# Create a table with the total number and 
PS_0.00 <- df_0_curated[df_0_curated$term_name == "Paraspeckle curated",]
PS_0.00$negative_log10_pvalue <- -log10(PS_0.00$adjusted_p_value)
PS_0.00$Number_clusters_with_term <- nrow(PS_0.00[PS_0.00$significant == TRUE,])

PS_0.05 <- df_0.05_curated[df_0.05_curated$term_name == "Paraspeckle curated",]
PS_0.05$negative_log10_pvalue <- -log10(PS_0.05$adjusted_p_value)
PS_0.05$Number_clusters_with_term <- nrow(PS_0.05[PS_0.05$significant == TRUE,])

PS_0.10 <- df_0.10_curated[df_0.10_curated$term_name == "Paraspeckle curated",]
PS_0.10$negative_log10_pvalue <- -log10(PS_0.10$adjusted_p_value)
PS_0.10$Number_clusters_with_term <- nrow(PS_0.10[PS_0.10$significant == TRUE,])

PS_0.15 <- df_0.15_curated[df_0.15_curated$term_name == "Paraspeckle curated",]
PS_0.15$negative_log10_pvalue <- -log10(PS_0.15$adjusted_p_value)
PS_0.15$Number_clusters_with_term <- nrow(PS_0.15[PS_0.15$significant == TRUE,])

PS_0.20 <- df_0.20_curated[df_0.20_curated$term_name == "Paraspeckle curated",]
PS_0.20$negative_log10_pvalue <- -log10(PS_0.20$adjusted_p_value)
PS_0.20$Number_clusters_with_term <- nrow(PS_0.20[PS_0.20$significant == TRUE,])

PS_0.25 <- df_0.25_curated[df_0.25_curated$term_name == "Paraspeckle curated",]
PS_0.25$negative_log10_pvalue <- -log10(PS_0.25$adjusted_p_value)
PS_0.25$Number_clusters_with_term <- nrow(PS_0.25[PS_0.25$significant == TRUE,])

PS_0.30 <- df_0.30_curated[df_0.30_curated$term_name == "Paraspeckle curated",]
PS_0.30$negative_log10_pvalue <- -log10(PS_0.30$adjusted_p_value)
PS_0.30$Number_clusters_with_term <- nrow(PS_0.30[PS_0.30$significant == TRUE,])

PS_0.35 <- df_0.35_curated[df_0.35_curated$term_name == "Paraspeckle curated",]
PS_0.35$negative_log10_pvalue <- -log10(PS_0.35$adjusted_p_value)
PS_0.35$Number_clusters_with_term <- nrow(PS_0.35[PS_0.35$significant == TRUE,])

PS_0.40 <- df_0.40_curated[df_0.40_curated$term_name == "Paraspeckle curated",]
PS_0.40$negative_log10_pvalue <- -log10(PS_0.40$adjusted_p_value)
PS_0.40$Number_clusters_with_term <- nrow(PS_0.40[PS_0.40$significant == TRUE,])

PS_0.45 <- df_0.45_curated[df_0.45_curated$term_name == "Paraspeckle curated",]
PS_0.45$negative_log10_pvalue <- -log10(PS_0.45$adjusted_p_value)
PS_0.45$Number_clusters_with_term <- nrow(PS_0.45[PS_0.45$significant == TRUE,])

PS_0.50 <- df_0.50_curated[df_0.50_curated$term_name == "Paraspeckle curated",]
PS_0.50$negative_log10_pvalue <- -log10(PS_0.50$adjusted_p_value)
PS_0.50$Number_clusters_with_term <- nrow(PS_0.50[PS_0.50$significant == TRUE,])

top_rank_PS <- rbind(PS_0.00[which.max(PS_0.00$negative_log10_pvalue),],
                     PS_0.05[which.max(PS_0.05$negative_log10_pvalue),],
                     PS_0.10[which.max(PS_0.10$negative_log10_pvalue),],
                     PS_0.15[which.max(PS_0.15$negative_log10_pvalue),],
                     PS_0.20[which.max(PS_0.20$negative_log10_pvalue),],
                     PS_0.25[which.max(PS_0.25$negative_log10_pvalue),],
                     PS_0.30[which.max(PS_0.30$negative_log10_pvalue),],
                     PS_0.35[which.max(PS_0.35$negative_log10_pvalue),],
                     PS_0.40[which.max(PS_0.40$negative_log10_pvalue),],
                     PS_0.45[which.max(PS_0.45$negative_log10_pvalue),],
                     PS_0.50[which.max(PS_0.50$negative_log10_pvalue),]
)

# NS stats (for top cluster)
# Create a table with the total number and 
NS_0.00 <- df_0_curated[df_0_curated$term_name == "Nuclear speckle curated",]
NS_0.00$negative_log10_pvalue <- -log10(NS_0.00$adjusted_p_value)
NS_0.00$Number_clusters_with_term <- nrow(NS_0.00[NS_0.00$significant == TRUE,])

NS_0.05 <- df_0.05_curated[df_0.05_curated$term_name == "Nuclear speckle curated",]
NS_0.05$negative_log10_pvalue <- -log10(NS_0.05$adjusted_p_value)
NS_0.05$Number_clusters_with_term <- nrow(NS_0.05[NS_0.05$significant == TRUE,])

NS_0.10 <- df_0.10_curated[df_0.10_curated$term_name == "Nuclear speckle curated",]
NS_0.10$negative_log10_pvalue <- -log10(NS_0.10$adjusted_p_value)
NS_0.10$Number_clusters_with_term <- nrow(NS_0.10[NS_0.10$significant == TRUE,])

NS_0.15 <- df_0.15_curated[df_0.15_curated$term_name == "Nuclear speckle curated",]
NS_0.15$negative_log10_pvalue <- -log10(NS_0.15$adjusted_p_value)
NS_0.15$Number_clusters_with_term <- nrow(NS_0.15[NS_0.15$significant == TRUE,])

NS_0.20 <- df_0.20_curated[df_0.20_curated$term_name == "Nuclear speckle curated",]
NS_0.20$negative_log10_pvalue <- -log10(NS_0.20$adjusted_p_value)
NS_0.20$Number_clusters_with_term <- nrow(NS_0.20[NS_0.20$significant == TRUE,])

NS_0.25 <- df_0.25_curated[df_0.25_curated$term_name == "Nuclear speckle curated",]
NS_0.25$negative_log10_pvalue <- -log10(NS_0.25$adjusted_p_value)
NS_0.25$Number_clusters_with_term <- nrow(NS_0.25[NS_0.25$significant == TRUE,])

NS_0.30 <- df_0.30_curated[df_0.30_curated$term_name == "Nuclear speckle curated",]
NS_0.30$negative_log10_pvalue <- -log10(NS_0.30$adjusted_p_value)
NS_0.30$Number_clusters_with_term <- nrow(NS_0.30[NS_0.30$significant == TRUE,])

NS_0.35 <- df_0.35_curated[df_0.35_curated$term_name == "Nuclear speckle curated",]
NS_0.35$negative_log10_pvalue <- -log10(NS_0.35$adjusted_p_value)
NS_0.35$Number_clusters_with_term <- nrow(NS_0.35[NS_0.35$significant == TRUE,])

NS_0.40 <- df_0.40_curated[df_0.40_curated$term_name == "Nuclear speckle curated",]
NS_0.40$negative_log10_pvalue <- -log10(NS_0.40$adjusted_p_value)
NS_0.40$Number_clusters_with_term <- nrow(NS_0.40[NS_0.40$significant == TRUE,])

NS_0.45 <- df_0.45_curated[df_0.45_curated$term_name == "Nuclear speckle curated",]
NS_0.45$negative_log10_pvalue <- -log10(NS_0.45$adjusted_p_value)
NS_0.45$Number_clusters_with_term <- nrow(NS_0.45[NS_0.45$significant == TRUE,])

NS_0.50 <- df_0.50_curated[df_0.50_curated$term_name == "Nuclear speckle curated",]
NS_0.50$negative_log10_pvalue <- -log10(NS_0.50$adjusted_p_value)
NS_0.50$Number_clusters_with_term <- nrow(NS_0.50[NS_0.50$significant == TRUE,])

top_rank_NS <- rbind(NS_0.00[which.max(NS_0.00$negative_log10_pvalue),],
                     NS_0.05[which.max(NS_0.05$negative_log10_pvalue),],
                     NS_0.10[which.max(NS_0.10$negative_log10_pvalue),],
                     NS_0.15[which.max(NS_0.15$negative_log10_pvalue),],
                     NS_0.20[which.max(NS_0.20$negative_log10_pvalue),],
                     NS_0.25[which.max(NS_0.25$negative_log10_pvalue),],
                     NS_0.30[which.max(NS_0.30$negative_log10_pvalue),],
                     NS_0.35[which.max(NS_0.35$negative_log10_pvalue),],
                     NS_0.40[which.max(NS_0.40$negative_log10_pvalue),],
                     NS_0.45[which.max(NS_0.45$negative_log10_pvalue),],
                     NS_0.50[which.max(NS_0.50$negative_log10_pvalue),]
)

# Nucleolus stats (for top cluster)
# Create a table with the total number and 
NUC_0.00 <- df_0_curated[df_0_curated$term_name == "Comprehensive nucleolus",]
NUC_0.00$negative_log10_pvalue <- -log10(NUC_0.00$adjusted_p_value)
NUC_0.00$Number_clusters_with_term <- nrow(NUC_0.00[NUC_0.00$significant == TRUE,])

NUC_0.05 <- df_0.05_curated[df_0.05_curated$term_name == "Comprehensive nucleolus",]
NUC_0.05$negative_log10_pvalue <- -log10(NUC_0.05$adjusted_p_value)
NUC_0.05$Number_clusters_with_term <- nrow(NUC_0.05[NUC_0.05$significant == TRUE,])

NUC_0.10 <- df_0.10_curated[df_0.10_curated$term_name == "Comprehensive nucleolus",]
NUC_0.10$negative_log10_pvalue <- -log10(NUC_0.10$adjusted_p_value)
NUC_0.10$Number_clusters_with_term <- nrow(NUC_0.10[NUC_0.10$significant == TRUE,])

NUC_0.15 <- df_0.15_curated[df_0.15_curated$term_name == "Comprehensive nucleolus",]
NUC_0.15$negative_log10_pvalue <- -log10(NUC_0.15$adjusted_p_value)
NUC_0.15$Number_clusters_with_term <- nrow(NUC_0.15[NUC_0.15$significant == TRUE,])

NUC_0.20 <- df_0.20_curated[df_0.20_curated$term_name == "Comprehensive nucleolus",]
NUC_0.20$negative_log10_pvalue <- -log10(NUC_0.20$adjusted_p_value)
NUC_0.20$Number_clusters_with_term <- nrow(NUC_0.20[NUC_0.20$significant == TRUE,])

NUC_0.25 <- df_0.25_curated[df_0.25_curated$term_name == "Comprehensive nucleolus",]
NUC_0.25$negative_log10_pvalue <- -log10(NUC_0.25$adjusted_p_value)
NUC_0.25$Number_clusters_with_term <- nrow(NUC_0.25[NUC_0.25$significant == TRUE,])

NUC_0.30 <- df_0.30_curated[df_0.30_curated$term_name == "Comprehensive nucleolus",]
NUC_0.30$negative_log10_pvalue <- -log10(NUC_0.30$adjusted_p_value)
NUC_0.30$Number_clusters_with_term <- nrow(NUC_0.30[NUC_0.30$significant == TRUE,])

NUC_0.35 <- df_0.35_curated[df_0.35_curated$term_name == "Comprehensive nucleolus",]
NUC_0.35$negative_log10_pvalue <- -log10(NUC_0.35$adjusted_p_value)
NUC_0.35$Number_clusters_with_term <- nrow(NUC_0.35[NUC_0.35$significant == TRUE,])

NUC_0.40 <- df_0.40_curated[df_0.40_curated$term_name == "Comprehensive nucleolus",]
NUC_0.40$negative_log10_pvalue <- -log10(NUC_0.40$adjusted_p_value)
NUC_0.40$Number_clusters_with_term <- nrow(NUC_0.40[NUC_0.40$significant == TRUE,])

NUC_0.45 <- df_0.45_curated[df_0.45_curated$term_name == "Comprehensive nucleolus",]
NUC_0.45$negative_log10_pvalue <- -log10(NUC_0.45$adjusted_p_value)
NUC_0.45$Number_clusters_with_term <- nrow(NUC_0.45[NUC_0.45$significant == TRUE,])

NUC_0.50 <- df_0.50_curated[df_0.50_curated$term_name == "Comprehensive nucleolus",]
NUC_0.50$negative_log10_pvalue <- -log10(NUC_0.50$adjusted_p_value)
NUC_0.50$Number_clusters_with_term <- nrow(NUC_0.50[NUC_0.50$significant == TRUE,])

top_rank_NUC <- rbind(NUC_0.00[which.max(NUC_0.00$negative_log10_pvalue),],
                      NUC_0.05[which.max(NUC_0.05$negative_log10_pvalue),],
                      NUC_0.10[which.max(NUC_0.10$negative_log10_pvalue),],
                      NUC_0.15[which.max(NUC_0.15$negative_log10_pvalue),],
                      NUC_0.20[which.max(NUC_0.20$negative_log10_pvalue),],
                      NUC_0.25[which.max(NUC_0.25$negative_log10_pvalue),],
                      NUC_0.30[which.max(NUC_0.30$negative_log10_pvalue),],
                      NUC_0.35[which.max(NUC_0.35$negative_log10_pvalue),],
                      NUC_0.40[which.max(NUC_0.40$negative_log10_pvalue),],
                      NUC_0.45[which.max(NUC_0.45$negative_log10_pvalue),],
                      NUC_0.50[which.max(NUC_0.50$negative_log10_pvalue),]
)

# CB stats (for top cluster)
# Create a table with the total number and 
CB_0.00 <- df_0_curated[df_0_curated$term_name == "Cajal body curated",]
CB_0.00$negative_log10_pvalue <- -log10(CB_0.00$adjusted_p_value)
CB_0.00$Number_clusters_with_term <- nrow(CB_0.00[CB_0.00$significant == TRUE,])

CB_0.05 <- df_0.05_curated[df_0.05_curated$term_name == "Cajal body curated",]
CB_0.05$negative_log10_pvalue <- -log10(CB_0.05$adjusted_p_value)
CB_0.05$Number_clusters_with_term <- nrow(CB_0.05[CB_0.05$significant == TRUE,])

CB_0.10 <- df_0.10_curated[df_0.10_curated$term_name == "Cajal body curated",]
CB_0.10$negative_log10_pvalue <- -log10(CB_0.10$adjusted_p_value)
CB_0.10$Number_clusters_with_term <- nrow(CB_0.10[CB_0.10$significant == TRUE,])

CB_0.15 <- df_0.15_curated[df_0.15_curated$term_name == "Cajal body curated",]
CB_0.15$negative_log10_pvalue <- -log10(CB_0.15$adjusted_p_value)
CB_0.15$Number_clusters_with_term <- nrow(CB_0.15[CB_0.15$significant == TRUE,])

CB_0.20 <- df_0.20_curated[df_0.20_curated$term_name == "Cajal body curated",]
CB_0.20$negative_log10_pvalue <- -log10(CB_0.20$adjusted_p_value)
CB_0.20$Number_clusters_with_term <- nrow(CB_0.20[CB_0.20$significant == TRUE,])

CB_0.25 <- df_0.25_curated[df_0.25_curated$term_name == "Cajal body curated",]
CB_0.25$negative_log10_pvalue <- -log10(CB_0.25$adjusted_p_value)
CB_0.25$Number_clusters_with_term <- nrow(CB_0.25[CB_0.25$significant == TRUE,])

CB_0.30 <- df_0.30_curated[df_0.30_curated$term_name == "Cajal body curated",]
CB_0.30$negative_log10_pvalue <- -log10(CB_0.30$adjusted_p_value)
CB_0.30$Number_clusters_with_term <- nrow(CB_0.30[CB_0.30$significant == TRUE,])

CB_0.35 <- df_0.35_curated[df_0.35_curated$term_name == "Cajal body curated",]
CB_0.35$negative_log10_pvalue <- -log10(CB_0.35$adjusted_p_value)
CB_0.35$Number_clusters_with_term <- nrow(CB_0.35[CB_0.35$significant == TRUE,])

CB_0.40 <- df_0.40_curated[df_0.40_curated$term_name == "Cajal body curated",]
CB_0.40$negative_log10_pvalue <- -log10(CB_0.40$adjusted_p_value)
CB_0.40$Number_clusters_with_term <- nrow(CB_0.40[CB_0.40$significant == TRUE,])

CB_0.45 <- df_0.45_curated[df_0.45_curated$term_name == "Cajal body curated",]
CB_0.45$negative_log10_pvalue <- -log10(CB_0.45$adjusted_p_value)
CB_0.45$Number_clusters_with_term <- nrow(CB_0.45[CB_0.45$significant == TRUE,])

CB_0.50 <- df_0.50_curated[df_0.50_curated$term_name == "Cajal body curated",]
CB_0.50$negative_log10_pvalue <- -log10(CB_0.50$adjusted_p_value)
CB_0.50$Number_clusters_with_term <- nrow(CB_0.50[CB_0.50$significant == TRUE,])

top_rank_CB <- rbind(CB_0.00[which.max(CB_0.00$negative_log10_pvalue),],
                     CB_0.05[which.max(CB_0.05$negative_log10_pvalue),],
                     CB_0.10[which.max(CB_0.10$negative_log10_pvalue),],
                     CB_0.15[which.max(CB_0.15$negative_log10_pvalue),],
                     CB_0.20[which.max(CB_0.20$negative_log10_pvalue),],
                     CB_0.25[which.max(CB_0.25$negative_log10_pvalue),],
                     CB_0.30[which.max(CB_0.30$negative_log10_pvalue),],
                     CB_0.35[which.max(CB_0.35$negative_log10_pvalue),],
                     CB_0.40[which.max(CB_0.40$negative_log10_pvalue),],
                     CB_0.45[which.max(CB_0.45$negative_log10_pvalue),],
                     CB_0.50[which.max(CB_0.50$negative_log10_pvalue),]
)



# PML stats (for top cluster)
# Create a table with the total number and 
PML_0.00 <- df_0_curated[df_0_curated$term_name == "PML body curated",]
PML_0.00$negative_log10_pvalue <- -log10(PML_0.00$adjusted_p_value)
PML_0.00$Number_clusters_with_term <- nrow(PML_0.00[PML_0.00$significant == TRUE,])

PML_0.05 <- df_0.05_curated[df_0.05_curated$term_name == "PML body curated",]
PML_0.05$negative_log10_pvalue <- -log10(PML_0.05$adjusted_p_value)
PML_0.05$Number_clusters_with_term <- nrow(PML_0.05[PML_0.05$significant == TRUE,])

PML_0.10 <- df_0.10_curated[df_0.10_curated$term_name == "PML body curated",]
PML_0.10$negative_log10_pvalue <- -log10(PML_0.10$adjusted_p_value)
PML_0.10$Number_clusters_with_term <- nrow(PML_0.10[PML_0.10$significant == TRUE,])

PML_0.15 <- df_0.15_curated[df_0.15_curated$term_name == "PML body curated",]
PML_0.15$negative_log10_pvalue <- -log10(PML_0.15$adjusted_p_value)
PML_0.15$Number_clusters_with_term <- nrow(PML_0.15[PML_0.15$significant == TRUE,])

PML_0.20 <- df_0.20_curated[df_0.20_curated$term_name == "PML body curated",]
PML_0.20$negative_log10_pvalue <- -log10(PML_0.20$adjusted_p_value)
PML_0.20$Number_clusters_with_term <- nrow(PML_0.20[PML_0.20$significant == TRUE,])

PML_0.25 <- df_0.25_curated[df_0.25_curated$term_name == "PML body curated",]
PML_0.25$negative_log10_pvalue <- -log10(PML_0.25$adjusted_p_value)
PML_0.25$Number_clusters_with_term <- nrow(PML_0.25[PML_0.25$significant == TRUE,])

PML_0.30 <- df_0.30_curated[df_0.30_curated$term_name == "PML body curated",]
PML_0.30$negative_log10_pvalue <- -log10(PML_0.30$adjusted_p_value)
PML_0.30$Number_clusters_with_term <- nrow(PML_0.30[PML_0.30$significant == TRUE,])

PML_0.35 <- df_0.35_curated[df_0.35_curated$term_name == "PML body curated",]
PML_0.35$negative_log10_pvalue <- -log10(PML_0.35$adjusted_p_value)
PML_0.35$Number_clusters_with_term <- nrow(PML_0.35[PML_0.35$significant == TRUE,])

PML_0.40 <- df_0.40_curated[df_0.40_curated$term_name == "PML body curated",]
PML_0.40$negative_log10_pvalue <- -log10(PML_0.40$adjusted_p_value)
PML_0.40$Number_clusters_with_term <- nrow(PML_0.40[PML_0.40$significant == TRUE,])

PML_0.45 <- df_0.45_curated[df_0.45_curated$term_name == "PML body curated",]
PML_0.45$negative_log10_pvalue <- -log10(PML_0.45$adjusted_p_value)
PML_0.45$Number_clusters_with_term <- nrow(PML_0.45[PML_0.45$significant == TRUE,])

PML_0.50 <- df_0.50_curated[df_0.50_curated$term_name == "PML body curated",]
PML_0.50$negative_log10_pvalue <- -log10(PML_0.50$adjusted_p_value)
PML_0.50$Number_clusters_with_term <- nrow(PML_0.50[PML_0.50$significant == TRUE,])

top_rank_PML <- rbind(PML_0.00[which.max(PML_0.00$negative_log10_pvalue),],
                      PML_0.05[which.max(PML_0.05$negative_log10_pvalue),],
                      PML_0.10[which.max(PML_0.10$negative_log10_pvalue),],
                      PML_0.15[which.max(PML_0.15$negative_log10_pvalue),],
                      PML_0.20[which.max(PML_0.20$negative_log10_pvalue),],
                      PML_0.25[which.max(PML_0.25$negative_log10_pvalue),],
                      PML_0.30[which.max(PML_0.30$negative_log10_pvalue),],
                      PML_0.35[which.max(PML_0.35$negative_log10_pvalue),],
                      PML_0.40[which.max(PML_0.40$negative_log10_pvalue),],
                      PML_0.45[which.max(PML_0.45$negative_log10_pvalue),],
                      PML_0.50[which.max(PML_0.50$negative_log10_pvalue),]
)


# Combine all top rank data into one file
df_combined_curated_top <- rbind(top_rank_PS,
                                 top_rank_NS,
                                 top_rank_NUC,
                                 top_rank_CB,
                                 top_rank_PML
)

write.table(df_combined_curated_top, "3_gProfile_Results_NMF19_curated_terms_combined_top10ctls_0.70_multi_percentfilter_toprank.txt", sep="\t", row.names=FALSE)

#--------
# Plot PS data
library(ggplot2)
library(ggrepel)

#Negative log10 p-values
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_pValue_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=negative_log10_pvalue, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Negative log10 pvalues) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("-log10(pvalue) (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,150))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Recall
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_Recall_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=recall, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Recall) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("recall (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,0.5))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Precision
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_Precision_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=precision, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Precision) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("precision (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,1.0))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Query size
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_QuerySize_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=query_size, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Query Size) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("Query size (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,200))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F1scores
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_F1score_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F1score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:F1-score)")+
  xlab("NMF percent column max score cutoff")+
  ylab("F1-score (2*Prec*Recall)/(Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F0.5 scores
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_F0.5score_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F0.5score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:F0.5-score)")+
  xlab("NMF percent column max score cutoff")+
  ylab("F0.5-score ((1+0.5^2)*Prec*Recall)/(0.5^2*Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F2scores
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_F2score_top10ctls_0.70_multi.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F2score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:F2-score)")+
  xlab("NMF percent column max score cutoff")+
  ylab("F2-score ((1+2^2)*Prec*Recall)/(2^2*Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

#F1score combined
F1scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Primary_cutoff)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Primary_cutoff == i,]
  df_fscores$average_F1score <- mean(df_fscores$F1score)
  avgFscore <- df_fscores[1,c(1,23)]
  F1scores <- rbind(avgFscore, F1scores)
}

# F1score average
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_avgF1score_top10ctls_0.70_multi.pdf", width = 7, height = 6)
ggplot(data=F1scores, aes(x=Primary_cutoff, y=average_F1score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("Average F1score (2*Prec*Recall)/(Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

#F0.5score combined
F0.5scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Primary_cutoff)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Primary_cutoff == i,]
  df_fscores$average_F0.5score <- mean(df_fscores$F0.5score)
  avgFscore <- df_fscores[1,c(1,23)]
  F0.5scores <- rbind(avgFscore, F0.5scores)
}

# F0.5score average
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_avgF0.5score_top10ctls_0.70_multi.pdf", width = 7, height = 6)
ggplot(data=F0.5scores, aes(x=Primary_cutoff, y=average_F0.5score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("Average F2score F0.5-score ((1+0.5^2)*Prec*Recall)/(0.5^2*Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

#F2score combined
F2scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Primary_cutoff)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Primary_cutoff == i,]
  df_fscores$average_F2score <- mean(df_fscores$F2score)
  avgFscore <- df_fscores[1,c(1,23)]
  F2scores <- rbind(avgFscore, F2scores)
}

# F2score average
pdf("Rplot_NMF_19ranks_score_percent_cutoff_testing_gProfiler_TopPvalueRank_avgF2score_top10ctls_0.70_multi.pdf", width = 7, height = 6)
ggplot(data=F2scores, aes(x=Primary_cutoff, y=average_F2score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF percent cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF percent column max score cutoff")+
  ylab("Average F2score F2-score ((1+2^2)*Prec*Recall)/(2^2*Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

Fscores_avg_combined <- cbind(F0.5scores,F1scores,F2scores)
Fscores_avg_combined <- Fscores_avg_combined[,c(1,2,4,6)]
write.table(Fscores_avg_combined, "4_gProfile_Results_NMF19_curated_terms_top10ctls_percentfilter_0.70_multi_toprank_avgFscores.txt", sep="\t", row.names=FALSE)
