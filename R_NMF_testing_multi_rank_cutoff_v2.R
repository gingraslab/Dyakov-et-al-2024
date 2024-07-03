# Outputting Precision/Recall and Enrichment plots to assess different NMF multi-rank cutoffs (multi)

##Part 1: Create columns with primary, secondary and tertiary ranks re-annotated with different multi-rank filters
##Part 2: Run gProfiler for each of the samples using Boris' custom terms
##Part 3: Calculate precision/recall and enrichments and plot across clusters

#---------
##Part 1: Create columns with primary, secondary and tertiary ranks re-annotated with different multi-rank filters
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

write.table(NMF_matrix, "0.5_nmf_prey_matrix_input_19ranks.csv", sep=",", row.names=FALSE)

#-----------
# Create NMF score matrices filtered by different multi rank filters for second/tertiary
# Matrix_0
NMF_matrix_0_multi_0.00 <- NMF_matrix
NMF_matrix_0_multi_0.00$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.00$Primary_score > 0, NMF_matrix_0_multi_0.00$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.00$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.00$Primary_score > 0 & NMF_matrix_0_multi_0.00$Secondary_score > 0 & NMF_matrix_0_multi_0.00$Multi_loc_score_secondary >= 0.00, NMF_matrix_0_multi_0.00$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.00$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.00$Primary_score > 0 & NMF_matrix_0_multi_0.00$Tertiary_score > 0 & NMF_matrix_0_multi_0.00$Multi_loc_score_tertiary >= 0.00, NMF_matrix_0_multi_0.00$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.00$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.00$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.00$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.00$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.00, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.05
NMF_matrix_0_multi_0.05 <- NMF_matrix
NMF_matrix_0_multi_0.05$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.05$Primary_score > 0, NMF_matrix_0_multi_0.05$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.05$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.05$Primary_score > 0 & NMF_matrix_0_multi_0.05$Secondary_score > 0 & NMF_matrix_0_multi_0.05$Multi_loc_score_secondary >= 0.05, NMF_matrix_0_multi_0.05$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.05$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.05$Primary_score > 0 & NMF_matrix_0_multi_0.05$Tertiary_score > 0 & NMF_matrix_0_multi_0.05$Multi_loc_score_tertiary >= 0.05, NMF_matrix_0_multi_0.05$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.05$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.05$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.05$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.05$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.05, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.05_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.10
NMF_matrix_0_multi_0.10 <- NMF_matrix
NMF_matrix_0_multi_0.10$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.10$Primary_score > 0, NMF_matrix_0_multi_0.10$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.10$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.10$Primary_score > 0 & NMF_matrix_0_multi_0.10$Secondary_score > 0 & NMF_matrix_0_multi_0.10$Multi_loc_score_secondary >= 0.10, NMF_matrix_0_multi_0.10$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.10$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.10$Primary_score > 0 & NMF_matrix_0_multi_0.10$Tertiary_score > 0 & NMF_matrix_0_multi_0.10$Multi_loc_score_tertiary >= 0.10, NMF_matrix_0_multi_0.10$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.10$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.10$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.10$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.10$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.10, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.10_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.15
NMF_matrix_0_multi_0.15 <- NMF_matrix
NMF_matrix_0_multi_0.15$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.15$Primary_score > 0, NMF_matrix_0_multi_0.15$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.15$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.15$Primary_score > 0 & NMF_matrix_0_multi_0.15$Secondary_score > 0 & NMF_matrix_0_multi_0.15$Multi_loc_score_secondary >= 0.15, NMF_matrix_0_multi_0.15$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.15$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.15$Primary_score > 0 & NMF_matrix_0_multi_0.15$Tertiary_score > 0 & NMF_matrix_0_multi_0.15$Multi_loc_score_tertiary >= 0.15, NMF_matrix_0_multi_0.15$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.15$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.15$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.15$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.15$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.15, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.15_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.20
NMF_matrix_0_multi_0.20 <- NMF_matrix
NMF_matrix_0_multi_0.20$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.20$Primary_score > 0, NMF_matrix_0_multi_0.20$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.20$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.20$Primary_score > 0 & NMF_matrix_0_multi_0.20$Secondary_score > 0 & NMF_matrix_0_multi_0.20$Multi_loc_score_secondary >= 0.20, NMF_matrix_0_multi_0.20$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.20$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.20$Primary_score > 0 & NMF_matrix_0_multi_0.20$Tertiary_score > 0 & NMF_matrix_0_multi_0.20$Multi_loc_score_tertiary >= 0.20, NMF_matrix_0_multi_0.20$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.20$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.20$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.20$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.20$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.20, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.20_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.25
NMF_matrix_0_multi_0.25 <- NMF_matrix
NMF_matrix_0_multi_0.25$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.25$Primary_score > 0, NMF_matrix_0_multi_0.25$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.25$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.25$Primary_score > 0 & NMF_matrix_0_multi_0.25$Secondary_score > 0 & NMF_matrix_0_multi_0.25$Multi_loc_score_secondary >= 0.25, NMF_matrix_0_multi_0.25$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.25$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.25$Primary_score > 0 & NMF_matrix_0_multi_0.25$Tertiary_score > 0 & NMF_matrix_0_multi_0.25$Multi_loc_score_tertiary >= 0.25, NMF_matrix_0_multi_0.25$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.25$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.25$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.25$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.25$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.25, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.25_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.30
NMF_matrix_0_multi_0.30 <- NMF_matrix
NMF_matrix_0_multi_0.30$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.30$Primary_score > 0, NMF_matrix_0_multi_0.30$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.30$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.30$Primary_score > 0 & NMF_matrix_0_multi_0.30$Secondary_score > 0 & NMF_matrix_0_multi_0.30$Multi_loc_score_secondary >= 0.30, NMF_matrix_0_multi_0.30$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.30$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.30$Primary_score > 0 & NMF_matrix_0_multi_0.30$Tertiary_score > 0 & NMF_matrix_0_multi_0.30$Multi_loc_score_tertiary >= 0.30, NMF_matrix_0_multi_0.30$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.30$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.30$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.30$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.30$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.30, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.30_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.35
NMF_matrix_0_multi_0.35 <- NMF_matrix
NMF_matrix_0_multi_0.35$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.35$Primary_score > 0, NMF_matrix_0_multi_0.35$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.35$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.35$Primary_score > 0 & NMF_matrix_0_multi_0.35$Secondary_score > 0 & NMF_matrix_0_multi_0.35$Multi_loc_score_secondary >= 0.35, NMF_matrix_0_multi_0.35$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.35$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.35$Primary_score > 0 & NMF_matrix_0_multi_0.35$Tertiary_score > 0 & NMF_matrix_0_multi_0.35$Multi_loc_score_tertiary >= 0.35, NMF_matrix_0_multi_0.35$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.35$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.35$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.35$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.35$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.35, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.35_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.40
NMF_matrix_0_multi_0.40 <- NMF_matrix
NMF_matrix_0_multi_0.40$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.40$Primary_score > 0, NMF_matrix_0_multi_0.40$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.40$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.40$Primary_score > 0 & NMF_matrix_0_multi_0.40$Secondary_score > 0 & NMF_matrix_0_multi_0.40$Multi_loc_score_secondary >= 0.40, NMF_matrix_0_multi_0.40$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.40$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.40$Primary_score > 0 & NMF_matrix_0_multi_0.40$Tertiary_score > 0 & NMF_matrix_0_multi_0.40$Multi_loc_score_tertiary >= 0.40, NMF_matrix_0_multi_0.40$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.40$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.40$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.40$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.40$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.40, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.40_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.45
NMF_matrix_0_multi_0.45 <- NMF_matrix
NMF_matrix_0_multi_0.45$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.45$Primary_score > 0, NMF_matrix_0_multi_0.45$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.45$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.45$Primary_score > 0 & NMF_matrix_0_multi_0.45$Secondary_score > 0 & NMF_matrix_0_multi_0.45$Multi_loc_score_secondary >= 0.45, NMF_matrix_0_multi_0.45$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.45$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.45$Primary_score > 0 & NMF_matrix_0_multi_0.45$Tertiary_score > 0 & NMF_matrix_0_multi_0.45$Multi_loc_score_tertiary >= 0.45, NMF_matrix_0_multi_0.45$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.45$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.45$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.45$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.45$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.45, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.45_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.50
NMF_matrix_0_multi_0.50 <- NMF_matrix
NMF_matrix_0_multi_0.50$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.50$Primary_score > 0, NMF_matrix_0_multi_0.50$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.50$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.50$Primary_score > 0 & NMF_matrix_0_multi_0.50$Secondary_score > 0 & NMF_matrix_0_multi_0.50$Multi_loc_score_secondary >= 0.50, NMF_matrix_0_multi_0.50$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.50$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.50$Primary_score > 0 & NMF_matrix_0_multi_0.50$Tertiary_score > 0 & NMF_matrix_0_multi_0.50$Multi_loc_score_tertiary >= 0.50, NMF_matrix_0_multi_0.50$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.50$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.50$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.50$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.50$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.50, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.50_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.55
NMF_matrix_0_multi_0.55 <- NMF_matrix
NMF_matrix_0_multi_0.55$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.55$Primary_score > 0, NMF_matrix_0_multi_0.55$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.55$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.55$Primary_score > 0 & NMF_matrix_0_multi_0.55$Secondary_score > 0 & NMF_matrix_0_multi_0.55$Multi_loc_score_secondary >= 0.55, NMF_matrix_0_multi_0.55$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.55$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.55$Primary_score > 0 & NMF_matrix_0_multi_0.55$Tertiary_score > 0 & NMF_matrix_0_multi_0.55$Multi_loc_score_tertiary >= 0.55, NMF_matrix_0_multi_0.55$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.55$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.55$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.55$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.55$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.55, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.55_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.60
NMF_matrix_0_multi_0.60 <- NMF_matrix
NMF_matrix_0_multi_0.60$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.60$Primary_score > 0, NMF_matrix_0_multi_0.60$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.60$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.60$Primary_score > 0 & NMF_matrix_0_multi_0.60$Secondary_score > 0 & NMF_matrix_0_multi_0.60$Multi_loc_score_secondary >= 0.60, NMF_matrix_0_multi_0.60$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.60$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.60$Primary_score > 0 & NMF_matrix_0_multi_0.60$Tertiary_score > 0 & NMF_matrix_0_multi_0.60$Multi_loc_score_tertiary >= 0.60, NMF_matrix_0_multi_0.60$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.60$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.60$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.60$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.60$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.60, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.60_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.65
NMF_matrix_0_multi_0.65 <- NMF_matrix
NMF_matrix_0_multi_0.65$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.65$Primary_score > 0, NMF_matrix_0_multi_0.65$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.65$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.65$Primary_score > 0 & NMF_matrix_0_multi_0.65$Secondary_score > 0 & NMF_matrix_0_multi_0.65$Multi_loc_score_secondary >= 0.65, NMF_matrix_0_multi_0.65$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.65$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.65$Primary_score > 0 & NMF_matrix_0_multi_0.65$Tertiary_score > 0 & NMF_matrix_0_multi_0.65$Multi_loc_score_tertiary >= 0.65, NMF_matrix_0_multi_0.65$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.65$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.65$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.65$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.65$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.65, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.65_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.70
NMF_matrix_0_multi_0.70 <- NMF_matrix
NMF_matrix_0_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.70$Primary_score > 0, NMF_matrix_0_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.70$Primary_score > 0 & NMF_matrix_0_multi_0.70$Secondary_score > 0 & NMF_matrix_0_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_0_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.70$Primary_score > 0 & NMF_matrix_0_multi_0.70$Tertiary_score > 0 & NMF_matrix_0_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_0_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.70$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.70, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.70_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.75
NMF_matrix_0_multi_0.75 <- NMF_matrix
NMF_matrix_0_multi_0.75$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.75$Primary_score > 0, NMF_matrix_0_multi_0.75$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.75$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.75$Primary_score > 0 & NMF_matrix_0_multi_0.75$Secondary_score > 0 & NMF_matrix_0_multi_0.75$Multi_loc_score_secondary >= 0.75, NMF_matrix_0_multi_0.75$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.75$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.75$Primary_score > 0 & NMF_matrix_0_multi_0.75$Tertiary_score > 0 & NMF_matrix_0_multi_0.75$Multi_loc_score_tertiary >= 0.75, NMF_matrix_0_multi_0.75$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.75$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.75$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.75$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.75$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.75, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.75_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.80
NMF_matrix_0_multi_0.80 <- NMF_matrix
NMF_matrix_0_multi_0.80$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.80$Primary_score > 0, NMF_matrix_0_multi_0.80$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.80$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.80$Primary_score > 0 & NMF_matrix_0_multi_0.80$Secondary_score > 0 & NMF_matrix_0_multi_0.80$Multi_loc_score_secondary >= 0.80, NMF_matrix_0_multi_0.80$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.80$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.80$Primary_score > 0 & NMF_matrix_0_multi_0.80$Tertiary_score > 0 & NMF_matrix_0_multi_0.80$Multi_loc_score_tertiary >= 0.80, NMF_matrix_0_multi_0.80$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.80$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.80$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.80$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.80$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.80, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.80_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.85
NMF_matrix_0_multi_0.85 <- NMF_matrix
NMF_matrix_0_multi_0.85$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.85$Primary_score > 0, NMF_matrix_0_multi_0.85$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.85$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.85$Primary_score > 0 & NMF_matrix_0_multi_0.85$Secondary_score > 0 & NMF_matrix_0_multi_0.85$Multi_loc_score_secondary >= 0.85, NMF_matrix_0_multi_0.85$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.85$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.85$Primary_score > 0 & NMF_matrix_0_multi_0.85$Tertiary_score > 0 & NMF_matrix_0_multi_0.85$Multi_loc_score_tertiary >= 0.85, NMF_matrix_0_multi_0.85$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.85$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.85$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.85$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.85$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.85, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.85_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.90
NMF_matrix_0_multi_0.90 <- NMF_matrix
NMF_matrix_0_multi_0.90$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.90$Primary_score > 0, NMF_matrix_0_multi_0.90$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.90$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.90$Primary_score > 0 & NMF_matrix_0_multi_0.90$Secondary_score > 0 & NMF_matrix_0_multi_0.90$Multi_loc_score_secondary >= 0.90, NMF_matrix_0_multi_0.90$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.90$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.90$Primary_score > 0 & NMF_matrix_0_multi_0.90$Tertiary_score > 0 & NMF_matrix_0_multi_0.90$Multi_loc_score_tertiary >= 0.90, NMF_matrix_0_multi_0.90$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.90$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.90$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.90$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.90$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.90, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.90_multi.csv", sep=",", row.names=FALSE)

# Matrix_0.95
NMF_matrix_0_multi_0.95 <- NMF_matrix
NMF_matrix_0_multi_0.95$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.95$Primary_score > 0, NMF_matrix_0_multi_0.95$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_0.95$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.95$Primary_score > 0 & NMF_matrix_0_multi_0.95$Secondary_score > 0 & NMF_matrix_0_multi_0.95$Multi_loc_score_secondary >= 0.95, NMF_matrix_0_multi_0.95$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_0.95$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_0.95$Primary_score > 0 & NMF_matrix_0_multi_0.95$Tertiary_score > 0 & NMF_matrix_0_multi_0.95$Multi_loc_score_tertiary >= 0.95, NMF_matrix_0_multi_0.95$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_0.95$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_0.95$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_0.95$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_0.95$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_0.95, "0.5_nmf_prey_matrix_input_19ranks_0_primary_0.95_multi.csv", sep=",", row.names=FALSE)

# Matrix_1.00
NMF_matrix_0_multi_1.00 <- NMF_matrix
NMF_matrix_0_multi_1.00$Primary_cluster_annotation <- ifelse(NMF_matrix_0_multi_1.00$Primary_score > 0, NMF_matrix_0_multi_1.00$Primary_cluster_annotation, NA)
NMF_matrix_0_multi_1.00$Secondary_cluster_annotation <- ifelse(NMF_matrix_0_multi_1.00$Primary_score > 0 & NMF_matrix_0_multi_1.00$Secondary_score > 0 & NMF_matrix_0_multi_1.00$Multi_loc_score_secondary >= 1.00, NMF_matrix_0_multi_1.00$Secondary_cluster_annotation, NA)
NMF_matrix_0_multi_1.00$Tertiary_cluster_annotation <- ifelse(NMF_matrix_0_multi_1.00$Primary_score > 0 & NMF_matrix_0_multi_1.00$Tertiary_score > 0 & NMF_matrix_0_multi_1.00$Multi_loc_score_tertiary >= 1.00, NMF_matrix_0_multi_1.00$Tertiary_cluster_annotation, NA)
NMF_matrix_0_multi_1.00$Multi_rank_cutoff <- 0
NMF_matrix_0_multi_1.00$Multi_rank <- ifelse(is.na(NMF_matrix_0_multi_1.00$Secondary_cluster_annotation) & is.na(NMF_matrix_0_multi_1.00$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_0_multi_1.00, "0.5_nmf_prey_matrix_input_19ranks_0_primary_1.00_multi.csv", sep=",", row.names=FALSE)


#-------
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
  query <- NMF_matrix_0_multi_0.00[which(NMF_matrix_0_multi_0.00$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.00$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.00$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0, 'gprofiler_Results_NMF19_customToken_top10ctls_0_NMFcutoff.xlsx', overwrite = TRUE)
df_0$parents <- as.character(df_0$parents)
df_0$Primary_cutoff <- 0
df_0 <- df_0[,c(17,1:16)]
df_0_curated <- df_0[df_0$term_name %in% curated_terms,]
df_0_curated$F1score <- (2*(df_0_curated$precision*df_0_curated$recall)) / (df_0_curated$precision + df_0_curated$recall)
df_0_curated$F0.5score <- ((1 + 0.5^2) * df_0_curated$precision * df_0_curated$recall) / (0.5^2 * df_0_curated$precision + df_0_curated$recall)
df_0_curated$F2score <- ((1 + 2^2) * df_0_curated$precision * df_0_curated$recall) / (2^2 * df_0_curated$precision + df_0_curated$recall)

write.table(df_0, "1_gProfile_Results_NMF19_customToken_top10ctls_0_multirank.txt", sep="\t", row.names=FALSE)


# 0.05 Filter
df_0.05 <- data.frame()
wb_0.05 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.05[which(NMF_matrix_0_multi_0.05$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.05$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.05$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.05, 'gprofiler_Results_NMF19_customToken_top10ctls_0.05_multirank.xlsx', overwrite = TRUE)
df_0.05$parents <- as.character(df_0.05$parents)
df_0.05$Primary_cutoff <- 0.05
df_0.05 <- df_0.05[,c(17,1:16)]
df_0.05_curated <- df_0.05[df_0.05$term_name %in% curated_terms,]
df_0.05_curated$F1score <- (2*(df_0.05_curated$precision*df_0.05_curated$recall)) / (df_0.05_curated$precision + df_0.05_curated$recall)
df_0.05_curated$F0.5score <- ((1 + 0.5^2) * df_0.05_curated$precision * df_0.05_curated$recall) / (0.5^2 * df_0.05_curated$precision + df_0.05_curated$recall)
df_0.05_curated$F2score <- ((1 + 2^2) * df_0.05_curated$precision * df_0.05_curated$recall) / (2^2 * df_0.05_curated$precision + df_0.05_curated$recall)

write.table(df_0.05, "1_gProfile_Results_NMF19_customToken_top10ctls_0.05_multirank.txt", sep="\t", row.names=FALSE)


# 0.10 Filter
df_0.10 <- data.frame()
wb_0.10 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.10[which(NMF_matrix_0_multi_0.10$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.10$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.10$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.10, 'gprofiler_Results_NMF19_customToken_top10ctls_0.10_multirank.xlsx', overwrite = TRUE)
df_0.10$parents <- as.character(df_0.10$parents)
df_0.10$Primary_cutoff <- 0.10
df_0.10 <- df_0.10[,c(17,1:16)]
df_0.10_curated <- df_0.10[df_0.10$term_name %in% curated_terms,]
df_0.10_curated$F1score <- (2*(df_0.10_curated$precision*df_0.10_curated$recall)) / (df_0.10_curated$precision + df_0.10_curated$recall)
df_0.10_curated$F0.5score <- ((1 + 0.5^2) * df_0.10_curated$precision * df_0.10_curated$recall) / (0.5^2 * df_0.10_curated$precision + df_0.10_curated$recall)
df_0.10_curated$F2score <- ((1 + 2^2) * df_0.10_curated$precision * df_0.10_curated$recall) / (2^2 * df_0.10_curated$precision + df_0.10_curated$recall)

write.table(df_0.10, "1_gProfile_Results_NMF19_customToken_top10ctls_0.10_multirank.txt", sep="\t", row.names=FALSE)


# 0.15 Filter
df_0.15 <- data.frame()
wb_0.15 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.15[which(NMF_matrix_0_multi_0.15$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.15$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.15$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.15, 'gprofiler_Results_NMF19_customToken_top10ctls_0.15_multirank.xlsx', overwrite = TRUE)
df_0.15$parents <- as.character(df_0.15$parents)
df_0.15$Primary_cutoff <- 0.15
df_0.15 <- df_0.15[,c(17,1:16)]
df_0.15_curated <- df_0.15[df_0.15$term_name %in% curated_terms,]
df_0.15_curated$F1score <- (2*(df_0.15_curated$precision*df_0.15_curated$recall)) / (df_0.15_curated$precision + df_0.15_curated$recall)
df_0.15_curated$F0.5score <- ((1 + 0.5^2) * df_0.15_curated$precision * df_0.15_curated$recall) / (0.5^2 * df_0.15_curated$precision + df_0.15_curated$recall)
df_0.15_curated$F2score <- ((1 + 2^2) * df_0.15_curated$precision * df_0.15_curated$recall) / (2^2 * df_0.15_curated$precision + df_0.15_curated$recall)

write.table(df_0.15, "1_gProfile_Results_NMF19_customToken_top10ctls_0.15_multirank.txt", sep="\t", row.names=FALSE)


# 0.20 Filter
df_0.20 <- data.frame()
wb_0.20 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.20[which(NMF_matrix_0_multi_0.20$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.20$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.20$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.20, 'gprofiler_Results_NMF19_customToken_top10ctls_0.20_multirank.xlsx', overwrite = TRUE)
df_0.20$parents <- as.character(df_0.20$parents)
df_0.20$Primary_cutoff <- 0.20
df_0.20 <- df_0.20[,c(17,1:16)]
df_0.20_curated <- df_0.20[df_0.20$term_name %in% curated_terms,]
df_0.20_curated$F1score <- (2*(df_0.20_curated$precision*df_0.20_curated$recall)) / (df_0.20_curated$precision + df_0.20_curated$recall)
df_0.20_curated$F0.5score <- ((1 + 0.5^2) * df_0.20_curated$precision * df_0.20_curated$recall) / (0.5^2 * df_0.20_curated$precision + df_0.20_curated$recall)
df_0.20_curated$F2score <- ((1 + 2^2) * df_0.20_curated$precision * df_0.20_curated$recall) / (2^2 * df_0.20_curated$precision + df_0.20_curated$recall)

write.table(df_0.20, "1_gProfile_Results_NMF19_customToken_top10ctls_0.20_multirank.txt", sep="\t", row.names=FALSE)


# 0.25 Filter
df_0.25 <- data.frame()
wb_0.25 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.25[which(NMF_matrix_0_multi_0.25$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.25$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.25$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.25, 'gprofiler_Results_NMF19_customToken_top10ctls_0.25_multirank.xlsx', overwrite = TRUE)
df_0.25$parents <- as.character(df_0.25$parents)
df_0.25$Primary_cutoff <- 0.25
df_0.25 <- df_0.25[,c(17,1:16)]
df_0.25_curated <- df_0.25[df_0.25$term_name %in% curated_terms,]
df_0.25_curated$F1score <- (2*(df_0.25_curated$precision*df_0.25_curated$recall)) / (df_0.25_curated$precision + df_0.25_curated$recall)
df_0.25_curated$F0.5score <- ((1 + 0.5^2) * df_0.25_curated$precision * df_0.25_curated$recall) / (0.5^2 * df_0.25_curated$precision + df_0.25_curated$recall)
df_0.25_curated$F2score <- ((1 + 2^2) * df_0.25_curated$precision * df_0.25_curated$recall) / (2^2 * df_0.25_curated$precision + df_0.25_curated$recall)

write.table(df_0.25, "1_gProfile_Results_NMF19_customToken_top10ctls_0.25_multirank.txt", sep="\t", row.names=FALSE)


# 0.30 Filter
df_0.30 <- data.frame()
wb_0.30 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.30[which(NMF_matrix_0_multi_0.30$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.30$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.30$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.30, 'gprofiler_Results_NMF19_customToken_top10ctls_0.30_multirank.xlsx', overwrite = TRUE)
df_0.30$parents <- as.character(df_0.30$parents)
df_0.30$Primary_cutoff <- 0.30
df_0.30 <- df_0.30[,c(17,1:16)]
df_0.30_curated <- df_0.30[df_0.30$term_name %in% curated_terms,]
df_0.30_curated$F1score <- (2*(df_0.30_curated$precision*df_0.30_curated$recall)) / (df_0.30_curated$precision + df_0.30_curated$recall)
df_0.30_curated$F0.5score <- ((1 + 0.5^2) * df_0.30_curated$precision * df_0.30_curated$recall) / (0.5^2 * df_0.30_curated$precision + df_0.30_curated$recall)
df_0.30_curated$F2score <- ((1 + 2^2) * df_0.30_curated$precision * df_0.30_curated$recall) / (2^2 * df_0.30_curated$precision + df_0.30_curated$recall)

write.table(df_0.30, "1_gProfile_Results_NMF19_customToken_top10ctls_0.30_multirank.txt", sep="\t", row.names=FALSE)


# 0.35 Filter
df_0.35 <- data.frame()
wb_0.35 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.35[which(NMF_matrix_0_multi_0.35$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.35$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.35$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.35, 'gprofiler_Results_NMF19_customToken_top10ctls_0.35_multirank.xlsx', overwrite = TRUE)
df_0.35$parents <- as.character(df_0.35$parents)
df_0.35$Primary_cutoff <- 0.35
df_0.35 <- df_0.35[,c(17,1:16)]
df_0.35_curated <- df_0.35[df_0.35$term_name %in% curated_terms,]
df_0.35_curated$F1score <- (2*(df_0.35_curated$precision*df_0.35_curated$recall)) / (df_0.35_curated$precision + df_0.35_curated$recall)
df_0.35_curated$F0.5score <- ((1 + 0.5^2) * df_0.35_curated$precision * df_0.35_curated$recall) / (0.5^2 * df_0.35_curated$precision + df_0.35_curated$recall)
df_0.35_curated$F2score <- ((1 + 2^2) * df_0.35_curated$precision * df_0.35_curated$recall) / (2^2 * df_0.35_curated$precision + df_0.35_curated$recall)

write.table(df_0.35, "1_gProfile_Results_NMF19_customToken_top10ctls_0.35_multirank.txt", sep="\t", row.names=FALSE)


# 0.40 Filter
df_0.40 <- data.frame()
wb_0.40 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.40[which(NMF_matrix_0_multi_0.40$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.40$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.40$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.40, 'gprofiler_Results_NMF19_customToken_top10ctls_0.40_multirank.xlsx', overwrite = TRUE)
df_0.40$parents <- as.character(df_0.40$parents)
df_0.40$Primary_cutoff <- 0.40
df_0.40 <- df_0.40[,c(17,1:16)]
df_0.40_curated <- df_0.40[df_0.40$term_name %in% curated_terms,]
df_0.40_curated$F1score <- (2*(df_0.40_curated$precision*df_0.40_curated$recall)) / (df_0.40_curated$precision + df_0.40_curated$recall)
df_0.40_curated$F0.5score <- ((1 + 0.5^2) * df_0.40_curated$precision * df_0.40_curated$recall) / (0.5^2 * df_0.40_curated$precision + df_0.40_curated$recall)
df_0.40_curated$F2score <- ((1 + 2^2) * df_0.40_curated$precision * df_0.40_curated$recall) / (2^2 * df_0.40_curated$precision + df_0.40_curated$recall)

write.table(df_0.40, "1_gProfile_Results_NMF19_customToken_top10ctls_0.40_multirank.txt", sep="\t", row.names=FALSE)


# 0.45 Filter
df_0.45 <- data.frame()
wb_0.45 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.45[which(NMF_matrix_0_multi_0.45$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.45$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.45$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.45, 'gprofiler_Results_NMF19_customToken_top10ctls_0.45_multirank.xlsx', overwrite = TRUE)
df_0.45$parents <- as.character(df_0.45$parents)
df_0.45$Primary_cutoff <- 0.45
df_0.45 <- df_0.45[,c(17,1:16)]
df_0.45_curated <- df_0.45[df_0.45$term_name %in% curated_terms,]
df_0.45_curated$F1score <- (2*(df_0.45_curated$precision*df_0.45_curated$recall)) / (df_0.45_curated$precision + df_0.45_curated$recall)
df_0.45_curated$F0.5score <- ((1 + 0.5^2) * df_0.45_curated$precision * df_0.45_curated$recall) / (0.5^2 * df_0.45_curated$precision + df_0.45_curated$recall)
df_0.45_curated$F2score <- ((1 + 2^2) * df_0.45_curated$precision * df_0.45_curated$recall) / (2^2 * df_0.45_curated$precision + df_0.45_curated$recall)

write.table(df_0.45, "1_gProfile_Results_NMF19_customToken_top10ctls_0.45_multirank.txt", sep="\t", row.names=FALSE)


# 0.50 Filter
df_0.50 <- data.frame()
wb_0.50 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.50[which(NMF_matrix_0_multi_0.50$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.50$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.50$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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

saveWorkbook(wb_0.50, 'gprofiler_Results_NMF19_customToken_top10ctls_0.50_multirank.xlsx', overwrite = TRUE)
df_0.50$parents <- as.character(df_0.50$parents)
df_0.50$Primary_cutoff <- 0.50
df_0.50 <- df_0.50[,c(17,1:16)]
df_0.50_curated <- df_0.50[df_0.50$term_name %in% curated_terms,]
df_0.50_curated$F1score <- (2*(df_0.50_curated$precision*df_0.50_curated$recall)) / (df_0.50_curated$precision + df_0.50_curated$recall)
df_0.50_curated$F0.5score <- ((1 + 0.5^2) * df_0.50_curated$precision * df_0.50_curated$recall) / (0.5^2 * df_0.50_curated$precision + df_0.50_curated$recall)
df_0.50_curated$F2score <- ((1 + 2^2) * df_0.50_curated$precision * df_0.50_curated$recall) / (2^2 * df_0.50_curated$precision + df_0.50_curated$recall)

write.table(df_0.50, "1_gProfile_Results_NMF19_customToken_top10ctls_0.50_multirank.txt", sep="\t", row.names=FALSE)


# 0.55 Filter
df_0.55 <- data.frame()
wb_0.55 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.55[which(NMF_matrix_0_multi_0.55$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.55$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.55$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.55, sheet_name)
  writeData(wb_0.55, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.55 <- rbind(df2, df_0.55)
}

saveWorkbook(wb_0.55, 'gprofiler_Results_NMF19_customToken_top10ctls_0.55_multirank.xlsx', overwrite = TRUE)
df_0.55$parents <- as.character(df_0.55$parents)
df_0.55$Primary_cutoff <- 0.55
df_0.55 <- df_0.55[,c(17,1:16)]
df_0.55_curated <- df_0.55[df_0.55$term_name %in% curated_terms,]
df_0.55_curated$F1score <- (2*(df_0.55_curated$precision*df_0.55_curated$recall)) / (df_0.55_curated$precision + df_0.55_curated$recall)
df_0.55_curated$F0.5score <- ((1 + 0.5^2) * df_0.55_curated$precision * df_0.55_curated$recall) / (0.5^2 * df_0.55_curated$precision + df_0.55_curated$recall)
df_0.55_curated$F2score <- ((1 + 2^2) * df_0.55_curated$precision * df_0.55_curated$recall) / (2^2 * df_0.55_curated$precision + df_0.55_curated$recall)

write.table(df_0.55, "1_gProfile_Results_NMF19_customToken_top10ctls_0.55_multirank.txt", sep="\t", row.names=FALSE)


# 0.60 Filter
df_0.60 <- data.frame()
wb_0.60 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.60[which(NMF_matrix_0_multi_0.60$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.60$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.60$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.60, sheet_name)
  writeData(wb_0.60, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.60 <- rbind(df2, df_0.60)
}

saveWorkbook(wb_0.60, 'gprofiler_Results_NMF19_customToken_top10ctls_0.60_multirank.xlsx', overwrite = TRUE)
df_0.60$parents <- as.character(df_0.60$parents)
df_0.60$Primary_cutoff <- 0.60
df_0.60 <- df_0.60[,c(17,1:16)]
df_0.60_curated <- df_0.60[df_0.60$term_name %in% curated_terms,]
df_0.60_curated$F1score <- (2*(df_0.60_curated$precision*df_0.60_curated$recall)) / (df_0.60_curated$precision + df_0.60_curated$recall)
df_0.60_curated$F0.5score <- ((1 + 0.5^2) * df_0.60_curated$precision * df_0.60_curated$recall) / (0.5^2 * df_0.60_curated$precision + df_0.60_curated$recall)
df_0.60_curated$F2score <- ((1 + 2^2) * df_0.60_curated$precision * df_0.60_curated$recall) / (2^2 * df_0.60_curated$precision + df_0.60_curated$recall)

write.table(df_0.60, "1_gProfile_Results_NMF19_customToken_top10ctls_0.60_multirank.txt", sep="\t", row.names=FALSE)

# 0.65 Filter
df_0.65 <- data.frame()
wb_0.65 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.65[which(NMF_matrix_0_multi_0.65$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.65$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.65$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.65, sheet_name)
  writeData(wb_0.65, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.65 <- rbind(df2, df_0.65)
}

saveWorkbook(wb_0.65, 'gprofiler_Results_NMF19_customToken_top10ctls_0.65_multirank.xlsx', overwrite = TRUE)
df_0.65$parents <- as.character(df_0.65$parents)
df_0.65$Primary_cutoff <- 0.65
df_0.65 <- df_0.65[,c(17,1:16)]
df_0.65_curated <- df_0.65[df_0.65$term_name %in% curated_terms,]
df_0.65_curated$F1score <- (2*(df_0.65_curated$precision*df_0.65_curated$recall)) / (df_0.65_curated$precision + df_0.65_curated$recall)
df_0.65_curated$F0.5score <- ((1 + 0.5^2) * df_0.65_curated$precision * df_0.65_curated$recall) / (0.5^2 * df_0.65_curated$precision + df_0.65_curated$recall)
df_0.65_curated$F2score <- ((1 + 2^2) * df_0.65_curated$precision * df_0.65_curated$recall) / (2^2 * df_0.65_curated$precision + df_0.65_curated$recall)

write.table(df_0.65, "1_gProfile_Results_NMF19_customToken_top10ctls_0.65_multirank.txt", sep="\t", row.names=FALSE)

# 0.70 Filter
df_0.70 <- data.frame()
wb_0.70 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.70[which(NMF_matrix_0_multi_0.70$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.70$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.70, sheet_name)
  writeData(wb_0.70, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.70 <- rbind(df2, df_0.70)
}

saveWorkbook(wb_0.70, 'gprofiler_Results_NMF19_customToken_top10ctls_0.70_multirank.xlsx', overwrite = TRUE)
df_0.70$parents <- as.character(df_0.70$parents)
df_0.70$Primary_cutoff <- 0.70
df_0.70 <- df_0.70[,c(17,1:16)]
df_0.70_curated <- df_0.70[df_0.70$term_name %in% curated_terms,]
df_0.70_curated$F1score <- (2*(df_0.70_curated$precision*df_0.70_curated$recall)) / (df_0.70_curated$precision + df_0.70_curated$recall)
df_0.70_curated$F0.5score <- ((1 + 0.5^2) * df_0.70_curated$precision * df_0.70_curated$recall) / (0.5^2 * df_0.70_curated$precision + df_0.70_curated$recall)
df_0.70_curated$F2score <- ((1 + 2^2) * df_0.70_curated$precision * df_0.70_curated$recall) / (2^2 * df_0.70_curated$precision + df_0.70_curated$recall)

write.table(df_0.70, "1_gProfile_Results_NMF19_customToken_top10ctls_0.70_multirank.txt", sep="\t", row.names=FALSE)

# 0.75 Filter
df_0.75 <- data.frame()
wb_0.75 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.75[which(NMF_matrix_0_multi_0.75$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.75$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.75$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.75, sheet_name)
  writeData(wb_0.75, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.75 <- rbind(df2, df_0.75)
}

saveWorkbook(wb_0.75, 'gprofiler_Results_NMF19_customToken_top10ctls_0.75_multirank.xlsx', overwrite = TRUE)
df_0.75$parents <- as.character(df_0.75$parents)
df_0.75$Primary_cutoff <- 0.75
df_0.75 <- df_0.75[,c(17,1:16)]
df_0.75_curated <- df_0.75[df_0.75$term_name %in% curated_terms,]
df_0.75_curated$F1score <- (2*(df_0.75_curated$precision*df_0.75_curated$recall)) / (df_0.75_curated$precision + df_0.75_curated$recall)
df_0.75_curated$F0.5score <- ((1 + 0.5^2) * df_0.75_curated$precision * df_0.75_curated$recall) / (0.5^2 * df_0.75_curated$precision + df_0.75_curated$recall)
df_0.75_curated$F2score <- ((1 + 2^2) * df_0.75_curated$precision * df_0.75_curated$recall) / (2^2 * df_0.75_curated$precision + df_0.75_curated$recall)

write.table(df_0.75, "1_gProfile_Results_NMF19_customToken_top10ctls_0.75_multirank.txt", sep="\t", row.names=FALSE)

# 0.80 Filter
df_0.80 <- data.frame()
wb_0.80 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.80[which(NMF_matrix_0_multi_0.80$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.80$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.80$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.80, sheet_name)
  writeData(wb_0.80, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.80 <- rbind(df2, df_0.80)
}

saveWorkbook(wb_0.80, 'gprofiler_Results_NMF19_customToken_top10ctls_0.80_multirank.xlsx', overwrite = TRUE)
df_0.80$parents <- as.character(df_0.80$parents)
df_0.80$Primary_cutoff <- 0.80
df_0.80 <- df_0.80[,c(17,1:16)]
df_0.80_curated <- df_0.80[df_0.80$term_name %in% curated_terms,]
df_0.80_curated$F1score <- (2*(df_0.80_curated$precision*df_0.80_curated$recall)) / (df_0.80_curated$precision + df_0.80_curated$recall)
df_0.80_curated$F0.5score <- ((1 + 0.5^2) * df_0.80_curated$precision * df_0.80_curated$recall) / (0.5^2 * df_0.80_curated$precision + df_0.80_curated$recall)
df_0.80_curated$F2score <- ((1 + 2^2) * df_0.80_curated$precision * df_0.80_curated$recall) / (2^2 * df_0.80_curated$precision + df_0.80_curated$recall)

write.table(df_0.80, "1_gProfile_Results_NMF19_customToken_top10ctls_0.80_multirank.txt", sep="\t", row.names=FALSE)

# 0.85 Filter
df_0.85 <- data.frame()
wb_0.85 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.85[which(NMF_matrix_0_multi_0.85$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.85$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.85$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.85, sheet_name)
  writeData(wb_0.85, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.85 <- rbind(df2, df_0.85)
}

saveWorkbook(wb_0.85, 'gprofiler_Results_NMF19_customToken_top10ctls_0.85_multirank.xlsx', overwrite = TRUE)
df_0.85$parents <- as.character(df_0.85$parents)
df_0.85$Primary_cutoff <- 0.85
df_0.85 <- df_0.85[,c(17,1:16)]
df_0.85_curated <- df_0.85[df_0.85$term_name %in% curated_terms,]
df_0.85_curated$F1score <- (2*(df_0.85_curated$precision*df_0.85_curated$recall)) / (df_0.85_curated$precision + df_0.85_curated$recall)
df_0.85_curated$F0.5score <- ((1 + 0.5^2) * df_0.85_curated$precision * df_0.85_curated$recall) / (0.5^2 * df_0.85_curated$precision + df_0.85_curated$recall)
df_0.85_curated$F2score <- ((1 + 2^2) * df_0.85_curated$precision * df_0.85_curated$recall) / (2^2 * df_0.85_curated$precision + df_0.85_curated$recall)

write.table(df_0.85, "1_gProfile_Results_NMF19_customToken_top10ctls_0.85_multirank.txt", sep="\t", row.names=FALSE)

# 0.90 Filter
df_0.90 <- data.frame()
wb_0.90 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.90[which(NMF_matrix_0_multi_0.90$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.90$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.90$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.90, sheet_name)
  writeData(wb_0.90, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.90 <- rbind(df2, df_0.90)
}

saveWorkbook(wb_0.90, 'gprofiler_Results_NMF19_customToken_top10ctls_0.90_multirank.xlsx', overwrite = TRUE)
df_0.90$parents <- as.character(df_0.90$parents)
df_0.90$Primary_cutoff <- 0.90
df_0.90 <- df_0.90[,c(17,1:16)]
df_0.90_curated <- df_0.90[df_0.90$term_name %in% curated_terms,]
df_0.90_curated$F1score <- (2*(df_0.90_curated$precision*df_0.90_curated$recall)) / (df_0.90_curated$precision + df_0.90_curated$recall)
df_0.90_curated$F0.5score <- ((1 + 0.5^2) * df_0.90_curated$precision * df_0.90_curated$recall) / (0.5^2 * df_0.90_curated$precision + df_0.90_curated$recall)
df_0.90_curated$F2score <- ((1 + 2^2) * df_0.90_curated$precision * df_0.90_curated$recall) / (2^2 * df_0.90_curated$precision + df_0.90_curated$recall)

write.table(df_0.90, "1_gProfile_Results_NMF19_customToken_top10ctls_0.90_multirank.txt", sep="\t", row.names=FALSE)

# 0.95 Filter
df_0.95 <- data.frame()
wb_0.95 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_0.95[which(NMF_matrix_0_multi_0.95$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_0.95$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_0.95$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_0.95, sheet_name)
  writeData(wb_0.95, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_0.95 <- rbind(df2, df_0.95)
}

saveWorkbook(wb_0.95, 'gprofiler_Results_NMF19_customToken_top10ctls_0.95_multirank.xlsx', overwrite = TRUE)
df_0.95$parents <- as.character(df_0.95$parents)
df_0.95$Primary_cutoff <- 0.95
df_0.95 <- df_0.95[,c(17,1:16)]
df_0.95_curated <- df_0.95[df_0.95$term_name %in% curated_terms,]
df_0.95_curated$F1score <- (2*(df_0.95_curated$precision*df_0.95_curated$recall)) / (df_0.95_curated$precision + df_0.95_curated$recall)
df_0.95_curated$F0.5score <- ((1 + 0.5^2) * df_0.95_curated$precision * df_0.95_curated$recall) / (0.5^2 * df_0.95_curated$precision + df_0.95_curated$recall)
df_0.95_curated$F2score <- ((1 + 2^2) * df_0.95_curated$precision * df_0.95_curated$recall) / (2^2 * df_0.95_curated$precision + df_0.95_curated$recall)

write.table(df_0.95, "1_gProfile_Results_NMF19_customToken_top10ctls_0.95_multirank.txt", sep="\t", row.names=FALSE)

# 1.00 Filter
df_1.00 <- data.frame()
wb_1.00 = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_0_multi_1.00[which(NMF_matrix_0_multi_1.00$Primary_cluster_annotation == i | 
                                           NMF_matrix_0_multi_1.00$Secondary_cluster_annotation == i |
                                           NMF_matrix_0_multi_1.00$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
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
  addWorksheet(wb_1.00, sheet_name)
  writeData(wb_1.00, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_1.00 <- rbind(df2, df_1.00)
}

saveWorkbook(wb_1.00, 'gprofiler_Results_NMF19_customToken_top10ctls_1.00_multirank.xlsx', overwrite = TRUE)
df_1.00$parents <- as.character(df_1.00$parents)
df_1.00$Primary_cutoff <- 1.00
df_1.00 <- df_1.00[,c(17,1:16)]
df_1.00_curated <- df_1.00[df_1.00$term_name %in% curated_terms,]
df_1.00_curated$F1score <- (2*(df_1.00_curated$precision*df_1.00_curated$recall)) / (df_1.00_curated$precision + df_1.00_curated$recall)
df_1.00_curated$F0.5score <- ((1 + 0.5^2) * df_1.00_curated$precision * df_1.00_curated$recall) / (0.5^2 * df_1.00_curated$precision + df_1.00_curated$recall)
df_1.00_curated$F2score <- ((1 + 2^2) * df_1.00_curated$precision * df_1.00_curated$recall) / (2^2 * df_1.00_curated$precision + df_1.00_curated$recall)

write.table(df_1.00, "1_gProfile_Results_NMF19_customToken_top10ctls_1.00_multirank.txt", sep="\t", row.names=FALSE)


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
                             df_0.50_curated,
                             df_0.55_curated,
                             df_0.60_curated,
                             df_0.65_curated,
                             df_0.70_curated,
                             df_0.75_curated,
                             df_0.80_curated,
                             df_0.85_curated,
                             df_0.90_curated,
                             df_0.95_curated,
                             df_1.00_curated
)

write.table(df_combined_curated, "2_gProfile_Results_NMF19_customToken_top10ctls_curated_ALL_multi_rank_cutoff.txt", sep="\t", row.names=FALSE)


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

PS_0.55 <- df_0.55_curated[df_0.55_curated$term_name == "Paraspeckle curated",]
PS_0.55$negative_log10_pvalue <- -log10(PS_0.55$adjusted_p_value)
PS_0.55$Number_clusters_with_term <- nrow(PS_0.55[PS_0.55$significant == TRUE,])

PS_0.60 <- df_0.60_curated[df_0.60_curated$term_name == "Paraspeckle curated",]
PS_0.60$negative_log10_pvalue <- -log10(PS_0.60$adjusted_p_value)
PS_0.60$Number_clusters_with_term <- nrow(PS_0.60[PS_0.60$significant == TRUE,])

PS_0.65 <- df_0.65_curated[df_0.65_curated$term_name == "Paraspeckle curated",]
PS_0.65$negative_log10_pvalue <- -log10(PS_0.65$adjusted_p_value)
PS_0.65$Number_clusters_with_term <- nrow(PS_0.65[PS_0.65$significant == TRUE,])

PS_0.70 <- df_0.70_curated[df_0.70_curated$term_name == "Paraspeckle curated",]
PS_0.70$negative_log10_pvalue <- -log10(PS_0.70$adjusted_p_value)
PS_0.70$Number_clusters_with_term <- nrow(PS_0.70[PS_0.70$significant == TRUE,])

PS_0.75 <- df_0.75_curated[df_0.75_curated$term_name == "Paraspeckle curated",]
PS_0.75$negative_log10_pvalue <- -log10(PS_0.75$adjusted_p_value)
PS_0.75$Number_clusters_with_term <- nrow(PS_0.75[PS_0.75$significant == TRUE,])

PS_0.80 <- df_0.80_curated[df_0.80_curated$term_name == "Paraspeckle curated",]
PS_0.80$negative_log10_pvalue <- -log10(PS_0.80$adjusted_p_value)
PS_0.80$Number_clusters_with_term <- nrow(PS_0.80[PS_0.80$significant == TRUE,])

PS_0.85 <- df_0.85_curated[df_0.85_curated$term_name == "Paraspeckle curated",]
PS_0.85$negative_log10_pvalue <- -log10(PS_0.85$adjusted_p_value)
PS_0.85$Number_clusters_with_term <- nrow(PS_0.85[PS_0.85$significant == TRUE,])

PS_0.90 <- df_0.90_curated[df_0.90_curated$term_name == "Paraspeckle curated",]
PS_0.90$negative_log10_pvalue <- -log10(PS_0.90$adjusted_p_value)
PS_0.90$Number_clusters_with_term <- nrow(PS_0.90[PS_0.90$significant == TRUE,])

PS_0.95 <- df_0.95_curated[df_0.95_curated$term_name == "Paraspeckle curated",]
PS_0.95$negative_log10_pvalue <- -log10(PS_0.95$adjusted_p_value)
PS_0.95$Number_clusters_with_term <- nrow(PS_0.95[PS_0.95$significant == TRUE,])

PS_1.00 <- df_1.00_curated[df_1.00_curated$term_name == "Paraspeckle curated",]
PS_1.00$negative_log10_pvalue <- -log10(PS_1.00$adjusted_p_value)
PS_1.00$Number_clusters_with_term <- nrow(PS_1.00[PS_1.00$significant == TRUE,])

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
                     PS_0.50[which.max(PS_0.50$negative_log10_pvalue),],
                     PS_0.55[which.max(PS_0.55$negative_log10_pvalue),],
                     PS_0.60[which.max(PS_0.60$negative_log10_pvalue),],
                     PS_0.65[which.max(PS_0.65$negative_log10_pvalue),],
                     PS_0.70[which.max(PS_0.70$negative_log10_pvalue),],
                     PS_0.75[which.max(PS_0.75$negative_log10_pvalue),],
                     PS_0.80[which.max(PS_0.80$negative_log10_pvalue),],
                     PS_0.85[which.max(PS_0.85$negative_log10_pvalue),],
                     PS_0.90[which.max(PS_0.90$negative_log10_pvalue),],
                     PS_0.95[which.max(PS_0.95$negative_log10_pvalue),],
                     PS_1.00[which.max(PS_1.00$negative_log10_pvalue),]
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

NS_0.55 <- df_0.55_curated[df_0.55_curated$term_name == "Nuclear speckle curated",]
NS_0.55$negative_log10_pvalue <- -log10(NS_0.55$adjusted_p_value)
NS_0.55$Number_clusters_with_term <- nrow(NS_0.55[NS_0.55$significant == TRUE,])

NS_0.60 <- df_0.60_curated[df_0.60_curated$term_name == "Nuclear speckle curated",]
NS_0.60$negative_log10_pvalue <- -log10(NS_0.60$adjusted_p_value)
NS_0.60$Number_clusters_with_term <- nrow(NS_0.60[NS_0.60$significant == TRUE,])

NS_0.65 <- df_0.65_curated[df_0.65_curated$term_name == "Nuclear speckle curated",]
NS_0.65$negative_log10_pvalue <- -log10(NS_0.65$adjusted_p_value)
NS_0.65$Number_clusters_with_term <- nrow(NS_0.65[NS_0.65$significant == TRUE,])

NS_0.70 <- df_0.70_curated[df_0.70_curated$term_name == "Nuclear speckle curated",]
NS_0.70$negative_log10_pvalue <- -log10(NS_0.70$adjusted_p_value)
NS_0.70$Number_clusters_with_term <- nrow(NS_0.70[NS_0.70$significant == TRUE,])

NS_0.75 <- df_0.75_curated[df_0.75_curated$term_name == "Nuclear speckle curated",]
NS_0.75$negative_log10_pvalue <- -log10(NS_0.75$adjusted_p_value)
NS_0.75$Number_clusters_with_term <- nrow(NS_0.75[NS_0.75$significant == TRUE,])

NS_0.80 <- df_0.80_curated[df_0.80_curated$term_name == "Nuclear speckle curated",]
NS_0.80$negative_log10_pvalue <- -log10(NS_0.80$adjusted_p_value)
NS_0.80$Number_clusters_with_term <- nrow(NS_0.80[NS_0.80$significant == TRUE,])

NS_0.85 <- df_0.85_curated[df_0.85_curated$term_name == "Nuclear speckle curated",]
NS_0.85$negative_log10_pvalue <- -log10(NS_0.85$adjusted_p_value)
NS_0.85$Number_clusters_with_term <- nrow(NS_0.85[NS_0.85$significant == TRUE,])

NS_0.90 <- df_0.90_curated[df_0.90_curated$term_name == "Nuclear speckle curated",]
NS_0.90$negative_log10_pvalue <- -log10(NS_0.90$adjusted_p_value)
NS_0.90$Number_clusters_with_term <- nrow(NS_0.90[NS_0.90$significant == TRUE,])

NS_0.95 <- df_0.95_curated[df_0.95_curated$term_name == "Nuclear speckle curated",]
NS_0.95$negative_log10_pvalue <- -log10(NS_0.95$adjusted_p_value)
NS_0.95$Number_clusters_with_term <- nrow(NS_0.95[NS_0.95$significant == TRUE,])

NS_1.00 <- df_1.00_curated[df_1.00_curated$term_name == "Nuclear speckle curated",]
NS_1.00$negative_log10_pvalue <- -log10(NS_1.00$adjusted_p_value)
NS_1.00$Number_clusters_with_term <- nrow(NS_1.00[NS_1.00$significant == TRUE,])

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
                     NS_0.50[which.max(NS_0.50$negative_log10_pvalue),],
                     NS_0.55[which.max(NS_0.55$negative_log10_pvalue),],
                     NS_0.60[which.max(NS_0.60$negative_log10_pvalue),],
                     NS_0.65[which.max(NS_0.65$negative_log10_pvalue),],
                     NS_0.70[which.max(NS_0.70$negative_log10_pvalue),],
                     NS_0.75[which.max(NS_0.75$negative_log10_pvalue),],
                     NS_0.80[which.max(NS_0.80$negative_log10_pvalue),],
                     NS_0.85[which.max(NS_0.85$negative_log10_pvalue),],
                     NS_0.90[which.max(NS_0.90$negative_log10_pvalue),],
                     NS_0.95[which.max(NS_0.95$negative_log10_pvalue),],
                     NS_1.00[which.max(NS_1.00$negative_log10_pvalue),]
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

NUC_0.55 <- df_0.55_curated[df_0.55_curated$term_name == "Comprehensive nucleolus",]
NUC_0.55$negative_log10_pvalue <- -log10(NUC_0.55$adjusted_p_value)
NUC_0.55$Number_clusters_with_term <- nrow(NUC_0.55[NUC_0.55$significant == TRUE,])

NUC_0.60 <- df_0.60_curated[df_0.60_curated$term_name == "Comprehensive nucleolus",]
NUC_0.60$negative_log10_pvalue <- -log10(NUC_0.60$adjusted_p_value)
NUC_0.60$Number_clusters_with_term <- nrow(NUC_0.60[NUC_0.60$significant == TRUE,])

NUC_0.65 <- df_0.65_curated[df_0.65_curated$term_name == "Comprehensive nucleolus",]
NUC_0.65$negative_log10_pvalue <- -log10(NUC_0.65$adjusted_p_value)
NUC_0.65$Number_clusters_with_term <- nrow(NUC_0.65[NUC_0.65$significant == TRUE,])

NUC_0.70 <- df_0.70_curated[df_0.70_curated$term_name == "Comprehensive nucleolus",]
NUC_0.70$negative_log10_pvalue <- -log10(NUC_0.70$adjusted_p_value)
NUC_0.70$Number_clusters_with_term <- nrow(NUC_0.70[NUC_0.70$significant == TRUE,])

NUC_0.75 <- df_0.75_curated[df_0.75_curated$term_name == "Comprehensive nucleolus",]
NUC_0.75$negative_log10_pvalue <- -log10(NUC_0.75$adjusted_p_value)
NUC_0.75$Number_clusters_with_term <- nrow(NUC_0.75[NUC_0.75$significant == TRUE,])

NUC_0.80 <- df_0.80_curated[df_0.80_curated$term_name == "Comprehensive nucleolus",]
NUC_0.80$negative_log10_pvalue <- -log10(NUC_0.80$adjusted_p_value)
NUC_0.80$Number_clusters_with_term <- nrow(NUC_0.80[NUC_0.80$significant == TRUE,])

NUC_0.85 <- df_0.85_curated[df_0.85_curated$term_name == "Comprehensive nucleolus",]
NUC_0.85$negative_log10_pvalue <- -log10(NUC_0.85$adjusted_p_value)
NUC_0.85$Number_clusters_with_term <- nrow(NUC_0.85[NUC_0.85$significant == TRUE,])

NUC_0.90 <- df_0.90_curated[df_0.90_curated$term_name == "Comprehensive nucleolus",]
NUC_0.90$negative_log10_pvalue <- -log10(NUC_0.90$adjusted_p_value)
NUC_0.90$Number_clusters_with_term <- nrow(NUC_0.90[NUC_0.90$significant == TRUE,])

NUC_0.95 <- df_0.95_curated[df_0.95_curated$term_name == "Comprehensive nucleolus",]
NUC_0.95$negative_log10_pvalue <- -log10(NUC_0.95$adjusted_p_value)
NUC_0.95$Number_clusters_with_term <- nrow(NUC_0.95[NUC_0.95$significant == TRUE,])

NUC_1.00 <- df_1.00_curated[df_1.00_curated$term_name == "Comprehensive nucleolus",]
NUC_1.00$negative_log10_pvalue <- -log10(NUC_1.00$adjusted_p_value)
NUC_1.00$Number_clusters_with_term <- nrow(NUC_1.00[NUC_1.00$significant == TRUE,])

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
                     NUC_0.50[which.max(NUC_0.50$negative_log10_pvalue),],
                     NUC_0.55[which.max(NUC_0.55$negative_log10_pvalue),],
                     NUC_0.60[which.max(NUC_0.60$negative_log10_pvalue),],
                     NUC_0.65[which.max(NUC_0.65$negative_log10_pvalue),],
                     NUC_0.70[which.max(NUC_0.70$negative_log10_pvalue),],
                     NUC_0.75[which.max(NUC_0.75$negative_log10_pvalue),],
                     NUC_0.80[which.max(NUC_0.80$negative_log10_pvalue),],
                     NUC_0.85[which.max(NUC_0.85$negative_log10_pvalue),],
                     NUC_0.90[which.max(NUC_0.90$negative_log10_pvalue),],
                     NUC_0.95[which.max(NUC_0.95$negative_log10_pvalue),],
                     NUC_1.00[which.max(NUC_1.00$negative_log10_pvalue),]
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

CB_0.55 <- df_0.55_curated[df_0.55_curated$term_name == "Cajal body curated",]
CB_0.55$negative_log10_pvalue <- -log10(CB_0.55$adjusted_p_value)
CB_0.55$Number_clusters_with_term <- nrow(CB_0.55[CB_0.55$significant == TRUE,])

CB_0.60 <- df_0.60_curated[df_0.60_curated$term_name == "Cajal body curated",]
CB_0.60$negative_log10_pvalue <- -log10(CB_0.60$adjusted_p_value)
CB_0.60$Number_clusters_with_term <- nrow(CB_0.60[CB_0.60$significant == TRUE,])

CB_0.65 <- df_0.65_curated[df_0.65_curated$term_name == "Cajal body curated",]
CB_0.65$negative_log10_pvalue <- -log10(CB_0.65$adjusted_p_value)
CB_0.65$Number_clusters_with_term <- nrow(CB_0.65[CB_0.65$significant == TRUE,])

CB_0.70 <- df_0.70_curated[df_0.70_curated$term_name == "Cajal body curated",]
CB_0.70$negative_log10_pvalue <- -log10(CB_0.70$adjusted_p_value)
CB_0.70$Number_clusters_with_term <- nrow(CB_0.70[CB_0.70$significant == TRUE,])

CB_0.75 <- df_0.75_curated[df_0.75_curated$term_name == "Cajal body curated",]
CB_0.75$negative_log10_pvalue <- -log10(CB_0.75$adjusted_p_value)
CB_0.75$Number_clusters_with_term <- nrow(CB_0.75[CB_0.75$significant == TRUE,])

CB_0.80 <- df_0.80_curated[df_0.80_curated$term_name == "Cajal body curated",]
CB_0.80$negative_log10_pvalue <- -log10(CB_0.80$adjusted_p_value)
CB_0.80$Number_clusters_with_term <- nrow(CB_0.80[CB_0.80$significant == TRUE,])

CB_0.85 <- df_0.85_curated[df_0.85_curated$term_name == "Cajal body curated",]
CB_0.85$negative_log10_pvalue <- -log10(CB_0.85$adjusted_p_value)
CB_0.85$Number_clusters_with_term <- nrow(CB_0.85[CB_0.85$significant == TRUE,])

CB_0.90 <- df_0.90_curated[df_0.90_curated$term_name == "Cajal body curated",]
CB_0.90$negative_log10_pvalue <- -log10(CB_0.90$adjusted_p_value)
CB_0.90$Number_clusters_with_term <- nrow(CB_0.90[CB_0.90$significant == TRUE,])

CB_0.95 <- df_0.95_curated[df_0.95_curated$term_name == "Cajal body curated",]
CB_0.95$negative_log10_pvalue <- -log10(CB_0.95$adjusted_p_value)
CB_0.95$Number_clusters_with_term <- nrow(CB_0.95[CB_0.95$significant == TRUE,])

CB_1.00 <- df_1.00_curated[df_1.00_curated$term_name == "Cajal body curated",]
CB_1.00$negative_log10_pvalue <- -log10(CB_1.00$adjusted_p_value)
CB_1.00$Number_clusters_with_term <- nrow(CB_1.00[CB_1.00$significant == TRUE,])

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
                     CB_0.50[which.max(CB_0.50$negative_log10_pvalue),],
                     CB_0.55[which.max(CB_0.55$negative_log10_pvalue),],
                     CB_0.60[which.max(CB_0.60$negative_log10_pvalue),],
                     CB_0.65[which.max(CB_0.65$negative_log10_pvalue),],
                     CB_0.70[which.max(CB_0.70$negative_log10_pvalue),],
                     CB_0.75[which.max(CB_0.75$negative_log10_pvalue),],
                     CB_0.80[which.max(CB_0.80$negative_log10_pvalue),],
                     CB_0.85[which.max(CB_0.85$negative_log10_pvalue),],
                     CB_0.90[which.max(CB_0.90$negative_log10_pvalue),],
                     CB_0.95[which.max(CB_0.95$negative_log10_pvalue),],
                     CB_1.00[which.max(CB_1.00$negative_log10_pvalue),]
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

PML_0.55 <- df_0.55_curated[df_0.55_curated$term_name == "PML body curated",]
PML_0.55$negative_log10_pvalue <- -log10(PML_0.55$adjusted_p_value)
PML_0.55$Number_clusters_with_term <- nrow(PML_0.55[PML_0.55$significant == TRUE,])

PML_0.60 <- df_0.60_curated[df_0.60_curated$term_name == "PML body curated",]
PML_0.60$negative_log10_pvalue <- -log10(PML_0.60$adjusted_p_value)
PML_0.60$Number_clusters_with_term <- nrow(PML_0.60[PML_0.60$significant == TRUE,])

PML_0.65 <- df_0.65_curated[df_0.65_curated$term_name == "PML body curated",]
PML_0.65$negative_log10_pvalue <- -log10(PML_0.65$adjusted_p_value)
PML_0.65$Number_clusters_with_term <- nrow(PML_0.65[PML_0.65$significant == TRUE,])

PML_0.70 <- df_0.70_curated[df_0.70_curated$term_name == "PML body curated",]
PML_0.70$negative_log10_pvalue <- -log10(PML_0.70$adjusted_p_value)
PML_0.70$Number_clusters_with_term <- nrow(PML_0.70[PML_0.70$significant == TRUE,])

PML_0.75 <- df_0.75_curated[df_0.75_curated$term_name == "PML body curated",]
PML_0.75$negative_log10_pvalue <- -log10(PML_0.75$adjusted_p_value)
PML_0.75$Number_clusters_with_term <- nrow(PML_0.75[PML_0.75$significant == TRUE,])

PML_0.80 <- df_0.80_curated[df_0.80_curated$term_name == "PML body curated",]
PML_0.80$negative_log10_pvalue <- -log10(PML_0.80$adjusted_p_value)
PML_0.80$Number_clusters_with_term <- nrow(PML_0.80[PML_0.80$significant == TRUE,])

PML_0.85 <- df_0.85_curated[df_0.85_curated$term_name == "PML body curated",]
PML_0.85$negative_log10_pvalue <- -log10(PML_0.85$adjusted_p_value)
PML_0.85$Number_clusters_with_term <- nrow(PML_0.85[PML_0.85$significant == TRUE,])

PML_0.90 <- df_0.90_curated[df_0.90_curated$term_name == "PML body curated",]
PML_0.90$negative_log10_pvalue <- -log10(PML_0.90$adjusted_p_value)
PML_0.90$Number_clusters_with_term <- nrow(PML_0.90[PML_0.90$significant == TRUE,])

PML_0.95 <- df_0.95_curated[df_0.95_curated$term_name == "PML body curated",]
PML_0.95$negative_log10_pvalue <- -log10(PML_0.95$adjusted_p_value)
PML_0.95$Number_clusters_with_term <- nrow(PML_0.95[PML_0.95$significant == TRUE,])

PML_1.00 <- df_1.00_curated[df_1.00_curated$term_name == "PML body curated",]
PML_1.00$negative_log10_pvalue <- -log10(PML_1.00$adjusted_p_value)
PML_1.00$Number_clusters_with_term <- nrow(PML_1.00[PML_1.00$significant == TRUE,])

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
                     PML_0.50[which.max(PML_0.50$negative_log10_pvalue),],
                     PML_0.55[which.max(PML_0.55$negative_log10_pvalue),],
                     PML_0.60[which.max(PML_0.60$negative_log10_pvalue),],
                     PML_0.65[which.max(PML_0.65$negative_log10_pvalue),],
                     PML_0.70[which.max(PML_0.70$negative_log10_pvalue),],
                     PML_0.75[which.max(PML_0.75$negative_log10_pvalue),],
                     PML_0.80[which.max(PML_0.80$negative_log10_pvalue),],
                     PML_0.85[which.max(PML_0.85$negative_log10_pvalue),],
                     PML_0.90[which.max(PML_0.90$negative_log10_pvalue),],
                     PML_0.95[which.max(PML_0.95$negative_log10_pvalue),],
                     PML_1.00[which.max(PML_1.00$negative_log10_pvalue),]
)


# Combine all top rank data into one file
df_combined_curated_top <- rbind(top_rank_PS,
                                 top_rank_NS,
                                 top_rank_NUC,
                                 top_rank_CB,
                                 top_rank_PML
)

write.table(df_combined_curated_top, "3_gProfile_Results_NMF19_curated_terms_combined_top10ctls_multi_rank_cutoff.txt", sep="\t", row.names=FALSE)

#--------
# Plot PS data
library(ggplot2)
library(ggrepel)

#Negative log10 p-values
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_pValue_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=negative_log10_pvalue, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF multi-score cutoff testing (gProfiler:Negative log10 pvalues) Top 10 ctl subtraction")+
  xlab("NMF multi-score cutoff")+
  ylab("-log10(pvalue) (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,150))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Recall
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_Recall_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=recall, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Recall) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
  ylab("recall (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,0.5))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Precision
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_Precision_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=precision, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Precision) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
  ylab("precision (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,1.0))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Query size
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_QuerySize_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=query_size, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Query Size) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
  ylab("Query size (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,350))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F1scores
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_F1score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F1score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:F1-score)")+
  xlab("NMF multi-rank cutoff")+
  ylab("F1-score (2*Prec*Recall)/(Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F0.5 scores
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_F0.5score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F0.5score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:F0.5-score)")+
  xlab("NMF multi-rank cutoff")+
  ylab("F0.5-score ((1+0.5^2)*Prec*Recall)/(0.5^2*Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F2scores
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_F2score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Primary_cutoff, y=F2score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:F2-score)")+
  xlab("NMF multi-rank cutoff")+
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
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_avgF1score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F1scores, aes(x=Primary_cutoff, y=average_F1score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
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
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_avgF0.5score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F0.5scores, aes(x=Primary_cutoff, y=average_F0.5score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
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
pdf("Rplot_NMF_19ranks_multi_score_cutoff_testing_gProfiler_TopPvalueRank_avgF2score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F2scores, aes(x=Primary_cutoff, y=average_F2score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF multi-rank cutoff testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("NMF multi-rank cutoff")+
  ylab("Average F2score F2-score ((1+2^2)*Prec*Recall)/(2^2*Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

Fscores_avg_combined <- cbind(F0.5scores,F1scores,F2scores)
Fscores_avg_combined <- Fscores_avg_combined[,c(1,2,4,6)]
write.table(Fscores_avg_combined, "4_gProfile_Results_NMF19_curated_terms_top10ctls_NMFcutoff_multirank_toprank_avgFscores.txt", sep="\t", row.names=FALSE)
