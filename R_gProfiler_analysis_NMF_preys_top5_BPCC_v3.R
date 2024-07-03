# Run gProfiler for Boris baits (separated by groups)

#install.packages("gprofiler2")

## Part1: Read in the NMF file and filtering according to criteria 
## Part2: Run gProfiler for each rank (primary, or multi-rank) and output results as a file to input into Prohits-viz (topn values)

#---------
## Part1: Read in the NMF file and filtering according to criteria 
#Read in NMF file and re-annotate
NMF_matrix <- read.delim("7013_cleaned_v2_NMF_19_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix <- NMF_matrix[,-(ncol(NMF_matrix))]

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

#write.table(NMF_matrix_col, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks.csv", sep=",", row.names=FALSE)

#-----------
# Re-annotate primary and multi-annotation columns filtered by different percent of max in column filters
NMF_matrix_col_0.05_multi_0.70 <- NMF_matrix_col
NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05, NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Secondary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Multi_loc_score_secondary >= 0.70, NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation <- ifelse(NMF_matrix_col_0.05_multi_0.70$Primary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Tertiary_ratio >= 0.05 & NMF_matrix_col_0.05_multi_0.70$Multi_loc_score_tertiary >= 0.70, NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation, NA)
NMF_matrix_col_0.05_multi_0.70$Multi_rank <- ifelse(is.na(NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation) & is.na(NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation), "single", "multi")
write.table(NMF_matrix_col_0.05_multi_0.70, "0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", row.names=FALSE)


#--------
## Part2: Run gProfiler for each rank (primary, or multi-rank) and output results as a file to input into Prohits-viz (topn values)
library("gprofiler2")
library(openxlsx)

ranks <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18")
#ranks <- c("X0","X1","X2")

df_NMF_all <- data.frame()
df_NMF_topn <- data.frame()
wb_NMF = createWorkbook()

for (i in ranks) {
  query <- NMF_matrix_col_0.05_multi_0.70[which(NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation == i | 
                                                  NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation == i |
                                                  NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  if (is.null(gprofiler$result)) {
    next  # Skip this iteration if there are no results
  }
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_NMF, sheet_name)
  writeData(wb_NMF, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NMF_all <- rbind(df2, df_NMF_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_NMF_topn <- rbind(df_topn, df_NMF_topn)
}

saveWorkbook(wb_NMF, 'gprofiler_Results_NMF19_top10ctls_0.05_percentfilter_0.70_multi_BPCC.xlsx', overwrite = TRUE)
df_NMF_all$parents <- as.character(df_NMF_all$parents)
df_NMF_topn$parents <- as.character(df_NMF_topn$parents)
write.table(df_NMF_all, "2_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_NMF_topn, "2_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Outputting all results (for viz)
df_NMF_noFilter_all <- data.frame()

for (i in ranks) {
  query <- NMF_matrix_col_0.05_multi_0.70[which(NMF_matrix_col_0.05_multi_0.70$Primary_cluster_annotation == i | 
                                                  NMF_matrix_col_0.05_multi_0.70$Secondary_cluster_annotation == i |
                                                  NMF_matrix_col_0.05_multi_0.70$Tertiary_cluster_annotation == i),] 
  if (nrow(query) == 0) {
    next  # Skip this iteration if the query is empty
  }
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  if (is.null(gprofiler$result)) {
    next  # Skip this iteration if there are no results
  }
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_NMF_noFilter_all <- rbind(df_notsign, df_NMF_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_NMF_noFilter_topn <- df_NMF_noFilter_all[df_NMF_noFilter_all$term_id %in% df_NMF_topn$term_id,]
df_NMF_noFilter_topn$parents <- as.character(df_NMF_noFilter_topn$parents)
df_NMF_noFilter_topn$signif_score <- ifelse(df_NMF_noFilter_topn$significant == TRUE, 1, 0)
df_NMF_noFilter_topn$negative_log10_padj <- -log10(df_NMF_noFilter_topn$adjusted_p_value)

write.table(df_NMF_noFilter_topn, "3_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_NMF_custom_all <- data.frame()
df_NMF_custom_curated <- data.frame()
wb_custom_NMF = createWorkbook()

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
                    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
  df_custom <- data.frame(gprofiler$result)
  colnames(df_custom)[3] <- "adjusted_p_value"
  df_custom$query <- i
  df_custom <- df_custom[order(df_custom$adjusted_p_value), ]
  df_custom$parents <- lapply(df_custom$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_custom_NMF, sheet_name)
  writeData(wb_custom_NMF, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_NMF_custom_all <- rbind(df2_custom, df_NMF_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NMF_custom_curated <- rbind(df_curated, df_NMF_custom_curated)
  
}

saveWorkbook(wb_custom_NMF, 'gprofiler_Results_NMF19_top10ctls_0.05_percentfilter_0.70_multi_BPCC_customtoken.xlsx', overwrite = TRUE)
df_NMF_custom_all$parents <- as.character(df_NMF_custom_all$parents)
df_NMF_custom_all <- df_NMF_custom_all[df_NMF_custom_all$significant == TRUE,]

df_NMF_custom_curated$parents <- as.character(df_NMF_custom_curated$parents)
df_NMF_custom_curated$signif_score <- ifelse(df_NMF_custom_curated$significant == TRUE, 1, 0)
df_NMF_custom_curated$negative_log10_padj <- -log10(df_NMF_custom_curated$adjusted_p_value)

write.table(df_NMF_custom_all, "3_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_customtoken_.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_NMF_custom_curated, "3_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_customtoken_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_NMF_combined <- rbind(df_NMF_noFilter_topn, df_NMF_custom_curated)
write.table(df_NMF_combined, "4_input_Prohits_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
library(stringr)
NMF_parents <- data.frame(str_split_fixed(df_NMF_combined$parents, "\\|", 4))
df_NMF_combined_temp <- cbind(df_NMF_combined,NMF_parents)
rows_to_remove <- apply(df_NMF_combined_temp[, 19:ncol(df_NMF_combined_temp)], 1, function(row) any(row %in% df_NMF_combined_temp$term_id))
df_NMF_combined_temp <- df_NMF_combined_temp[!rows_to_remove, ]
df_NMF_combined_temp_2 <- df_NMF_combined_temp[,-c(19:ncol(df_NMF_combined_temp))]
test <- df_NMF_combined[rows_to_remove, ]

write.table(df_NMF_combined_temp_2, "5_input_Prohits_NMF19_0.05_percentfilter_0.70_multi_gProfileroutput_BPCC_wNonSign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)

