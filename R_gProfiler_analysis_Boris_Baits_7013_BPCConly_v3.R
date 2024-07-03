# Run gProfiler for Boris baits (separated by groups)

#install.packages("gprofiler2")

## Part1: Read in SAINT file and output gene lists of signif (BFDR < 0.01) hits
## Part2: Run gProfiler on each bait and combine into excel file to input into dotplot tool

#--------
## Part1: Read in SAINT file and output gene lists of signif (BFDR < 0.01) hits 
# Input SAINT file
SAINT_input <- read.delim("SAINT_7013_cleaned_v2.txt", sep="\t", header=TRUE)
bait_input <- read.delim("Supplementary Table 2 - BirA bait table - v2_cleaned.csv", sep=",", header=TRUE)

# Create column with gene name 
baitGene <- data.frame(do.call('rbind', strsplit(as.character(SAINT_input$Bait),'_',fixed=TRUE)))
SAINT_input <- cbind(SAINT_input,baitGene)
SAINT_input <- SAINT_input[,c(1,23,2:22)]
colnames(SAINT_input)[2] <- "BaitGene"

# Output file to help edit bait_input file
bait_uniq <- unique(SAINT_input[,c(1,2)])
write.table(bait_uniq, "0.5_Boris_baits_unique_list.txt", sep="\t", row.names=FALSE)

# Test
#test <- SAINT_input[SAINT_input$Bait == "RNPS1_CBF" & SAINT_input$BFDR <= 0.01,]
#write.table(test, "0.5_gProfiler_bait_RNPS1_sign_test.txt", sep="\t", row.names=FALSE)

## Part2: Run gProfiler on each bait and combine into excel file to input into dotplot tool
library("gprofiler2")
library(openxlsx)

#-------
# NS baits
baits_NS <- bait_input[bait_input$Primary.annotation == "Nuclear speckle",]
namelist_NS <- baits_NS$Bait.name.in.SAINT.results.6099..supp.table.3.
#namelist_NS <- c("MFAP1_CBF","RNPS1_CBF","SRSF1_CBF")
write.table(baits_NS, "1_Boris_baits_NS_list.txt", sep="\t", row.names=FALSE)

df_NS_all <- data.frame()
df_NS_topn <- data.frame()
wb_NS = createWorkbook()

for (i in namelist_NS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_NS, sheet_name)
  writeData(wb_NS, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NS_all <- rbind(df2, df_NS_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_NS_topn <- rbind(df_topn, df_NS_topn)
}

saveWorkbook(wb_NS, 'gprofiler_NS_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_NS_all$parents <- as.character(df_NS_all$parents)
df_NS_topn$parents <- as.character(df_NS_topn$parents)
write.table(df_NS_all, "2_Boris_baits_NS_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_NS_topn, "2_Boris_baits_NS_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_NS_noFilter_all <- data.frame()

for (i in namelist_NS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_NS_noFilter_all <- rbind(df_notsign, df_NS_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_NS_noFilter_topn <- df_NS_noFilter_all[df_NS_noFilter_all$term_id %in% df_NS_topn$term_id,]
df_NS_noFilter_topn$parents <- as.character(df_NS_noFilter_topn$parents)
df_NS_noFilter_topn$signif_score <- ifelse(df_NS_noFilter_topn$significant == TRUE, 1, 0)
df_NS_noFilter_topn$negative_log10_padj <- -log10(df_NS_noFilter_topn$adjusted_p_value)

write.table(df_NS_noFilter_topn, "3_Boris_baits_NS_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_NS_custom_all <- data.frame()
df_NS_custom_curated <- data.frame()
wb_custom_NS = createWorkbook()

for (i in namelist_NS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_NS, sheet_name)
  writeData(wb_custom_NS, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_NS_custom_all <- rbind(df2_custom, df_NS_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NS_custom_curated <- rbind(df_curated, df_NS_custom_curated)
  
}

saveWorkbook(wb_custom_NS, 'gprofiler_NS_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_NS_custom_all$parents <- as.character(df_NS_custom_all$parents)
df_NS_custom_all <- df_NS_custom_all[df_NS_custom_all$significant == TRUE,]
df_NS_custom_curated$parents <- as.character(df_NS_custom_curated$parents)
df_NS_custom_curated$signif_score <- ifelse(df_NS_custom_curated$significant == TRUE, 1, 0)
df_NS_custom_curated$negative_log10_padj <- -log10(df_NS_custom_curated$adjusted_p_value)

write.table(df_NS_custom_all, "2_Boris_baits_NS_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_NS_custom_curated, "3_Boris_baits_NS_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_NS_combined <- rbind(df_NS_noFilter_topn, df_NS_custom_curated)
write.table(df_NS_combined, "4_input_Prohits_Boris_baits_NS_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
library(stringr)
NS_parents <- data.frame(str_split_fixed(df_NS_combined$parents, "\\|", 4))
df_NS_combined_temp <- cbind(df_NS_combined,NS_parents)
rows_to_remove <- apply(df_NS_combined_temp[, 19:ncol(df_NS_combined_temp)], 1, function(row) any(row %in% df_NS_combined_temp$term_id))
df_NS_combined_temp <- df_NS_combined_temp[!rows_to_remove, ]
df_NS_combined_temp_2 <- df_NS_combined_temp[,-c(19:ncol(df_NS_combined_temp))]
#test <- df_NS_combined[rows_to_remove, ]

write.table(df_NS_combined_temp_2, "5_input_Prohits_Boris_baits_NS_gProfileroutput_BPCC_curated_wNonSign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)


#-------
# PS baits
baits_PS <- bait_input[bait_input$Primary.annotation == "Paraspeckle",]
namelist_PS <- baits_PS$Bait.name.in.SAINT.results.6099..supp.table.3.

write.table(baits_PS, "1_Boris_baits_PS_list.txt", sep="\t", row.names=FALSE)

df_PS_all <- data.frame()
df_PS_topn <- data.frame()
wb_PS = createWorkbook()

for (i in namelist_PS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_PS, sheet_name)
  writeData(wb_PS, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_PS_all <- rbind(df2, df_PS_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_PS_topn <- rbind(df_topn, df_PS_topn)
}

saveWorkbook(wb_PS, 'gprofiler_PS_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_PS_all$parents <- as.character(df_PS_all$parents)
df_PS_topn$parents <- as.character(df_PS_topn$parents)
write.table(df_PS_all, "2_Boris_baits_PS_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_PS_topn, "2_Boris_baits_PS_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_PS_noFilter_all <- data.frame()

for (i in namelist_PS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_PS_noFilter_all <- rbind(df_notsign, df_PS_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_PS_noFilter_topn <- df_PS_noFilter_all[df_PS_noFilter_all$term_id %in% df_PS_topn$term_id,]
df_PS_noFilter_topn$parents <- as.character(df_PS_noFilter_topn$parents)
df_PS_noFilter_topn$signif_score <- ifelse(df_PS_noFilter_topn$significant == TRUE, 1, 0)
df_PS_noFilter_topn$negative_log10_padj <- -log10(df_PS_noFilter_topn$adjusted_p_value)

write.table(df_PS_noFilter_topn, "3_Boris_baits_PS_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_PS_custom_all <- data.frame()
df_PS_custom_curated <- data.frame()
wb_custom_PS = createWorkbook()

for (i in namelist_PS) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_PS, sheet_name)
  writeData(wb_custom_PS, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_PS_custom_all <- rbind(df2_custom, df_PS_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_PS_custom_curated <- rbind(df_curated, df_PS_custom_curated)
  
}

saveWorkbook(wb_custom_PS, 'gprofiler_PS_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_PS_custom_all$parents <- as.character(df_PS_custom_all$parents)
df_PS_custom_all <- df_PS_custom_all[df_PS_custom_all$significant == TRUE,]
df_PS_custom_curated$parents <- as.character(df_PS_custom_curated$parents)
df_PS_custom_curated$signif_score <- ifelse(df_PS_custom_curated$significant == TRUE, 1, 0)
df_PS_custom_curated$negative_log10_padj <- -log10(df_PS_custom_curated$adjusted_p_value)

write.table(df_PS_custom_all, "2_Boris_baits_PS_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_PS_custom_curated, "3_Boris_baits_PS_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_PS_combined <- rbind(df_PS_noFilter_topn, df_PS_custom_curated)
write.table(df_PS_combined, "4_input_Prohits_Boris_baits_PS_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
PS_parents <- data.frame(str_split_fixed(df_PS_combined$parents, "\\|", 4))
df_PS_combined_temp <- cbind(df_PS_combined,PS_parents)
rows_to_remove <- apply(df_PS_combined_temp[, 19:ncol(df_PS_combined_temp)], 1, function(row) any(row %in% df_PS_combined_temp$term_id))
df_PS_combined_temp <- df_PS_combined_temp[!rows_to_remove, ]
df_PS_combined_temp_2 <- df_PS_combined_temp[,-c(19:ncol(df_PS_combined_temp))]

write.table(df_PS_combined_temp_2, "5_input_Prohits_Boris_baits_PS_gProfileroutput_BPCC_curated_wNoPSign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)


#-------
# CB baits
baits_CB <- bait_input[bait_input$Primary.annotation == "Cajal body",]
namelist_CB <- baits_CB$Bait.name.in.SAINT.results.6099..supp.table.3.
write.table(baits_CB, "1_Boris_baits_CB_list.txt", sep="\t", row.names=FALSE)

df_CB_all <- data.frame()
df_CB_topn <- data.frame()
wb_CB = createWorkbook()

for (i in namelist_CB) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_CB, sheet_name)
  writeData(wb_CB, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_CB_all <- rbind(df2, df_CB_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_CB_topn <- rbind(df_topn, df_CB_topn)
}

saveWorkbook(wb_CB, 'gprofiler_CB_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_CB_all$parents <- as.character(df_CB_all$parents)
df_CB_topn$parents <- as.character(df_CB_topn$parents)
write.table(df_CB_all, "2_Boris_baits_CB_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_CB_topn, "2_Boris_baits_CB_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_CB_noFilter_all <- data.frame()

for (i in namelist_CB) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_CB_noFilter_all <- rbind(df_notsign, df_CB_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_CB_noFilter_topn <- df_CB_noFilter_all[df_CB_noFilter_all$term_id %in% df_CB_topn$term_id,]
df_CB_noFilter_topn$parents <- as.character(df_CB_noFilter_topn$parents)
df_CB_noFilter_topn$signif_score <- ifelse(df_CB_noFilter_topn$significant == TRUE, 1, 0)
df_CB_noFilter_topn$negative_log10_padj <- -log10(df_CB_noFilter_topn$adjusted_p_value)

write.table(df_CB_noFilter_topn, "3_Boris_baits_CB_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_CB_custom_all <- data.frame()
df_CB_custom_curated <- data.frame()
wb_custom_CB = createWorkbook()

for (i in namelist_CB) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_CB, sheet_name)
  writeData(wb_custom_CB, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_CB_custom_all <- rbind(df2_custom, df_CB_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_CB_custom_curated <- rbind(df_curated, df_CB_custom_curated)
  
}

saveWorkbook(wb_custom_CB, 'gprofiler_CB_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_CB_custom_all$parents <- as.character(df_CB_custom_all$parents)
df_CB_custom_all <- df_CB_custom_all[df_CB_custom_all$significant == TRUE,]
df_CB_custom_curated$parents <- as.character(df_CB_custom_curated$parents)
df_CB_custom_curated$signif_score <- ifelse(df_CB_custom_curated$significant == TRUE, 1, 0)
df_CB_custom_curated$negative_log10_padj <- -log10(df_CB_custom_curated$adjusted_p_value)

write.table(df_CB_custom_all, "2_Boris_baits_CB_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_CB_custom_curated, "3_Boris_baits_CB_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_CB_combined <- rbind(df_CB_noFilter_topn, df_CB_custom_curated)
write.table(df_CB_combined, "4_input_Prohits_Boris_baits_CB_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
CB_parents <- data.frame(str_split_fixed(df_CB_combined$parents, "\\|", 4))
df_CB_combined_temp <- cbind(df_CB_combined,CB_parents)
rows_to_remove <- apply(df_CB_combined_temp[, 19:ncol(df_CB_combined_temp)], 1, function(row) any(row %in% df_CB_combined_temp$term_id))
df_CB_combined_temp <- df_CB_combined_temp[!rows_to_remove, ]
df_CB_combined_temp_2 <- df_CB_combined_temp[,-c(19:ncol(df_CB_combined_temp))]

write.table(df_CB_combined_temp_2, "5_input_Prohits_Boris_baits_CB_gProfileroutput_BPCC_curated_wNoCBign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)


#-------
# Nucleolus baits
baits_NUC <- bait_input[bait_input$Primary.annotation == "Nucleolus",]
namelist_NUC <- baits_NUC$Bait.name.in.SAINT.results.6099..supp.table.3.
write.table(baits_NUC, "1_Boris_baits_Nucleolus_list.txt", sep="\t", row.names=FALSE)

df_NUC_all <- data.frame()
df_NUC_topn <- data.frame()
wb_NUC = createWorkbook()

for (i in namelist_NUC) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_NUC, sheet_name)
  writeData(wb_NUC, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NUC_all <- rbind(df2, df_NUC_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_NUC_topn <- rbind(df_topn, df_NUC_topn)
}

saveWorkbook(wb_NUC, 'gprofiler_NUC_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_NUC_all$parents <- as.character(df_NUC_all$parents)
df_NUC_topn$parents <- as.character(df_NUC_topn$parents)
write.table(df_NUC_all, "2_Boris_baits_NUC_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_NUC_topn, "2_Boris_baits_NUC_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_NUC_noFilter_all <- data.frame()

for (i in namelist_NUC) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_NUC_noFilter_all <- rbind(df_notsign, df_NUC_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_NUC_noFilter_topn <- df_NUC_noFilter_all[df_NUC_noFilter_all$term_id %in% df_NUC_topn$term_id,]
df_NUC_noFilter_topn$parents <- as.character(df_NUC_noFilter_topn$parents)
df_NUC_noFilter_topn$signif_score <- ifelse(df_NUC_noFilter_topn$significant == TRUE, 1, 0)
df_NUC_noFilter_topn$negative_log10_padj <- -log10(df_NUC_noFilter_topn$adjusted_p_value)

write.table(df_NUC_noFilter_topn, "3_Boris_baits_NUC_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_NUC_custom_all <- data.frame()
df_NUC_custom_curated <- data.frame()
wb_custom_NUC = createWorkbook()

for (i in namelist_NUC) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_NUC, sheet_name)
  writeData(wb_custom_NUC, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_NUC_custom_all <- rbind(df2_custom, df_NUC_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_NUC_custom_curated <- rbind(df_curated, df_NUC_custom_curated)
  
}

saveWorkbook(wb_custom_NUC, 'gprofiler_NUC_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_NUC_custom_all$parents <- as.character(df_NUC_custom_all$parents)
df_NUC_custom_all <- df_NUC_custom_all[df_NUC_custom_all$significant == TRUE,]
df_NUC_custom_curated$parents <- as.character(df_NUC_custom_curated$parents)
df_NUC_custom_curated$signif_score <- ifelse(df_NUC_custom_curated$significant == TRUE, 1, 0)
df_NUC_custom_curated$negative_log10_padj <- -log10(df_NUC_custom_curated$adjusted_p_value)

write.table(df_NUC_custom_all, "2_Boris_baits_NUC_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_NUC_custom_curated, "3_Boris_baits_NUC_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_NUC_combined <- rbind(df_NUC_noFilter_topn, df_NUC_custom_curated)
write.table(df_NUC_combined, "4_input_Prohits_Boris_baits_NUC_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
NUC_parents <- data.frame(str_split_fixed(df_NUC_combined$parents, "\\|", 4))
df_NUC_combined_temp <- cbind(df_NUC_combined,NUC_parents)
rows_to_remove <- apply(df_NUC_combined_temp[, 19:ncol(df_NUC_combined_temp)], 1, function(row) any(row %in% df_NUC_combined_temp$term_id))
df_NUC_combined_temp <- df_NUC_combined_temp[!rows_to_remove, ]
df_NUC_combined_temp_2 <- df_NUC_combined_temp[,-c(19:ncol(df_NUC_combined_temp))]

write.table(df_NUC_combined_temp_2, "5_input_Prohits_Boris_baits_NUC_gProfileroutput_BPCC_curated_wNoNUCign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)


#-------
# PML body
baits_PML <- bait_input[bait_input$Primary.annotation == "PML body",]
namelist_PML <- baits_PML$Bait.name.in.SAINT.results.6099..supp.table.3.
write.table(baits_PML, "1_Boris_baits_PML_list.txt", sep="\t", row.names=FALSE)

df_PML_all <- data.frame()
df_PML_topn <- data.frame()
wb_PML = createWorkbook()

for (i in namelist_PML) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_PML, sheet_name)
  writeData(wb_PML, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_PML_all <- rbind(df2, df_PML_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_PML_topn <- rbind(df_topn, df_PML_topn)
}

saveWorkbook(wb_PML, 'gprofiler_PML_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_PML_all$parents <- as.character(df_PML_all$parents)
df_PML_topn$parents <- as.character(df_PML_topn$parents)
write.table(df_PML_all, "2_Boris_baits_PML_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_PML_topn, "2_Boris_baits_PML_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_PML_noFilter_all <- data.frame()

for (i in namelist_PML) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_PML_noFilter_all <- rbind(df_notsign, df_PML_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_PML_noFilter_topn <- df_PML_noFilter_all[df_PML_noFilter_all$term_id %in% df_PML_topn$term_id,]
df_PML_noFilter_topn$parents <- as.character(df_PML_noFilter_topn$parents)
df_PML_noFilter_topn$signif_score <- ifelse(df_PML_noFilter_topn$significant == TRUE, 1, 0)
df_PML_noFilter_topn$negative_log10_padj <- -log10(df_PML_noFilter_topn$adjusted_p_value)

write.table(df_PML_noFilter_topn, "3_Boris_baits_PML_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_PML_custom_all <- data.frame()
df_PML_custom_curated <- data.frame()
wb_custom_PML = createWorkbook()

for (i in namelist_PML) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_PML, sheet_name)
  writeData(wb_custom_PML, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_PML_custom_all <- rbind(df2_custom, df_PML_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_PML_custom_curated <- rbind(df_curated, df_PML_custom_curated)
  
}

saveWorkbook(wb_custom_PML, 'gprofiler_PML_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_PML_custom_all$parents <- as.character(df_PML_custom_all$parents)
df_PML_custom_all <- df_PML_custom_all[df_PML_custom_all$significant == TRUE,]
df_PML_custom_curated$parents <- as.character(df_PML_custom_curated$parents)
df_PML_custom_curated$signif_score <- ifelse(df_PML_custom_curated$significant == TRUE, 1, 0)
df_PML_custom_curated$negative_log10_padj <- -log10(df_PML_custom_curated$adjusted_p_value)

write.table(df_PML_custom_all, "2_Boris_baits_PML_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_PML_custom_curated, "3_Boris_baits_PML_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_PML_combined <- rbind(df_PML_noFilter_topn, df_PML_custom_curated)
write.table(df_PML_combined, "4_input_Prohits_Boris_baits_PML_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
PML_parents <- data.frame(str_split_fixed(df_PML_combined$parents, "\\|", 4))
df_PML_combined_temp <- cbind(df_PML_combined,PML_parents)
rows_to_remove <- apply(df_PML_combined_temp[, 19:ncol(df_PML_combined_temp)], 1, function(row) any(row %in% df_PML_combined_temp$term_id))
df_PML_combined_temp <- df_PML_combined_temp[!rows_to_remove, ]
df_PML_combined_temp_2 <- df_PML_combined_temp[,-c(19:ncol(df_PML_combined_temp))]

write.table(df_PML_combined_temp_2, "5_input_Prohits_Boris_baits_PML_gProfileroutput_BPCC_curated_wNoPMLign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)


#-------
# Other
baits_other <- bait_input[!(bait_input$Primary.annotation == "Nuclear speckle" |
                            bait_input$Primary.annotation == "Paraspeckle" |
                            bait_input$Primary.annotation == "Nucleolus" |
                            bait_input$Primary.annotation == "Cajal body" |
                            bait_input$Primary.annotation == "PML body"),]
namelist_other <- baits_other$Bait.name.in.SAINT.results.6099..supp.table.3.
write.table(baits_other, "1_Boris_baits_other_list.txt", sep="\t", row.names=FALSE)

df_other_all <- data.frame()
df_other_topn <- data.frame()
wb_other = createWorkbook()

for (i in namelist_other) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_other, sheet_name)
  writeData(wb_other, sheet_name, df)
  df2 <- df[(grepl("^GO:", df$term_id) | 
               df$term_id == "Dyakov:Nuclear speckle curated" |
               df$term_id == "Dyakov:Cajal body curated" |
               df$term_id == "Dyakov:Paraspeckle curated" |
               df$term_id == "Dyakov:PML body curated" |
               df$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_other_all <- rbind(df2, df_other_all)
  df_topn <- df2[df2$term_size <= 500 & df2$term_size >= 25 & df2$highlighted == TRUE,]
  df_topn <- head(df_topn, 5)
  df_other_topn <- rbind(df_topn, df_other_topn)
}

saveWorkbook(wb_other, 'gprofiler_other_Boris_individual_baits_hsapiens_BPCC.xlsx', overwrite = TRUE)
df_other_all$parents <- as.character(df_other_all$parents)
df_other_topn$parents <- as.character(df_other_topn$parents)
write.table(df_other_all, "2_Boris_baits_other_gProfileroutput_BPCC.txt", sep="\t", row.names=FALSE)
write.table(df_other_topn, "2_Boris_baits_other_gProfileroutput_BPCC_highlighted_top5.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Outputting all results (for viz)
df_other_noFilter_all <- data.frame()

for (i in namelist_other) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
                    organism = "hsapiens", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df_notsign <- data.frame(gprofiler$result)
  colnames(df_notsign)[3] <- "adjusted_p_value"
  df_notsign$query <- i
  df_notsign <- df_notsign[order(df_notsign$adjusted_p_value), ]
  df_notsign$parents <- lapply(df_notsign$parents, function(x){paste0(x,collapse='|')})
  df_other_noFilter_all <- rbind(df_notsign, df_other_noFilter_all)
}

# Filter by topn list to get all overlapping 
df_other_noFilter_topn <- df_other_noFilter_all[df_other_noFilter_all$term_id %in% df_other_topn$term_id,]
df_other_noFilter_topn$parents <- as.character(df_other_noFilter_topn$parents)
df_other_noFilter_topn$signif_score <- ifelse(df_other_noFilter_topn$significant == TRUE, 1, 0)
df_other_noFilter_topn$negative_log10_padj <- -log10(df_other_noFilter_topn$adjusted_p_value)

write.table(df_other_noFilter_topn, "3_Boris_baits_other_gProfileroutput_BPCC_highlighted_top5_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)


# Rerun with custom term and add to 
df_other_custom_all <- data.frame()
df_other_custom_curated <- data.frame()
wb_custom_other = createWorkbook()

for (i in namelist_other) {
  gene <- SAINT_input[SAINT_input$Bait == i & SAINT_input$BFDR <= 0.01,]
  if (nrow(gene) == 0) {
    next  # Skip this iteration if the query is empty
  }
  query <- gene$PreyGene
  gprofiler <- gost(query = query,
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
  addWorksheet(wb_custom_other, sheet_name)
  writeData(wb_custom_other, sheet_name, df_custom)
  df2_custom <- df_custom[grepl("^GO:|^Dyakov:", df_custom$term_id),]
  df_other_custom_all <- rbind(df2_custom, df_other_custom_all)
  df_curated <- df_custom[(df_custom$term_id == "Dyakov:Nuclear speckle curated" |
                             df_custom$term_id == "Dyakov:Cajal body curated" |
                             df_custom$term_id == "Dyakov:Paraspeckle curated" |
                             df_custom$term_id == "Dyakov:PML body curated" |
                             df_custom$term_id == "Dyakov:Comprehensive nucleolus"), ]
  df_other_custom_curated <- rbind(df_curated, df_other_custom_curated)
  
}

saveWorkbook(wb_custom_other, 'gprofiler_other_Boris_individual_baits_customtoken_BPCC.xlsx', overwrite = TRUE)
df_other_custom_all$parents <- as.character(df_other_custom_all$parents)
df_other_custom_all <- df_other_custom_all[df_other_custom_all$significant == TRUE,]
df_other_custom_curated$parents <- as.character(df_other_custom_curated$parents)
df_other_custom_curated$signif_score <- ifelse(df_other_custom_curated$significant == TRUE, 1, 0)
df_other_custom_curated$negative_log10_padj <- -log10(df_other_custom_curated$adjusted_p_value)

write.table(df_other_custom_all, "2_Boris_baits_other_gProfileroutput_BPCC_customToken.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_other_custom_curated, "3_Boris_baits_other_gProfileroutput_BPCC_customToken_curated_wNonSign.txt", sep="\t", row.names=FALSE, quote=FALSE)

# rBind results to get input
df_other_combined <- rbind(df_other_noFilter_topn, df_other_custom_curated)
write.table(df_other_combined, "4_input_Prohits_Boris_baits_other_gProfileroutput_BPCC_curated_wNonSign_wCurated.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out any rows that are children terms of term already in file
other_parents <- data.frame(str_split_fixed(df_other_combined$parents, "\\|", 4))
df_other_combined_temp <- cbind(df_other_combined,other_parents)
rows_to_remove <- apply(df_other_combined_temp[, 19:ncol(df_other_combined_temp)], 1, function(row) any(row %in% df_other_combined_temp$term_id))
df_other_combined_temp <- df_other_combined_temp[!rows_to_remove, ]
df_other_combined_temp_2 <- df_other_combined_temp[,-c(19:ncol(df_other_combined_temp))]

write.table(df_other_combined_temp_2, "5_input_Prohits_Boris_baits_other_gProfileroutput_BPCC_curated_wNootherign_wCurated_removeChildren.txt", sep="\t", row.names=FALSE, quote=FALSE)

