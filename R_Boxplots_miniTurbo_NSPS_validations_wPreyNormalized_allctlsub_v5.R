# Plot boxplots showing normalized SPCs across NS/PS preys (from NMF results)

## Part1: Create matrix of MSPIT_SPCs for NS/PS curated preys across predicted NS/PS miniTurbo preys (normalize each prey 0 to 1 relative to the max score)
## Part2: Filter for only baits with signficant hits and annotated in NMF
## Part3: Plot boxplots showing avgSPC for NS and PS preys (subtracting avgSPC in controls + normalized from 0 to 1 (divide by max))
## Part4: Plot barplots showing number of NS/PS preys recovered for each bait

#----------
## Part1: Create matrix of MSPIT_SPCs across all miniTurbo preys (filter for NS/PS primary annotation?)
SPC_input <- read.delim("SAINT_miniTurbo_Raktan_results_ID_7015_list_tmp.txt", sep="\t", header=TRUE)

# Read in NMF matrix and output list of NS and PS proteins
NMF_matrix <- read.delim("input_0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", header=TRUE)

NMF_matrix_NS <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X12" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X12" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X12"),] 

NMF_matrix_PS <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X0" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X0" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X0"),] 

NMF_both <- intersect(NMF_matrix_NS$PreyGene, NMF_matrix_PS$PreyGene)

# Create matrix with the control counts and calculate average, add to SAINT file
library(stringr)
control_counts <- data.frame(str_split_fixed(SPC_input$ctrlCounts, "\\|", 27))
control_counts <- as.data.frame(sapply(control_counts, as.numeric))
num_controls <- ncol(control_counts)

# Take the average of all controls 
control_counts$AvgSpec_controls <- apply(control_counts[, c(1:ncol(control_counts))], 1, FUN=mean, na.rm=TRUE)
control_counts$AvgSpec_controls <- as.numeric(sprintf(control_counts$AvgSpec_controls, fmt = '%#.1f'))

# add to SAINT file (if negative set AvgSpec at 0)
SPC_input$AvgSpec_controls <- control_counts$AvgSpec_controls
SPC_input$AvgSpec_minusCtls <- SPC_input$AvgSpec - SPC_input$AvgSpec_controls
SPC_input$AvgSpec_minusCtls <- ifelse(SPC_input$AvgSpec_minusCtls >= 0, SPC_input$AvgSpec_minusCtls, 0)
SPC_input$Prey_name_ID <- paste0(SPC_input$PreyGene, "_", SPC_input$Prey)

# Add annotation for NS and PS based on NMF results
SPC_input_PS <- SPC_input[SPC_input$PreyGene %in% NMF_matrix_PS$PreyGene,]
SPC_input_PS <- SPC_input_PS[!(SPC_input_PS$PreyGene %in% NMF_both),]
SPC_input_PS$PSNS_annotation <- "PS"

SPC_input_NS <- SPC_input[SPC_input$PreyGene %in% NMF_matrix_NS$PreyGene,]
SPC_input_NS <- SPC_input_NS[!(SPC_input_NS$PreyGene %in% NMF_both),]
SPC_input_NS$PSNS_annotation <- "NS"

SPC_input_both <- SPC_input[SPC_input$PreyGene %in% NMF_both,]
SPC_input_both$PSNS_annotation <- "Both"

SPC_input_combined <- rbind(SPC_input_PS,SPC_input_NS,SPC_input_both) 
  
# Combine and expand to create matrix with SPC for each prey across all baits 
SPC_input_SPC <- SPC_input_combined[, c(1,25,24,26)]
SPC_input_BFDR <- SPC_input_combined[, c(1,25,16,26)]

# Expand (make horizontal) and annotate all "NA" in SPC as 0 (no values)
library(tidyr)
SPC_input_SPC_wide <- spread(SPC_input_SPC, Bait, AvgSpec_minusCtls)
SPC_input_SPC_wide[is.na(SPC_input_SPC_wide)] <- 0
SPC_input_BFDR_wide <- spread(SPC_input_BFDR, Bait, BFDR)
SPC_input_BFDR_wide[is.na(SPC_input_BFDR_wide)] <- 1

write.table(SPC_input_SPC_wide, "1_miniTurbo_7015_validations_avgSPC_subALLctls_SPC_NSPS_table.txt", sep="\t", row.names=FALSE)
write.table(SPC_input_BFDR_wide, "1_miniTurbo_SAINT_7015_validations_avgSPC_subALLctls_BFDR_NSPS_table.txt", sep="\t", row.names=FALSE)

#---------
## Part2: Filter data and normalize from 0 to 1 across preys
# Create list of baits that overlap NMF data
baits_SAINT <- data.frame(str_split_fixed(SPC_input$Bait, "_", 2))
baits_SAINT$Bait <- SPC_input$Bait
colnames(baits_SAINT)[1:2] <- c("BaitGene","BaitTag")
baits_SAINT_uniq <- unique(baits_SAINT)

# Data frame of preys with NS annotations
baits_in_NMF <- baits_SAINT_uniq[baits_SAINT_uniq$BaitGene %in% NMF_matrix_NS$PreyGene | baits_SAINT_uniq$BaitGene %in% NMF_matrix_PS$PreyGene,]

SPC_input_SPC_wide_filtered <- SPC_input_SPC_wide[,baits_in_NMF$Bait]
SPC_input_SPC_wide_filtered <- cbind(SPC_input_SPC_wide[,c(1,2)],SPC_input_SPC_wide_filtered)

SPC_input_BFDR_wide_filtered <- SPC_input_BFDR_wide[,baits_in_NMF$Bait]
SPC_input_BFDR_wide_filtered <- cbind(SPC_input_BFDR_wide[,c(1,2)],SPC_input_BFDR_wide_filtered)

# Manually remove any rows without significant hits
columns_to_remove <- c("PPP1CA_CmT", "RPRD2_CmT", "RSBN1L_CmT", "TASOR_CmT")
SPC_input_SPC_wide_filtered <- SPC_input_SPC_wide_filtered[,!(names(SPC_input_SPC_wide_filtered) %in% columns_to_remove)]

# Output file with dotplot info for selected baits (as preys)
NMF_matrix_Prohits_input <- NMF_matrix[NMF_matrix$PreyGene %in% baits_SAINT_uniq$BaitGene,]
write.table(NMF_matrix_Prohits_input, "0.5_Prohits_input_NSPS_miniTurbo_Raktan_Baits_NMF_scores.txt", sep="\t", row.names=FALSE)


#-------
## Part3: Normalize the avgSPCs relative to the row MAX value and plot values
# Normalize each avgSPC row between 0 to 1 relative to the MAX value 
normSPC_input <- SPC_input_SPC_wide_filtered
normSPC_input_preys <- normSPC_input$Prey_name_ID
normSPC_input_NSPS <- normSPC_input$PSNS_annotation
normSPC_input$Max_avgSPC <- apply(normSPC_input[,c(3:ncol(normSPC_input))], 1, FUN=max, na.rm = TRUE)
normSPC_input_2 <- normSPC_input[,c(3:(ncol(normSPC_input)-1))] / normSPC_input$Max_avgSPC 

normSPC_input_3 <- cbind(normSPC_input_preys,normSPC_input_NSPS,normSPC_input_2)
colnames(normSPC_input_3)[1:2] <- c("PreyGene_ID","PSNS_annotation")

# Remove any rows that are all 0 (NaN)
normSPC_input_3 <- normSPC_input_3[!normSPC_input_3$DDX17_CmT == "NaN",]
write.table(normSPC_input_3, "2_miniTurbo_7015_validations_normalized_avgSPC_subALLctls_NSPS_table.txt", sep="\t", row.names=FALSE)

#---------
# Re-plot the boxplots with the row normalized avgSPCs (relative to highest Max value across all baits)
# Create input matrices
library("ggplot2")
library("reshape2")

# Manually specify order
bait_order <- c("PreyGene_ID",
               "PSNS_annotation",
               "DDX17_CmT",
               "RNMT_CmT",
               "STRBP_CmT",
               "RBMX_CmT",
               "NONO_NmT",
               "RPRD1B_CmT",
               "RBFOX2_CmT",
               "PTBP3_CmT",
               "MBNL1_CmT",
               "LARP7_CmT",
               "SRSF7_CmT",
               "PQBP1_CmT",
               "MFAP1_CmT",
               "DHX38_CmT",
               "RBM15_CmT",
               "SAP18_CmT",
               "PPIL4_NmT",
               "RBM15B_CmT",
               "PAXBP1_CmT",
               "CHD2_CmT"
)

input_normSPC <- normSPC_input_3[,bait_order]
input_normSPC <- as.matrix(input_normSPC)
row.names(input_normSPC) <- input_normSPC[,1]
input_normSPC <- input_normSPC[,-c(1,2)]
class(input_normSPC)<-"numeric"

input_normSPC_2 <- melt(input_normSPC)
colnames(input_normSPC_2) <- c("PreyName_ID","Bait","normalized_avgSPC")
input_normSPC_2$NSPS_annotation <- NA
input_normSPC_2$NSPS_annotation[normSPC_input_3$PSNS_annotation == "PS"] <- "PS" 
input_normSPC_2$NSPS_annotation[normSPC_input_3$PSNS_annotation == "NS"] <- "NS" 

input_normSPC_2 <- input_normSPC_2[!is.na(input_normSPC_2$NSPS_annotation),]

# Plot "merged"
pdf("Rplot_7015_mT_validations_NSPS_boxplots_NMFpreys_normSPC_allctlsub.pdf", width = 13, height = 5)
ggplot(input_normSPC_2, aes(x=Bait, y=normalized_avgSPC))+
  geom_boxplot(aes(fill=NSPS_annotation), width=0.8, outlier.alpha=FALSE, position=position_dodge2(padding = 0.3, preserve = c("total", "single")))+
  scale_fill_manual(values=c(NS="#7C96C2", PS="#F180BA"))+
  ylim(0,1.2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title="Relative avgSPC for NS and PS preys across validation baits (prey normalized 0 to 1 across baits (lowest=0 to highest=1), ctl subtraction[ALL]")+
  xlab("Baits")+
  ylab("Relative avgSPC across baits")+
  #theme_classic()
  theme_bw()
dev.off()

pdf("Rplot_7015_mT_validations_NSPS_boxplots_NMFpreys_normSPC_allctlsub_wOutliers.pdf", width = 13, height = 5)
ggplot(input_normSPC_2, aes(x=Bait, y=normalized_avgSPC))+
  geom_boxplot(aes(fill=NSPS_annotation), width=0.8, outlier.alpha=TRUE, position=position_dodge2(padding = 0.3, preserve = c("total", "single")))+
  scale_fill_manual(values=c(NS="#7C96C2", PS="#F180BA"))+
  ylim(0,1.2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title="Relative avgSPC for NS and PS preys across validation baits (prey normalized 0 to 1 across baits (lowest=0 to highest=1), ctl subtraction[ALL])")+
  xlab("Baits")+
  ylab("Relative avgSPC across baits")+
  #theme_classic()
  theme_bw()
dev.off()

# Separated
pdf("Rplot_7015_mT_validations_NSPS_boxplots_NMFpreys_normSPC_allctlsub_individual.pdf", width = 20, height = 20)
ggplot(input_normSPC_2, aes(x=Bait, y=normalized_avgSPC))+
  geom_boxplot(aes(fill=NSPS_annotation), width=0.5, outlier.alpha=TRUE)+
  scale_fill_manual(values=c(NS="#7C96C2", PS="#F180BA"))+
  ylim(0,1.2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title="Relative avgSPC for NS and PS preys across validation baits (prey normalized 0 to 1 across baits (lowest=0 to highest=1), ctl subtraction[ALL])")+
  xlab("Baits")+
  ylab("Relative avgSPC across baits")+
  #theme_classic()
  theme_bw()+
  facet_wrap(~Bait,scales = "free")
dev.off()

# Statistical tests
normSPC_input_PS <- normSPC_input_3[normSPC_input_3$PSNS_annotation == "PS",]
names(normSPC_input_PS) <- paste(names(normSPC_input_PS), "_PS", sep = "")

normSPC_input_NS <- normSPC_input_3[normSPC_input_3$PSNS_annotation == "NS",]
names(normSPC_input_NS) <- paste(names(normSPC_input_NS), "_NS", sep = "")

normSPC_input_PSNS_forsign <- merge(normSPC_input_PS,normSPC_input_NS,by.x = "PreyGene_ID_PS" , by.y = "PreyGene_ID_NS", all.x=TRUE, all.y=TRUE)

# Run t.test or wilcox.test for each combination
column_pairs <- list(c("RBMX_CmT_NS", "RBMX_CmT_PS"), 
                     c("NONO_NmT_NS", "NONO_NmT_PS"),
                     c("DDX17_CmT_NS", "DDX17_CmT_PS"),
                     c("MBNL1_CmT_NS", "MBNL1_CmT_PS"),
                     c("PTBP3_CmT_NS", "PTBP3_CmT_PS"),
                     c("RNMT_CmT_NS", "RNMT_CmT_PS"),
                     c("RPRD1B_CmT_NS", "RPRD1B_CmT_PS"),
                     c("RBFOX2_CmT_NS", "RBFOX2_CmT_PS"),
                     c("LARP7_CmT_NS", "LARP7_CmT_PS"),
                     c("STRBP_CmT_NS", "STRBP_CmT_PS"),
                     c("SRSF7_CmT_NS", "SRSF7_CmT_PS"),
                     c("SAP18_CmT_NS", "SAP18_CmT_PS"),
                     c("PQBP1_CmT_NS", "PQBP1_CmT_PS"),
                     c("RBM15B_CmT_NS", "RBM15B_CmT_PS"),
                     c("RBM15_CmT_NS", "RBM15_CmT_PS"),
                     c("DHX38_CmT_NS", "DHX38_CmT_PS"),
                     c("PAXBP1_CmT_NS", "PAXBP1_CmT_PS"),
                     c("PPIL4_NmT_NS", "PPIL4_NmT_PS"),
                     c("MFAP1_CmT_NS", "MFAP1_CmT_PS"),
                     c("CHD2_CmT_NS", "CHD2_CmT_PS")
                     )

# Function to perform t-test for each column pair
perform_t_test <- function(df, col1, col2) {
  t_test_result <- t.test(df[[col1]], df[[col2]])
  wilcox_result <- wilcox.test(df[[col1]], df[[col2]])
  return(list(t_test_result = t_test_result, wilcox_result = wilcox_result))
}

# Perform t-test for each column pair and store results in a dataframe
results_df <- data.frame(
  Column1 = character(length(column_pairs)),
  Column2 = character(length(column_pairs)),
  mean_Column1 = character(length(column_pairs)),
  mean_Column2 = character(length(column_pairs)),
  t.test_p_value = numeric(length(column_pairs)),
  wilcox_p_value = numeric(length(column_pairs))
)

for (i in seq_along(column_pairs)) {
  test_results <- perform_t_test(normSPC_input_PSNS_forsign, column_pairs[[i]][1], column_pairs[[i]][2])
  results_df[i, "Column1"] <- column_pairs[[i]][1]
  results_df[i, "Column2"] <- column_pairs[[i]][2]
  results_df[i, "mean_Column1"] <- test_results$t_test_result$estimate[["mean of x"]]
  results_df[i, "mean_Column2"] <- test_results$t_test_result$estimate[["mean of y"]]
  results_df[i, "t.test_p_value"] <- test_results$t_test_result$p.value
  results_df[i, "wilcox_p_value"] <- test_results$wilcox_result$p.value
} 

write.table(results_df, "4_miniTurbo_7015_validations_normalized_avgSPC_subALLctls_NSPS_significance.txt", sep="\t", row.names=FALSE)

#----------
## Part4: Plot barplots showing % of NS/PS preys recovered for each bait
# Filter column by value
SPC_input_BFDR_PS_wide <- SPC_input_BFDR_wide[SPC_input_BFDR_wide$PSNS_annotation == "PS",]
SPC_input_BFDR_NS_wide <- SPC_input_BFDR_wide[SPC_input_BFDR_wide$PSNS_annotation == "NS",]

# PS
BFDR_PS_0.01_RBMX_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RBMX_CmT <= 0.01),]
BFDR_PS_0.01_NONO_NmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$NONO_NmT <= 0.01),]
BFDR_PS_0.01_DDX17_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$DDX17_CmT <= 0.01),]
BFDR_PS_0.01_MBNL1_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$MBNL1_CmT <= 0.01),]
BFDR_PS_0.01_PTBP3_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$PTBP3_CmT <= 0.01),]
BFDR_PS_0.01_RNMT_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RNMT_CmT <= 0.01),]
BFDR_PS_0.01_RPRD1B_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RPRD1B_CmT <= 0.01),]
BFDR_PS_0.01_RBFOX2_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RBFOX2_CmT <= 0.01),]
BFDR_PS_0.01_LARP7_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$LARP7_CmT <= 0.01),]
BFDR_PS_0.01_STRBP_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$STRBP_CmT <= 0.01),]

BFDR_PS_0.01_SRSF7_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$SRSF7_CmT <= 0.01),]
BFDR_PS_0.01_SAP18_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$SAP18_CmT <= 0.01),]
BFDR_PS_0.01_PQBP1_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$PQBP1_CmT <= 0.01),]
BFDR_PS_0.01_RBM15B_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RBM15B_CmT <= 0.01),]
BFDR_PS_0.01_RBM15_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$RBM15_CmT <= 0.01),]
BFDR_PS_0.01_DHX38_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$DHX38_CmT <= 0.01),]
BFDR_PS_0.01_PAXBP1_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$PAXBP1_CmT <= 0.01),]
BFDR_PS_0.01_PPIL4_NmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$PPIL4_NmT <= 0.01),]
BFDR_PS_0.01_MFAP1_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$MFAP1_CmT <= 0.01),]
BFDR_PS_0.01_CHD2_CmT <- SPC_input_BFDR_PS_wide[which(SPC_input_BFDR_PS_wide$CHD2_CmT <= 0.01),]

# NS
BFDR_NS_0.01_RBMX_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RBMX_CmT <= 0.01),]
BFDR_NS_0.01_NONO_NmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$NONO_NmT <= 0.01),]
BFDR_NS_0.01_DDX17_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$DDX17_CmT <= 0.01),]
BFDR_NS_0.01_MBNL1_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$MBNL1_CmT <= 0.01),]
BFDR_NS_0.01_PTBP3_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$PTBP3_CmT <= 0.01),]
BFDR_NS_0.01_RNMT_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RNMT_CmT <= 0.01),]
BFDR_NS_0.01_RPRD1B_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RPRD1B_CmT <= 0.01),]
BFDR_NS_0.01_RBFOX2_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RBFOX2_CmT <= 0.01),]
BFDR_NS_0.01_LARP7_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$LARP7_CmT <= 0.01),]
BFDR_NS_0.01_STRBP_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$STRBP_CmT <= 0.01),]

BFDR_NS_0.01_SRSF7_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$SRSF7_CmT <= 0.01),]
BFDR_NS_0.01_SAP18_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$SAP18_CmT <= 0.01),]
BFDR_NS_0.01_PQBP1_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$PQBP1_CmT <= 0.01),]
BFDR_NS_0.01_RBM15B_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RBM15B_CmT <= 0.01),]
BFDR_NS_0.01_RBM15_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$RBM15_CmT <= 0.01),]
BFDR_NS_0.01_DHX38_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$DHX38_CmT <= 0.01),]
BFDR_NS_0.01_PAXBP1_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$PAXBP1_CmT <= 0.01),]
BFDR_NS_0.01_PPIL4_NmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$PPIL4_NmT <= 0.01),]
BFDR_NS_0.01_MFAP1_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$MFAP1_CmT <= 0.01),]
BFDR_NS_0.01_CHD2_CmT <- SPC_input_BFDR_NS_wide[which(SPC_input_BFDR_NS_wide$CHD2_CmT <= 0.01),]


# Create matrices filtered by BFDR to get number of preys identified
mT_baits <- c("RBMX_CmT",
                "NONO_NmT",
                "DDX17_CmT",
                "MBNL1_CmT",
                "PTBP3_CmT",
                "RNMT_CmT",
                "RPRD1B_CmT",
                "RBFOX2_CmT",
                "LARP7_CmT",
                "STRBP_CmT",
                "SRSF7_CmT",
                "SAP18_CmT",
                "PQBP1_CmT",
                "RBM15B_CmT",
                "RBM15_CmT",
                "DHX38_CmT",
                "PAXBP1_CmT",
                "PPIL4_NmT",
                "MFAP1_CmT",
                "CHD2_CmT"
)


NMF_localization <- c("PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "PS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS",
                           "NS"
)

# Calculate percent recovery of known preys
PS_preys_recovered <- c(length(BFDR_PS_0.01_RBMX_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_NONO_NmT$Prey_name_ID),
                        length(BFDR_PS_0.01_DDX17_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_MBNL1_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_PTBP3_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_RNMT_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_RPRD1B_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_RBFOX2_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_LARP7_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_STRBP_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_SRSF7_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_SAP18_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_PQBP1_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_RBM15B_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_RBM15_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_DHX38_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_PAXBP1_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_PPIL4_NmT$Prey_name_ID),
                        length(BFDR_PS_0.01_MFAP1_CmT$Prey_name_ID),
                        length(BFDR_PS_0.01_CHD2_CmT$Prey_name_ID)
)

PS_total_preys_NMF <- c(length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID),
                    length(SPC_input_BFDR_PS_wide$Prey_name_ID)
)

NS_preys_recovered <- c(length(BFDR_NS_0.01_RBMX_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_NONO_NmT$Prey_name_ID),
                        length(BFDR_NS_0.01_DDX17_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_MBNL1_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_PTBP3_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_RNMT_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_RPRD1B_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_RBFOX2_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_LARP7_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_STRBP_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_SRSF7_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_SAP18_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_PQBP1_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_RBM15B_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_RBM15_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_DHX38_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_PAXBP1_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_PPIL4_NmT$Prey_name_ID),
                        length(BFDR_NS_0.01_MFAP1_CmT$Prey_name_ID),
                        length(BFDR_NS_0.01_CHD2_CmT$Prey_name_ID)
)

NS_total_preys_NMF <- c(length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID),
                        length(SPC_input_BFDR_NS_wide$Prey_name_ID)
)



# Create input for barplots
df_number_recovery <- data.frame(mT_baits, NMF_localization, PS_total_preys_NMF, PS_preys_recovered, NS_total_preys_NMF, NS_preys_recovered)
df_number_recovery[, 3:ncol(df_number_recovery)] <- lapply(3:ncol(df_number_recovery), function(x) as.numeric(df_number_recovery[[x]]))
df_number_recovery$PS_percent_recovery <- df_number_recovery$PS_preys_recovered / df_number_recovery$PS_total_preys_NMF*100
df_number_recovery$NS_percent_recovery <- df_number_recovery$NS_preys_recovered / df_number_recovery$NS_total_preys_NMF*100

write.table(df_number_recovery, "3_NSPS_number_preys_recovered.txt", sep="\t", row.names=FALSE)

input_prey_recovery <- df_number_recovery[,c(1,2,4,6)]
input_prey_recovery$mT_baits <- as.factor(input_prey_recovery$mT_baits)
input_prey_recovery_2 <- melt(input_prey_recovery)
input_prey_recovery_2$variable <- factor(input_prey_recovery_2$variable, levels = c("NS_preys_recovered", "PS_preys_recovered"))

# Plot number of NS and PS preys recovered
pdf("Rplot_NSPSmerged_barplots_NumberPreysRecovered_allctlsub.pdf", width = 13, height = 5)
ggplot(input_prey_recovery_2, aes(x=factor(mT_baits,bait_order), y=value, fill=variable))+
  geom_bar(stat="identity", width = 0.8, position=position_dodge(), colour="black")+ 
  scale_fill_manual(values=c(PS_preys_recovered="#F180BA", NS_preys_recovered="#7C96C2"))+
  ylim(0,100)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title="Number of NS and PS preys recovered (NMF annotations) across validation baits (BFDR <= 0.01")+
  xlab("Baits")+
  ylab("Number of preys recovered in rank (PS=X0, NS=X12)")+
  #theme_classic()
  theme_bw()
dev.off()
