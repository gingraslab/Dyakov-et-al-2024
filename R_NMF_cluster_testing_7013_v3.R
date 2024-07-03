# Outputting Precision/Recall and Enrichment plots to assess NMF cluster numbers

##Part 1: Create columns with primary, secondary and tertiary ranks (use primary annotation only)
##Part 2: Run gProfiler for each of the samples using Boris' custom terms
##Part 3: Calculate precision/recall and enrichments and plot across clusters

#---------
##Part 1: Create columns with primary, secondary and tertiary ranks (primary annotation only)
#Read in NMF files and re-annotate
Boris_curated <- read.delim("Boris_curated_lists_wIF_2.csv", sep=",",header = TRUE)

NMF_matrix_10 <- read.delim("7013_cleaned_v2_NMF_10_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_10 <- NMF_matrix_10[,-(ncol(NMF_matrix_10))]
NMF_matrix_11 <- read.delim("7013_cleaned_v2_NMF_11_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_11 <- NMF_matrix_11[,-(ncol(NMF_matrix_11))]
NMF_matrix_12 <- read.delim("7013_cleaned_v2_NMF_12_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_12 <- NMF_matrix_12[,-(ncol(NMF_matrix_12))]
NMF_matrix_13 <- read.delim("7013_cleaned_v2_NMF_13_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_13 <- NMF_matrix_13[,-(ncol(NMF_matrix_13))]
NMF_matrix_14 <- read.delim("7013_cleaned_v2_NMF_14_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_14 <- NMF_matrix_14[,-(ncol(NMF_matrix_14))]
NMF_matrix_15 <- read.delim("7013_cleaned_v2_NMF_15_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_15 <- NMF_matrix_15[,-(ncol(NMF_matrix_15))]
NMF_matrix_16 <- read.delim("7013_cleaned_v2_NMF_16_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_16 <- NMF_matrix_16[,-(ncol(NMF_matrix_16))]
NMF_matrix_17 <- read.delim("7013_cleaned_v2_NMF_17_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_17 <- NMF_matrix_17[,-(ncol(NMF_matrix_17))]
NMF_matrix_18 <- read.delim("7013_cleaned_v2_NMF_18_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_18 <- NMF_matrix_18[,-(ncol(NMF_matrix_18))]
NMF_matrix_19 <- read.delim("7013_cleaned_v2_NMF_19_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_19 <- NMF_matrix_19[,-(ncol(NMF_matrix_19))]
NMF_matrix_20 <- read.delim("7013_cleaned_v2_NMF_20_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_20 <- NMF_matrix_20[,-(ncol(NMF_matrix_20))]
NMF_matrix_21 <- read.delim("7013_cleaned_v2_NMF_21_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_21 <- NMF_matrix_21[,-(ncol(NMF_matrix_21))]
NMF_matrix_22 <- read.delim("7013_cleaned_v2_NMF_22_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_22 <- NMF_matrix_22[,-(ncol(NMF_matrix_22))]
NMF_matrix_23 <- read.delim("7013_cleaned_v2_NMF_23_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_23 <- NMF_matrix_23[,-(ncol(NMF_matrix_23))]
NMF_matrix_24 <- read.delim("7013_cleaned_v2_NMF_24_preys_matrix_results-top10ctrls.csv", sep=",",header = TRUE)
NMF_matrix_24 <- NMF_matrix_24[,-(ncol(NMF_matrix_24))]


# Annotate the top, second and third ranked scores for each file
# Create functions to output 2nd and third largest values
func_second <- function(x){
  as.numeric(sort(x, TRUE)[2])
}  

func_third <- function(x){
  as.numeric(sort(x, TRUE)[3])
}  

# Add column with primary score
NMF_matrix_10$Primary_score <- apply(NMF_matrix_10[, c(2:ncol(NMF_matrix_10))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_11$Primary_score <- apply(NMF_matrix_11[, c(2:ncol(NMF_matrix_11))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_12$Primary_score <- apply(NMF_matrix_12[, c(2:ncol(NMF_matrix_12))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_13$Primary_score <- apply(NMF_matrix_13[, c(2:ncol(NMF_matrix_13))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_14$Primary_score <- apply(NMF_matrix_14[, c(2:ncol(NMF_matrix_14))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_15$Primary_score <- apply(NMF_matrix_15[, c(2:ncol(NMF_matrix_15))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_16$Primary_score <- apply(NMF_matrix_16[, c(2:ncol(NMF_matrix_16))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_17$Primary_score <- apply(NMF_matrix_17[, c(2:ncol(NMF_matrix_17))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_18$Primary_score <- apply(NMF_matrix_18[, c(2:ncol(NMF_matrix_18))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_19$Primary_score <- apply(NMF_matrix_19[, c(2:ncol(NMF_matrix_19))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_20$Primary_score <- apply(NMF_matrix_20[, c(2:ncol(NMF_matrix_20))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_21$Primary_score <- apply(NMF_matrix_21[, c(2:ncol(NMF_matrix_21))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_22$Primary_score <- apply(NMF_matrix_22[, c(2:ncol(NMF_matrix_22))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_23$Primary_score <- apply(NMF_matrix_23[, c(2:ncol(NMF_matrix_23))], 1, FUN=max, na.rm=TRUE)
NMF_matrix_24$Primary_score <- apply(NMF_matrix_24[, c(2:ncol(NMF_matrix_24))], 1, FUN=max, na.rm=TRUE)

# Add column with secondary score
NMF_matrix_10$Secondary_score <- apply(NMF_matrix_10[, c(2:(ncol(NMF_matrix_10)-1))], 1, FUN=func_second)
NMF_matrix_11$Secondary_score <- apply(NMF_matrix_11[, c(2:(ncol(NMF_matrix_11)-1))], 1, FUN=func_second)
NMF_matrix_12$Secondary_score <- apply(NMF_matrix_12[, c(2:(ncol(NMF_matrix_12)-1))], 1, FUN=func_second)
NMF_matrix_13$Secondary_score <- apply(NMF_matrix_13[, c(2:(ncol(NMF_matrix_13)-1))], 1, FUN=func_second)
NMF_matrix_14$Secondary_score <- apply(NMF_matrix_14[, c(2:(ncol(NMF_matrix_14)-1))], 1, FUN=func_second)
NMF_matrix_15$Secondary_score <- apply(NMF_matrix_15[, c(2:(ncol(NMF_matrix_15)-1))], 1, FUN=func_second)
NMF_matrix_16$Secondary_score <- apply(NMF_matrix_16[, c(2:(ncol(NMF_matrix_16)-1))], 1, FUN=func_second)
NMF_matrix_17$Secondary_score <- apply(NMF_matrix_17[, c(2:(ncol(NMF_matrix_17)-1))], 1, FUN=func_second)
NMF_matrix_18$Secondary_score <- apply(NMF_matrix_18[, c(2:(ncol(NMF_matrix_18)-1))], 1, FUN=func_second)
NMF_matrix_19$Secondary_score <- apply(NMF_matrix_19[, c(2:(ncol(NMF_matrix_19)-1))], 1, FUN=func_second)
NMF_matrix_20$Secondary_score <- apply(NMF_matrix_20[, c(2:(ncol(NMF_matrix_20)-1))], 1, FUN=func_second)
NMF_matrix_21$Secondary_score <- apply(NMF_matrix_21[, c(2:(ncol(NMF_matrix_21)-1))], 1, FUN=func_second)
NMF_matrix_22$Secondary_score <- apply(NMF_matrix_22[, c(2:(ncol(NMF_matrix_22)-1))], 1, FUN=func_second)
NMF_matrix_23$Secondary_score <- apply(NMF_matrix_23[, c(2:(ncol(NMF_matrix_23)-1))], 1, FUN=func_second)
NMF_matrix_24$Secondary_score <- apply(NMF_matrix_24[, c(2:(ncol(NMF_matrix_24)-1))], 1, FUN=func_second)

# Add column with tertiary score
NMF_matrix_10$Tertiary_score <- apply(NMF_matrix_10[, c(2:(ncol(NMF_matrix_10)-2))], 1, FUN=func_third)
NMF_matrix_11$Tertiary_score <- apply(NMF_matrix_11[, c(2:(ncol(NMF_matrix_11)-2))], 1, FUN=func_third)
NMF_matrix_12$Tertiary_score <- apply(NMF_matrix_12[, c(2:(ncol(NMF_matrix_12)-2))], 1, FUN=func_third)
NMF_matrix_13$Tertiary_score <- apply(NMF_matrix_13[, c(2:(ncol(NMF_matrix_13)-2))], 1, FUN=func_third)
NMF_matrix_14$Tertiary_score <- apply(NMF_matrix_14[, c(2:(ncol(NMF_matrix_14)-2))], 1, FUN=func_third)
NMF_matrix_15$Tertiary_score <- apply(NMF_matrix_15[, c(2:(ncol(NMF_matrix_15)-2))], 1, FUN=func_third)
NMF_matrix_16$Tertiary_score <- apply(NMF_matrix_16[, c(2:(ncol(NMF_matrix_16)-2))], 1, FUN=func_third)
NMF_matrix_17$Tertiary_score <- apply(NMF_matrix_17[, c(2:(ncol(NMF_matrix_17)-2))], 1, FUN=func_third)
NMF_matrix_18$Tertiary_score <- apply(NMF_matrix_18[, c(2:(ncol(NMF_matrix_18)-2))], 1, FUN=func_third)
NMF_matrix_19$Tertiary_score <- apply(NMF_matrix_19[, c(2:(ncol(NMF_matrix_19)-2))], 1, FUN=func_third)
NMF_matrix_20$Tertiary_score <- apply(NMF_matrix_20[, c(2:(ncol(NMF_matrix_20)-2))], 1, FUN=func_third)
NMF_matrix_21$Tertiary_score <- apply(NMF_matrix_21[, c(2:(ncol(NMF_matrix_21)-2))], 1, FUN=func_third)
NMF_matrix_22$Tertiary_score <- apply(NMF_matrix_22[, c(2:(ncol(NMF_matrix_22)-2))], 1, FUN=func_third)
NMF_matrix_23$Tertiary_score <- apply(NMF_matrix_23[, c(2:(ncol(NMF_matrix_23)-2))], 1, FUN=func_third)
NMF_matrix_24$Tertiary_score <- apply(NMF_matrix_24[, c(2:(ncol(NMF_matrix_24)-2))], 1, FUN=func_third)

NMF_matrix_10[,-1] <- round(NMF_matrix_10[,-1], digits=2)
NMF_matrix_11[,-1] <- round(NMF_matrix_11[,-1], digits=2)
NMF_matrix_12[,-1] <- round(NMF_matrix_12[,-1], digits=2)
NMF_matrix_13[,-1] <- round(NMF_matrix_13[,-1], digits=2)
NMF_matrix_14[,-1] <- round(NMF_matrix_14[,-1], digits=2)
NMF_matrix_15[,-1] <- round(NMF_matrix_15[,-1], digits=2)
NMF_matrix_16[,-1] <- round(NMF_matrix_16[,-1], digits=2)
NMF_matrix_17[,-1] <- round(NMF_matrix_17[,-1], digits=2)
NMF_matrix_18[,-1] <- round(NMF_matrix_18[,-1], digits=2)
NMF_matrix_19[,-1] <- round(NMF_matrix_19[,-1], digits=2)
NMF_matrix_20[,-1] <- round(NMF_matrix_20[,-1], digits=2)
NMF_matrix_21[,-1] <- round(NMF_matrix_21[,-1], digits=2)
NMF_matrix_22[,-1] <- round(NMF_matrix_22[,-1], digits=2)
NMF_matrix_23[,-1] <- round(NMF_matrix_23[,-1], digits=2)
NMF_matrix_24[,-1] <- round(NMF_matrix_24[,-1], digits=2)

# Annotate column with primary, secondary and tertiary annotations
# Create function that outputs the matching column name 
maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]

library(dplyr)
# NMF 10
NMF_matrix_10 <- NMF_matrix_10 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_10)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_10 <- NMF_matrix_10 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_10)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_10 <- NMF_matrix_10 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_10)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_10$Secondary_cluster_annotation[NMF_matrix_10$Secondary_score == 0] <- NA
NMF_matrix_10$Tertiary_cluster_annotation[NMF_matrix_10$Tertiary_score == 0] <- NA

# NMF 11
NMF_matrix_11 <- NMF_matrix_11 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_11)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_11 <- NMF_matrix_11 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_11)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_11 <- NMF_matrix_11 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_11)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_11$Secondary_cluster_annotation[NMF_matrix_11$Secondary_score == 0] <- NA
NMF_matrix_11$Tertiary_cluster_annotation[NMF_matrix_11$Tertiary_score == 0] <- NA

# NMF 12
NMF_matrix_12 <- NMF_matrix_12 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_12)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_12 <- NMF_matrix_12 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_12)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_12 <- NMF_matrix_12 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_12)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_12$Secondary_cluster_annotation[NMF_matrix_12$Secondary_score == 0] <- NA
NMF_matrix_12$Tertiary_cluster_annotation[NMF_matrix_12$Tertiary_score == 0] <- NA

# NMF 13
NMF_matrix_13 <- NMF_matrix_13 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_13)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_13 <- NMF_matrix_13 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_13)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_13 <- NMF_matrix_13 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_13)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_13$Secondary_cluster_annotation[NMF_matrix_13$Secondary_score == 0] <- NA
NMF_matrix_13$Tertiary_cluster_annotation[NMF_matrix_13$Tertiary_score == 0] <- NA

# NMF 14
NMF_matrix_14 <- NMF_matrix_14 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_14)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_14 <- NMF_matrix_14 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_14)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_14 <- NMF_matrix_14 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_14)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_14$Secondary_cluster_annotation[NMF_matrix_14$Secondary_score == 0] <- NA
NMF_matrix_14$Tertiary_cluster_annotation[NMF_matrix_14$Tertiary_score == 0] <- NA

# NMF 15
NMF_matrix_15 <- NMF_matrix_15 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_15)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_15 <- NMF_matrix_15 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_15)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_15 <- NMF_matrix_15 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_15)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_15$Secondary_cluster_annotation[NMF_matrix_15$Secondary_score == 0] <- NA
NMF_matrix_15$Tertiary_cluster_annotation[NMF_matrix_15$Tertiary_score == 0] <- NA

# NMF 16
NMF_matrix_16 <- NMF_matrix_16 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_16)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_16 <- NMF_matrix_16 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_16)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_16 <- NMF_matrix_16 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_16)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_16$Secondary_cluster_annotation[NMF_matrix_16$Secondary_score == 0] <- NA
NMF_matrix_16$Tertiary_cluster_annotation[NMF_matrix_16$Tertiary_score == 0] <- NA

# NMF 17
NMF_matrix_17 <- NMF_matrix_17 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_17)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_17 <- NMF_matrix_17 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_17)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_17 <- NMF_matrix_17 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_17)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_17$Secondary_cluster_annotation[NMF_matrix_17$Secondary_score == 0] <- NA
NMF_matrix_17$Tertiary_cluster_annotation[NMF_matrix_17$Tertiary_score == 0] <- NA

# NMF 18
NMF_matrix_18 <- NMF_matrix_18 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_18)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_18 <- NMF_matrix_18 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_18)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_18 <- NMF_matrix_18 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_18)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_18$Secondary_cluster_annotation[NMF_matrix_18$Secondary_score == 0] <- NA
NMF_matrix_18$Tertiary_cluster_annotation[NMF_matrix_18$Tertiary_score == 0] <- NA

# NMF 19
NMF_matrix_19 <- NMF_matrix_19 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_19)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_19 <- NMF_matrix_19 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_19)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_19 <- NMF_matrix_19 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_19)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_19$Secondary_cluster_annotation[NMF_matrix_19$Secondary_score == 0] <- NA
NMF_matrix_19$Tertiary_cluster_annotation[NMF_matrix_19$Tertiary_score == 0] <- NA

# NMF 20
NMF_matrix_20 <- NMF_matrix_20 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_20)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_20 <- NMF_matrix_20 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_20)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_20 <- NMF_matrix_20 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_20)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_20$Secondary_cluster_annotation[NMF_matrix_20$Secondary_score == 0] <- NA
NMF_matrix_20$Tertiary_cluster_annotation[NMF_matrix_20$Tertiary_score == 0] <- NA

# NMF 21
NMF_matrix_21 <- NMF_matrix_21 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_21)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_21 <- NMF_matrix_21 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_21)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_21 <- NMF_matrix_21 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_21)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_21$Secondary_cluster_annotation[NMF_matrix_21$Secondary_score == 0] <- NA
NMF_matrix_21$Tertiary_cluster_annotation[NMF_matrix_21$Tertiary_score == 0] <- NA

# NMF 22
NMF_matrix_22 <- NMF_matrix_22 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_22)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_22 <- NMF_matrix_22 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_22)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_22 <- NMF_matrix_22 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_22)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_22$Secondary_cluster_annotation[NMF_matrix_22$Secondary_score == 0] <- NA
NMF_matrix_22$Tertiary_cluster_annotation[NMF_matrix_22$Tertiary_score == 0] <- NA

# NMF 23
NMF_matrix_23 <- NMF_matrix_23 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_23)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_23 <- NMF_matrix_23 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_23)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_23 <- NMF_matrix_23 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_23)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_23$Secondary_cluster_annotation[NMF_matrix_23$Secondary_score == 0] <- NA
NMF_matrix_23$Tertiary_cluster_annotation[NMF_matrix_23$Tertiary_score == 0] <- NA

# NMF 24
NMF_matrix_24 <- NMF_matrix_24 %>% 
  mutate(Primary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_24)-3))], 1, function(x) names(x)[maxn(1)(x)])) 
NMF_matrix_24 <- NMF_matrix_24 %>% 
  mutate(Secondary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_24)-4))], 1, function(x) names(x)[maxn(2)(x)])) 
NMF_matrix_24 <- NMF_matrix_24 %>% 
  mutate(Tertiary_cluster_annotation=apply(.[c(2:(ncol(NMF_matrix_24)-5))], 1, function(x) names(x)[maxn(3)(x)])) 

NMF_matrix_24$Secondary_cluster_annotation[NMF_matrix_24$Secondary_score == 0] <- NA
NMF_matrix_24$Tertiary_cluster_annotation[NMF_matrix_24$Tertiary_score == 0] <- NA

# Output matrices
write.table(NMF_matrix_10, "0.5_nmf_prey_matrix_results_10_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_11, "0.5_nmf_prey_matrix_results_11_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_12, "0.5_nmf_prey_matrix_results_12_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_13, "0.5_nmf_prey_matrix_results_13_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_14, "0.5_nmf_prey_matrix_results_14_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_15, "0.5_nmf_prey_matrix_results_15_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_16, "0.5_nmf_prey_matrix_results_16_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_17, "0.5_nmf_prey_matrix_results_17_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_18, "0.5_nmf_prey_matrix_results_18_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_19, "0.5_nmf_prey_matrix_results_19_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_20, "0.5_nmf_prey_matrix_results_20_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_21, "0.5_nmf_prey_matrix_results_21_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_22, "0.5_nmf_prey_matrix_results_22_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_23, "0.5_nmf_prey_matrix_results_23_top10ctls.csv", sep=",", row.names=FALSE)
write.table(NMF_matrix_24, "0.5_nmf_prey_matrix_results_24_top10ctls.csv", sep=",", row.names=FALSE)

#---------
##Part 2: Run gProfiler for each of the samples using Boris' custom terms
library("gprofiler2")
library(openxlsx)

curated_terms <- c("Nuclear speckle curated",
           "Paraspeckle curated", 
           "Comprehensive nucleolus",
           "Cajal body curated",
           "PML body curated"
)

# gProfiler NMF 10
#ranks_10 <- unique(NMF_matrix_10$Primary_cluster_annotation) 
ranks_10 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9")
df_rank_10 <- data.frame()
wb_10 = createWorkbook()

for (i in ranks_10) {
  query <- NMF_matrix_10[NMF_matrix_10$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_10, sheet_name)
  writeData(wb_10, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_10 <- rbind(df2, df_rank_10)
}

saveWorkbook(wb_10, 'gprofiler_Results_NMF_10_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_10$parents <- as.character(df_rank_10$parents)
df_rank_10$Total_ranks <- "10"
df_rank_10 <- df_rank_10[,c(15,1:13)]
df_rank_10_curated <- df_rank_10[df_rank_10$term_name %in% curated_terms,]
df_rank_10_curated$F1score <- (2*(df_rank_10_curated$precision*df_rank_10_curated$recall)) / (df_rank_10_curated$precision + df_rank_10_curated$recall)
df_rank_10_curated$F0.5score <- ((1 + 0.5^2) * df_rank_10_curated$precision * df_rank_10_curated$recall) / (0.5^2 * df_rank_10_curated$precision + df_rank_10_curated$recall)
df_rank_10_curated$F2score <- ((1 + 2^2) * df_rank_10_curated$precision * df_rank_10_curated$recall) / (2^2 * df_rank_10_curated$precision + df_rank_10_curated$recall)

write.table(df_rank_10, "1_gProfile_Results_NMF_10_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 11
#ranks_11 <- unique(NMF_matrix_11$Primary_cluster_annotation) 
ranks_11 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
df_rank_11 <- data.frame()
wb_11 = createWorkbook()

for (i in ranks_11) {
  query <- NMF_matrix_11[NMF_matrix_11$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_11, sheet_name)
  writeData(wb_11, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_11 <- rbind(df2, df_rank_11)
}

saveWorkbook(wb_11, 'gprofiler_Results_NMF_11_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_11$parents <- as.character(df_rank_11$parents)
df_rank_11$Total_ranks <- "11"
df_rank_11 <- df_rank_11[,c(15,1:13)]
df_rank_11_curated <- df_rank_11[df_rank_11$term_name %in% curated_terms,]
df_rank_11_curated$F1score <- (2*(df_rank_11_curated$precision*df_rank_11_curated$recall)) / (df_rank_11_curated$precision + df_rank_11_curated$recall)
df_rank_11_curated$F0.5score <- ((1 + 0.5^2) * df_rank_11_curated$precision * df_rank_11_curated$recall) / (0.5^2 * df_rank_11_curated$precision + df_rank_11_curated$recall)
df_rank_11_curated$F2score <- ((1 + 2^2) * df_rank_11_curated$precision * df_rank_11_curated$recall) / (2^2 * df_rank_11_curated$precision + df_rank_11_curated$recall)

write.table(df_rank_11, "1_gProfile_Results_NMF_11_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 12
#ranks_12 <- unique(NMF_matrix_12$Primary_cluster_annotation) 
ranks_12 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11")
df_rank_12 <- data.frame()
wb_12 = createWorkbook()

for (i in ranks_12) {
  query <- NMF_matrix_12[NMF_matrix_12$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_12, sheet_name)
  writeData(wb_12, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_12 <- rbind(df2, df_rank_12)
}

saveWorkbook(wb_12, 'gprofiler_Results_NMF_12_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_12$parents <- as.character(df_rank_12$parents)
df_rank_12$Total_ranks <- "12"
df_rank_12 <- df_rank_12[,c(15,1:13)]
df_rank_12_curated <- df_rank_12[df_rank_12$term_name %in% curated_terms,]
df_rank_12_curated$F1score <- (2*(df_rank_12_curated$precision*df_rank_12_curated$recall)) / (df_rank_12_curated$precision + df_rank_12_curated$recall)
df_rank_12_curated$F0.5score <- ((1 + 0.5^2) * df_rank_12_curated$precision * df_rank_12_curated$recall) / (0.5^2 * df_rank_12_curated$precision + df_rank_12_curated$recall)
df_rank_12_curated$F2score <- ((1 + 2^2) * df_rank_12_curated$precision * df_rank_12_curated$recall) / (2^2 * df_rank_12_curated$precision + df_rank_12_curated$recall)

write.table(df_rank_12, "1_gProfile_Results_NMF_12_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 13
ranks_13 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12")
df_rank_13 <- data.frame()
wb_13 = createWorkbook()

for (i in ranks_13) {
  query <- NMF_matrix_13[NMF_matrix_13$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_13, sheet_name)
  writeData(wb_13, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_13 <- rbind(df2, df_rank_13)
}

saveWorkbook(wb_13, 'gprofiler_Results_NMF_13_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_13$parents <- as.character(df_rank_13$parents)
df_rank_13$Total_ranks <- "13"
df_rank_13 <- df_rank_13[,c(15,1:13)]
df_rank_13_curated <- df_rank_13[df_rank_13$term_name %in% curated_terms,]
df_rank_13_curated$F1score <- (2*(df_rank_13_curated$precision*df_rank_13_curated$recall)) / (df_rank_13_curated$precision + df_rank_13_curated$recall)
df_rank_13_curated$F0.5score <- ((1 + 0.5^2) * df_rank_13_curated$precision * df_rank_13_curated$recall) / (0.5^2 * df_rank_13_curated$precision + df_rank_13_curated$recall)
df_rank_13_curated$F2score <- ((1 + 2^2) * df_rank_13_curated$precision * df_rank_13_curated$recall) / (2^2 * df_rank_13_curated$precision + df_rank_13_curated$recall)

write.table(df_rank_13, "1_gProfile_Results_NMF_13_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 14
ranks_14 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13")
df_rank_14 <- data.frame()
wb_14 = createWorkbook()

for (i in ranks_14) {
  query <- NMF_matrix_14[NMF_matrix_14$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_14, sheet_name)
  writeData(wb_14, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_14 <- rbind(df2, df_rank_14)
}

saveWorkbook(wb_14, 'gprofiler_Results_NMF_14_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_14$parents <- as.character(df_rank_14$parents)
df_rank_14$Total_ranks <- "14"
df_rank_14 <- df_rank_14[,c(15,1:13)]
df_rank_14_curated <- df_rank_14[df_rank_14$term_name %in% curated_terms,]
df_rank_14_curated$F1score <- (2*(df_rank_14_curated$precision*df_rank_14_curated$recall)) / (df_rank_14_curated$precision + df_rank_14_curated$recall)
df_rank_14_curated$F0.5score <- ((1 + 0.5^2) * df_rank_14_curated$precision * df_rank_14_curated$recall) / (0.5^2 * df_rank_14_curated$precision + df_rank_14_curated$recall)
df_rank_14_curated$F2score <- ((1 + 2^2) * df_rank_14_curated$precision * df_rank_14_curated$recall) / (2^2 * df_rank_14_curated$precision + df_rank_14_curated$recall)

write.table(df_rank_14, "1_gProfile_Results_NMF_14_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 15
ranks_15 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14")
df_rank_15 <- data.frame()
wb_15 = createWorkbook()

for (i in ranks_15) {
  query <- NMF_matrix_15[NMF_matrix_15$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_15, sheet_name)
  writeData(wb_15, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_15 <- rbind(df2, df_rank_15)
}

saveWorkbook(wb_15, 'gprofiler_Results_NMF_15_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_15$parents <- as.character(df_rank_15$parents)
df_rank_15$Total_ranks <- "15"
df_rank_15 <- df_rank_15[,c(15,1:13)]
df_rank_15_curated <- df_rank_15[df_rank_15$term_name %in% curated_terms,]
df_rank_15_curated$F1score <- (2*(df_rank_15_curated$precision*df_rank_15_curated$recall)) / (df_rank_15_curated$precision + df_rank_15_curated$recall)
df_rank_15_curated$F0.5score <- ((1 + 0.5^2) * df_rank_15_curated$precision * df_rank_15_curated$recall) / (0.5^2 * df_rank_15_curated$precision + df_rank_15_curated$recall)
df_rank_15_curated$F2score <- ((1 + 2^2) * df_rank_15_curated$precision * df_rank_15_curated$recall) / (2^2 * df_rank_15_curated$precision + df_rank_15_curated$recall)

write.table(df_rank_15, "1_gProfile_Results_NMF_15_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 16
ranks_16 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15")
df_rank_16 <- data.frame()
wb_16 = createWorkbook()

for (i in ranks_16) {
  query <- NMF_matrix_16[NMF_matrix_16$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_16, sheet_name)
  writeData(wb_16, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_16 <- rbind(df2, df_rank_16)
}

saveWorkbook(wb_16, 'gprofiler_Results_NMF_16_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_16$parents <- as.character(df_rank_16$parents)
df_rank_16$Total_ranks <- "16"
df_rank_16 <- df_rank_16[,c(15,1:13)]
df_rank_16_curated <- df_rank_16[df_rank_16$term_name %in% curated_terms,]
df_rank_16_curated$F1score <- (2*(df_rank_16_curated$precision*df_rank_16_curated$recall)) / (df_rank_16_curated$precision + df_rank_16_curated$recall)
df_rank_16_curated$F0.5score <- ((1 + 0.5^2) * df_rank_16_curated$precision * df_rank_16_curated$recall) / (0.5^2 * df_rank_16_curated$precision + df_rank_16_curated$recall)
df_rank_16_curated$F2score <- ((1 + 2^2) * df_rank_16_curated$precision * df_rank_16_curated$recall) / (2^2 * df_rank_16_curated$precision + df_rank_16_curated$recall)

write.table(df_rank_16, "1_gProfile_Results_NMF_16_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 17
ranks_17 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16")
df_rank_17 <- data.frame()
wb_17 = createWorkbook()

for (i in ranks_17) {
  query <- NMF_matrix_17[NMF_matrix_17$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_17, sheet_name)
  writeData(wb_17, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_17 <- rbind(df2, df_rank_17)
}

saveWorkbook(wb_17, 'gprofiler_Results_NMF_17_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_17$parents <- as.character(df_rank_17$parents)
df_rank_17$Total_ranks <- "17"
df_rank_17 <- df_rank_17[,c(15,1:13)]
df_rank_17_curated <- df_rank_17[df_rank_17$term_name %in% curated_terms,]
df_rank_17_curated$F1score <- (2*(df_rank_17_curated$precision*df_rank_17_curated$recall)) / (df_rank_17_curated$precision + df_rank_17_curated$recall)
df_rank_17_curated$F0.5score <- ((1 + 0.5^2) * df_rank_17_curated$precision * df_rank_17_curated$recall) / (0.5^2 * df_rank_17_curated$precision + df_rank_17_curated$recall)
df_rank_17_curated$F2score <- ((1 + 2^2) * df_rank_17_curated$precision * df_rank_17_curated$recall) / (2^2 * df_rank_17_curated$precision + df_rank_17_curated$recall)

write.table(df_rank_17, "1_gProfile_Results_NMF_17_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 18
ranks_18 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17")
df_rank_18 <- data.frame()
wb_18 = createWorkbook()

for (i in ranks_18) {
  query <- NMF_matrix_18[NMF_matrix_18$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_18, sheet_name)
  writeData(wb_18, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_18 <- rbind(df2, df_rank_18)
}

saveWorkbook(wb_18, 'gprofiler_Results_NMF_18_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_18$parents <- as.character(df_rank_18$parents)
df_rank_18$Total_ranks <- "18"
df_rank_18 <- df_rank_18[,c(15,1:13)]
df_rank_18_curated <- df_rank_18[df_rank_18$term_name %in% curated_terms,]
df_rank_18_curated$F1score <- (2*(df_rank_18_curated$precision*df_rank_18_curated$recall)) / (df_rank_18_curated$precision + df_rank_18_curated$recall)
df_rank_18_curated$F0.5score <- ((1 + 0.5^2) * df_rank_18_curated$precision * df_rank_18_curated$recall) / (0.5^2 * df_rank_18_curated$precision + df_rank_18_curated$recall)
df_rank_18_curated$F2score <- ((1 + 2^2) * df_rank_18_curated$precision * df_rank_18_curated$recall) / (2^2 * df_rank_18_curated$precision + df_rank_18_curated$recall)

write.table(df_rank_18, "1_gProfile_Results_NMF_18_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 19
ranks_19 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18")
df_rank_19 <- data.frame()
wb_19 = createWorkbook()

for (i in ranks_19) {
  query <- NMF_matrix_19[NMF_matrix_19$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_19, sheet_name)
  writeData(wb_19, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_19 <- rbind(df2, df_rank_19)
}

saveWorkbook(wb_19, 'gprofiler_Results_NMF_19_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_19$parents <- as.character(df_rank_19$parents)
df_rank_19$Total_ranks <- "19"
df_rank_19 <- df_rank_19[,c(15,1:13)]
df_rank_19_curated <- df_rank_19[df_rank_19$term_name %in% curated_terms,]
df_rank_19_curated$F1score <- (2*(df_rank_19_curated$precision*df_rank_19_curated$recall)) / (df_rank_19_curated$precision + df_rank_19_curated$recall)
df_rank_19_curated$F0.5score <- ((1 + 0.5^2) * df_rank_19_curated$precision * df_rank_19_curated$recall) / (0.5^2 * df_rank_19_curated$precision + df_rank_19_curated$recall)
df_rank_19_curated$F2score <- ((1 + 2^2) * df_rank_19_curated$precision * df_rank_19_curated$recall) / (2^2 * df_rank_19_curated$precision + df_rank_19_curated$recall)

write.table(df_rank_19, "1_gProfile_Results_NMF_19_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 20
ranks_20 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19")
df_rank_20 <- data.frame()
wb_20 = createWorkbook()

for (i in ranks_20) {
  query <- NMF_matrix_20[NMF_matrix_20$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_20, sheet_name)
  writeData(wb_20, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_20 <- rbind(df2, df_rank_20)
}

saveWorkbook(wb_20, 'gprofiler_Results_NMF_20_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_20$parents <- as.character(df_rank_20$parents)
df_rank_20$Total_ranks <- "20"
df_rank_20 <- df_rank_20[,c(15,1:13)]
df_rank_20_curated <- df_rank_20[df_rank_20$term_name %in% curated_terms,]
df_rank_20_curated$F1score <- (2*(df_rank_20_curated$precision*df_rank_20_curated$recall)) / (df_rank_20_curated$precision + df_rank_20_curated$recall)
df_rank_20_curated$F0.5score <- ((1 + 0.5^2) * df_rank_20_curated$precision * df_rank_20_curated$recall) / (0.5^2 * df_rank_20_curated$precision + df_rank_20_curated$recall)
df_rank_20_curated$F2score <- ((1 + 2^2) * df_rank_20_curated$precision * df_rank_20_curated$recall) / (2^2 * df_rank_20_curated$precision + df_rank_20_curated$recall)

write.table(df_rank_20, "1_gProfile_Results_NMF_20_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 21
ranks_21 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20")
df_rank_21 <- data.frame()
wb_21 = createWorkbook()

for (i in ranks_21) {
  query <- NMF_matrix_21[NMF_matrix_21$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_21, sheet_name)
  writeData(wb_21, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_21 <- rbind(df2, df_rank_21)
}

saveWorkbook(wb_21, 'gprofiler_Results_NMF_21_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_21$parents <- as.character(df_rank_21$parents)
df_rank_21$Total_ranks <- "21"
df_rank_21 <- df_rank_21[,c(15,1:13)]
df_rank_21_curated <- df_rank_21[df_rank_21$term_name %in% curated_terms,]
df_rank_21_curated$F1score <- (2*(df_rank_21_curated$precision*df_rank_21_curated$recall)) / (df_rank_21_curated$precision + df_rank_21_curated$recall)
df_rank_21_curated$F0.5score <- ((1 + 0.5^2) * df_rank_21_curated$precision * df_rank_21_curated$recall) / (0.5^2 * df_rank_21_curated$precision + df_rank_21_curated$recall)
df_rank_21_curated$F2score <- ((1 + 2^2) * df_rank_21_curated$precision * df_rank_21_curated$recall) / (2^2 * df_rank_21_curated$precision + df_rank_21_curated$recall)

write.table(df_rank_21, "1_gProfile_Results_NMF_21_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 22
ranks_22 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21")
df_rank_22 <- data.frame()
wb_22 = createWorkbook()

for (i in ranks_22) {
  query <- NMF_matrix_22[NMF_matrix_22$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_22, sheet_name)
  writeData(wb_22, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_22 <- rbind(df2, df_rank_22)
}

saveWorkbook(wb_22, 'gprofiler_Results_NMF_22_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_22$parents <- as.character(df_rank_22$parents)
df_rank_22$Total_ranks <- "22"
df_rank_22 <- df_rank_22[,c(15,1:13)]
df_rank_22_curated <- df_rank_22[df_rank_22$term_name %in% curated_terms,]
df_rank_22_curated$F1score <- (2*(df_rank_22_curated$precision*df_rank_22_curated$recall)) / (df_rank_22_curated$precision + df_rank_22_curated$recall)
df_rank_22_curated$F0.5score <- ((1 + 0.5^2) * df_rank_22_curated$precision * df_rank_22_curated$recall) / (0.5^2 * df_rank_22_curated$precision + df_rank_22_curated$recall)
df_rank_22_curated$F2score <- ((1 + 2^2) * df_rank_22_curated$precision * df_rank_22_curated$recall) / (2^2 * df_rank_22_curated$precision + df_rank_22_curated$recall)

write.table(df_rank_22, "1_gProfile_Results_NMF_22_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 23
ranks_23 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22")
df_rank_23 <- data.frame()
wb_23 = createWorkbook()

for (i in ranks_23) {
  query <- NMF_matrix_23[NMF_matrix_23$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_23, sheet_name)
  writeData(wb_23, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_23 <- rbind(df2, df_rank_23)
}

saveWorkbook(wb_23, 'gprofiler_Results_NMF_23_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_23$parents <- as.character(df_rank_23$parents)
df_rank_23$Total_ranks <- "23"
df_rank_23 <- df_rank_23[,c(15,1:13)]
df_rank_23_curated <- df_rank_23[df_rank_23$term_name %in% curated_terms,]
df_rank_23_curated$F1score <- (2*(df_rank_23_curated$precision*df_rank_23_curated$recall)) / (df_rank_23_curated$precision + df_rank_23_curated$recall)
df_rank_23_curated$F0.5score <- ((1 + 0.5^2) * df_rank_23_curated$precision * df_rank_23_curated$recall) / (0.5^2 * df_rank_23_curated$precision + df_rank_23_curated$recall)
df_rank_23_curated$F2score <- ((1 + 2^2) * df_rank_23_curated$precision * df_rank_23_curated$recall) / (2^2 * df_rank_23_curated$precision + df_rank_23_curated$recall)

write.table(df_rank_23, "1_gProfile_Results_NMF_23_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)


# gProfiler NMF 24
ranks_24 <- c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23")
df_rank_24 <- data.frame()
wb_24 = createWorkbook()

for (i in ranks_24) {
  query <- NMF_matrix_24[NMF_matrix_24$Primary_cluster_annotation == i,]
  gene <- query$PreyGene
  gprofiler <- gost(query = gene,
                    organism = "gp__2aO3_xfUn_hEA", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "gSCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC"), as_short_link = FALSE, highlight = TRUE)
  df <- data.frame(gprofiler$result)
  colnames(df)[3] <- "adjusted_p_value"
  df$query <- i
  df <- df[order(df$adjusted_p_value), ]
  df$parents <- lapply(df$parents, function(x){paste0(x,collapse='|')})
  sheet_name = paste(i)
  addWorksheet(wb_24, sheet_name)
  writeData(wb_24, sheet_name, df)
  df2 <- df[df$source == "hsapiens.GO:MF.name" | df$source == "hsapiens.GO:BP.name" | df$source == "hsapiens.GO:CC.name" | df$source == "Dyakov_Curated_NB_lists_20231201"| df$source == "gprofiler_full_hsapiens.name_PLUS_DYAKOV_20231211",]
  df_rank_24 <- rbind(df2, df_rank_24)
}

saveWorkbook(wb_24, 'gprofiler_Results_NMF_24_ranks_customToken_top10ctls.xlsx', overwrite = TRUE)
df_rank_24$parents <- as.character(df_rank_24$parents)
df_rank_24$Total_ranks <- "24"
df_rank_24 <- df_rank_24[,c(15,1:13)]
df_rank_24_curated <- df_rank_24[df_rank_24$term_name %in% curated_terms,]
df_rank_24_curated$F1score <- (2*(df_rank_24_curated$precision*df_rank_24_curated$recall)) / (df_rank_24_curated$precision + df_rank_24_curated$recall)
df_rank_24_curated$F0.5score <- ((1 + 0.5^2) * df_rank_24_curated$precision * df_rank_24_curated$recall) / (0.5^2 * df_rank_24_curated$precision + df_rank_24_curated$recall)
df_rank_24_curated$F2score <- ((1 + 2^2) * df_rank_24_curated$precision * df_rank_24_curated$recall) / (2^2 * df_rank_24_curated$precision + df_rank_24_curated$recall)

write.table(df_rank_24, "1_gProfile_Results_NMF_24_ranks_customToken_top10ctls.txt", sep="\t", row.names=FALSE)

# Combine all the gProfiler curated results into one
df_combined_curated <- rbind(df_rank_10_curated,
                             df_rank_11_curated,
                             df_rank_12_curated,
                             df_rank_13_curated,
                             df_rank_14_curated,
                             df_rank_15_curated,
                             df_rank_16_curated,
                             df_rank_17_curated,
                             df_rank_18_curated,
                             df_rank_19_curated,
                             df_rank_20_curated,
                             df_rank_21_curated,
                             df_rank_22_curated,
                             df_rank_23_curated,
                             df_rank_24_curated
                             )

write.table(df_combined_curated, "2_gProfile_Results_NMF_curated_terms_combined_top10ctls.txt", sep="\t", row.names=FALSE)

#---------
## Part 3: Calculate precision/recall and enrichment and plot across clusters
# PS stats (for top cluster)
# Create a table with the total number and 
PS_10 <- df_rank_10_curated[df_rank_10_curated$term_name == "Paraspeckle curated",]
PS_10$negative_log10_pvalue <- -log10(PS_10$adjusted_p_value)
PS_10$Number_clusters_with_term <- nrow(PS_10[PS_10$significant == TRUE,])

PS_11 <- df_rank_11_curated[df_rank_11_curated$term_name == "Paraspeckle curated",]
PS_11$negative_log10_pvalue <- -log10(PS_11$adjusted_p_value)
PS_11$Number_clusters_with_term <- nrow(PS_11[PS_11$significant == TRUE,])

PS_12 <- df_rank_12_curated[df_rank_12_curated$term_name == "Paraspeckle curated",]
PS_12$negative_log10_pvalue <- -log10(PS_12$adjusted_p_value)
PS_12$Number_clusters_with_term <- nrow(PS_12[PS_12$significant == TRUE,])

PS_13 <- df_rank_13_curated[df_rank_13_curated$term_name == "Paraspeckle curated",]
PS_13$negative_log10_pvalue <- -log10(PS_13$adjusted_p_value)
PS_13$Number_clusters_with_term <- nrow(PS_13[PS_13$significant == TRUE,])

PS_14 <- df_rank_14_curated[df_rank_14_curated$term_name == "Paraspeckle curated",]
PS_14$negative_log10_pvalue <- -log10(PS_14$adjusted_p_value)
PS_14$Number_clusters_with_term <- nrow(PS_14[PS_14$significant == TRUE,])

PS_15 <- df_rank_15_curated[df_rank_15_curated$term_name == "Paraspeckle curated",]
PS_15$negative_log10_pvalue <- -log10(PS_15$adjusted_p_value)
PS_15$Number_clusters_with_term <- nrow(PS_15[PS_15$significant == TRUE,])

PS_16 <- df_rank_16_curated[df_rank_16_curated$term_name == "Paraspeckle curated",]
PS_16$negative_log10_pvalue <- -log10(PS_16$adjusted_p_value)
PS_16$Number_clusters_with_term <- nrow(PS_16[PS_16$significant == TRUE,])

PS_17 <- df_rank_17_curated[df_rank_17_curated$term_name == "Paraspeckle curated",]
PS_17$negative_log10_pvalue <- -log10(PS_17$adjusted_p_value)
PS_17$Number_clusters_with_term <- nrow(PS_17[PS_17$significant == TRUE,])

PS_18 <- df_rank_18_curated[df_rank_18_curated$term_name == "Paraspeckle curated",]
PS_18$negative_log10_pvalue <- -log10(PS_18$adjusted_p_value)
PS_18$Number_clusters_with_term <- nrow(PS_18[PS_18$significant == TRUE,])

PS_19 <- df_rank_19_curated[df_rank_19_curated$term_name == "Paraspeckle curated",]
PS_19$negative_log10_pvalue <- -log10(PS_19$adjusted_p_value)
PS_19$Number_clusters_with_term <- nrow(PS_19[PS_19$significant == TRUE,])

PS_20 <- df_rank_20_curated[df_rank_20_curated$term_name == "Paraspeckle curated",]
PS_20$negative_log10_pvalue <- -log10(PS_20$adjusted_p_value)
PS_20$Number_clusters_with_term <- nrow(PS_20[PS_20$significant == TRUE,])

PS_21 <- df_rank_21_curated[df_rank_21_curated$term_name == "Paraspeckle curated",]
PS_21$negative_log10_pvalue <- -log10(PS_21$adjusted_p_value)
PS_21$Number_clusters_with_term <- nrow(PS_21[PS_21$significant == TRUE,])

PS_22 <- df_rank_22_curated[df_rank_22_curated$term_name == "Paraspeckle curated",]
PS_22$negative_log10_pvalue <- -log10(PS_22$adjusted_p_value)
PS_22$Number_clusters_with_term <- nrow(PS_22[PS_22$significant == TRUE,])

PS_23 <- df_rank_23_curated[df_rank_23_curated$term_name == "Paraspeckle curated",]
PS_23$negative_log10_pvalue <- -log10(PS_23$adjusted_p_value)
PS_23$Number_clusters_with_term <- nrow(PS_23[PS_23$significant == TRUE,])

PS_24 <- df_rank_24_curated[df_rank_24_curated$term_name == "Paraspeckle curated",]
PS_24$negative_log10_pvalue <- -log10(PS_24$adjusted_p_value)
PS_24$Number_clusters_with_term <- nrow(PS_24[PS_24$significant == TRUE,])

top_rank_PS <- rbind(PS_10[which.max(PS_10$negative_log10_pvalue),],
                     PS_11[which.max(PS_11$negative_log10_pvalue),],
                     PS_12[which.max(PS_12$negative_log10_pvalue),],
                     PS_13[which.max(PS_13$negative_log10_pvalue),],
                     PS_14[which.max(PS_14$negative_log10_pvalue),],
                     PS_15[which.max(PS_15$negative_log10_pvalue),],
                     PS_16[which.max(PS_16$negative_log10_pvalue),],
                     PS_17[which.max(PS_17$negative_log10_pvalue),],
                     PS_18[which.max(PS_18$negative_log10_pvalue),],
                     PS_19[which.max(PS_19$negative_log10_pvalue),],
                     PS_20[which.max(PS_20$negative_log10_pvalue),],
                     PS_21[which.max(PS_21$negative_log10_pvalue),],
                     PS_22[which.max(PS_22$negative_log10_pvalue),],
                     PS_23[which.max(PS_23$negative_log10_pvalue),],
                     PS_24[which.max(PS_24$negative_log10_pvalue),]
                     )

# NS stats (for top cluster)
# Create a table with the total number and 
NS_10 <- df_rank_10_curated[df_rank_10_curated$term_name == "Nuclear speckle curated",]
NS_10$negative_log10_pvalue <- -log10(NS_10$adjusted_p_value)
NS_10$Number_clusters_with_term <- nrow(NS_10[NS_10$significant == TRUE,])

NS_11 <- df_rank_11_curated[df_rank_11_curated$term_name == "Nuclear speckle curated",]
NS_11$negative_log10_pvalue <- -log10(NS_11$adjusted_p_value)
NS_11$Number_clusters_with_term <- nrow(NS_11[NS_11$significant == TRUE,])

NS_12 <- df_rank_12_curated[df_rank_12_curated$term_name == "Nuclear speckle curated",]
NS_12$negative_log10_pvalue <- -log10(NS_12$adjusted_p_value)
NS_12$Number_clusters_with_term <- nrow(NS_12[NS_12$significant == TRUE,])

NS_13 <- df_rank_13_curated[df_rank_13_curated$term_name == "Nuclear speckle curated",]
NS_13$negative_log10_pvalue <- -log10(NS_13$adjusted_p_value)
NS_13$Number_clusters_with_term <- nrow(NS_13[NS_13$significant == TRUE,])

NS_14 <- df_rank_14_curated[df_rank_14_curated$term_name == "Nuclear speckle curated",]
NS_14$negative_log10_pvalue <- -log10(NS_14$adjusted_p_value)
NS_14$Number_clusters_with_term <- nrow(NS_14[NS_14$significant == TRUE,])

NS_15 <- df_rank_15_curated[df_rank_15_curated$term_name == "Nuclear speckle curated",]
NS_15$negative_log10_pvalue <- -log10(NS_15$adjusted_p_value)
NS_15$Number_clusters_with_term <- nrow(NS_15[NS_15$significant == TRUE,])

NS_16 <- df_rank_16_curated[df_rank_16_curated$term_name == "Nuclear speckle curated",]
NS_16$negative_log10_pvalue <- -log10(NS_16$adjusted_p_value)
NS_16$Number_clusters_with_term <- nrow(NS_16[NS_16$significant == TRUE,])

NS_17 <- df_rank_17_curated[df_rank_17_curated$term_name == "Nuclear speckle curated",]
NS_17$negative_log10_pvalue <- -log10(NS_17$adjusted_p_value)
NS_17$Number_clusters_with_term <- nrow(NS_17[NS_17$significant == TRUE,])

NS_18 <- df_rank_18_curated[df_rank_18_curated$term_name == "Nuclear speckle curated",]
NS_18$negative_log10_pvalue <- -log10(NS_18$adjusted_p_value)
NS_18$Number_clusters_with_term <- nrow(NS_18[NS_18$significant == TRUE,])

NS_19 <- df_rank_19_curated[df_rank_19_curated$term_name == "Nuclear speckle curated",]
NS_19$negative_log10_pvalue <- -log10(NS_19$adjusted_p_value)
NS_19$Number_clusters_with_term <- nrow(NS_19[NS_19$significant == TRUE,])

NS_20 <- df_rank_20_curated[df_rank_20_curated$term_name == "Nuclear speckle curated",]
NS_20$negative_log10_pvalue <- -log10(NS_20$adjusted_p_value)
NS_20$Number_clusters_with_term <- nrow(NS_20[NS_20$significant == TRUE,])

NS_21 <- df_rank_21_curated[df_rank_21_curated$term_name == "Nuclear speckle curated",]
NS_21$negative_log10_pvalue <- -log10(NS_21$adjusted_p_value)
NS_21$Number_clusters_with_term <- nrow(NS_21[NS_21$significant == TRUE,])

NS_22 <- df_rank_22_curated[df_rank_22_curated$term_name == "Nuclear speckle curated",]
NS_22$negative_log10_pvalue <- -log10(NS_22$adjusted_p_value)
NS_22$Number_clusters_with_term <- nrow(NS_22[NS_22$significant == TRUE,])

NS_23 <- df_rank_23_curated[df_rank_23_curated$term_name == "Nuclear speckle curated",]
NS_23$negative_log10_pvalue <- -log10(NS_23$adjusted_p_value)
NS_23$Number_clusters_with_term <- nrow(NS_23[NS_23$significant == TRUE,])

NS_24 <- df_rank_24_curated[df_rank_24_curated$term_name == "Nuclear speckle curated",]
NS_24$negative_log10_pvalue <- -log10(NS_24$adjusted_p_value)
NS_24$Number_clusters_with_term <- nrow(NS_24[NS_24$significant == TRUE,])

top_rank_NS <- rbind(NS_10[which.max(NS_10$negative_log10_pvalue),],
                     NS_11[which.max(NS_11$negative_log10_pvalue),],
                     NS_12[which.max(NS_12$negative_log10_pvalue),],
                     NS_13[which.max(NS_13$negative_log10_pvalue),],
                     NS_14[which.max(NS_14$negative_log10_pvalue),],
                     NS_15[which.max(NS_15$negative_log10_pvalue),],
                     NS_16[which.max(NS_16$negative_log10_pvalue),],
                     NS_17[which.max(NS_17$negative_log10_pvalue),],
                     NS_18[which.max(NS_18$negative_log10_pvalue),],
                     NS_19[which.max(NS_19$negative_log10_pvalue),],
                     NS_20[which.max(NS_20$negative_log10_pvalue),],
                     NS_21[which.max(NS_21$negative_log10_pvalue),],
                     NS_22[which.max(NS_22$negative_log10_pvalue),],
                     NS_23[which.max(NS_23$negative_log10_pvalue),],
                     NS_24[which.max(NS_24$negative_log10_pvalue),]
)

# Nucleolus stats (for top cluster)
# Create a table with the total number and 
NUC_10 <- df_rank_10_curated[df_rank_10_curated$term_name == "Comprehensive nucleolus",]
NUC_10$negative_log10_pvalue <- -log10(NUC_10$adjusted_p_value)
NUC_10$Number_clusters_with_term <- nrow(NUC_10[NUC_10$significant == TRUE,])

NUC_11 <- df_rank_11_curated[df_rank_11_curated$term_name == "Comprehensive nucleolus",]
NUC_11$negative_log10_pvalue <- -log10(NUC_11$adjusted_p_value)
NUC_11$Number_clusters_with_term <- nrow(NUC_11[NUC_11$significant == TRUE,])

NUC_12 <- df_rank_12_curated[df_rank_12_curated$term_name == "Comprehensive nucleolus",]
NUC_12$negative_log10_pvalue <- -log10(NUC_12$adjusted_p_value)
NUC_12$Number_clusters_with_term <- nrow(NUC_12[NUC_12$significant == TRUE,])

NUC_13 <- df_rank_13_curated[df_rank_13_curated$term_name == "Comprehensive nucleolus",]
NUC_13$negative_log10_pvalue <- -log10(NUC_13$adjusted_p_value)
NUC_13$Number_clusters_with_term <- nrow(NUC_13[NUC_13$significant == TRUE,])

NUC_14 <- df_rank_14_curated[df_rank_14_curated$term_name == "Comprehensive nucleolus",]
NUC_14$negative_log10_pvalue <- -log10(NUC_14$adjusted_p_value)
NUC_14$Number_clusters_with_term <- nrow(NUC_14[NUC_14$significant == TRUE,])

NUC_15 <- df_rank_15_curated[df_rank_15_curated$term_name == "Comprehensive nucleolus",]
NUC_15$negative_log10_pvalue <- -log10(NUC_15$adjusted_p_value)
NUC_15$Number_clusters_with_term <- nrow(NUC_15[NUC_15$significant == TRUE,])

NUC_16 <- df_rank_16_curated[df_rank_16_curated$term_name == "Comprehensive nucleolus",]
NUC_16$negative_log10_pvalue <- -log10(NUC_16$adjusted_p_value)
NUC_16$Number_clusters_with_term <- nrow(NUC_16[NUC_16$significant == TRUE,])

NUC_17 <- df_rank_17_curated[df_rank_17_curated$term_name == "Comprehensive nucleolus",]
NUC_17$negative_log10_pvalue <- -log10(NUC_17$adjusted_p_value)
NUC_17$Number_clusters_with_term <- nrow(NUC_17[NUC_17$significant == TRUE,])

NUC_18 <- df_rank_18_curated[df_rank_18_curated$term_name == "Comprehensive nucleolus",]
NUC_18$negative_log10_pvalue <- -log10(NUC_18$adjusted_p_value)
NUC_18$Number_clusters_with_term <- nrow(NUC_18[NUC_18$significant == TRUE,])

NUC_19 <- df_rank_19_curated[df_rank_19_curated$term_name == "Comprehensive nucleolus",]
NUC_19$negative_log10_pvalue <- -log10(NUC_19$adjusted_p_value)
NUC_19$Number_clusters_with_term <- nrow(NUC_19[NUC_19$significant == TRUE,])

NUC_20 <- df_rank_20_curated[df_rank_20_curated$term_name == "Comprehensive nucleolus",]
NUC_20$negative_log10_pvalue <- -log10(NUC_20$adjusted_p_value)
NUC_20$Number_clusters_with_term <- nrow(NUC_20[NUC_20$significant == TRUE,])

NUC_21 <- df_rank_21_curated[df_rank_21_curated$term_name == "Comprehensive nucleolus",]
NUC_21$negative_log10_pvalue <- -log10(NUC_21$adjusted_p_value)
NUC_21$Number_clusters_with_term <- nrow(NUC_21[NUC_21$significant == TRUE,])

NUC_22 <- df_rank_22_curated[df_rank_22_curated$term_name == "Comprehensive nucleolus",]
NUC_22$negative_log10_pvalue <- -log10(NUC_22$adjusted_p_value)
NUC_22$Number_clusters_with_term <- nrow(NUC_22[NUC_22$significant == TRUE,])

NUC_23 <- df_rank_23_curated[df_rank_23_curated$term_name == "Comprehensive nucleolus",]
NUC_23$negative_log10_pvalue <- -log10(NUC_23$adjusted_p_value)
NUC_23$Number_clusters_with_term <- nrow(NUC_23[NUC_23$significant == TRUE,])

NUC_24 <- df_rank_24_curated[df_rank_24_curated$term_name == "Comprehensive nucleolus",]
NUC_24$negative_log10_pvalue <- -log10(NUC_24$adjusted_p_value)
NUC_24$Number_clusters_with_term <- nrow(NUC_24[NUC_24$significant == TRUE,])

top_rank_NUC <- rbind(NUC_10[which.max(NUC_10$negative_log10_pvalue),],
                     NUC_11[which.max(NUC_11$negative_log10_pvalue),],
                     NUC_12[which.max(NUC_12$negative_log10_pvalue),],
                     NUC_13[which.max(NUC_13$negative_log10_pvalue),],
                     NUC_14[which.max(NUC_14$negative_log10_pvalue),],
                     NUC_15[which.max(NUC_15$negative_log10_pvalue),],
                     NUC_16[which.max(NUC_16$negative_log10_pvalue),],
                     NUC_17[which.max(NUC_17$negative_log10_pvalue),],
                     NUC_18[which.max(NUC_18$negative_log10_pvalue),],
                     NUC_19[which.max(NUC_19$negative_log10_pvalue),],
                     NUC_20[which.max(NUC_20$negative_log10_pvalue),],
                     NUC_21[which.max(NUC_21$negative_log10_pvalue),],
                     NUC_22[which.max(NUC_22$negative_log10_pvalue),],
                     NUC_23[which.max(NUC_23$negative_log10_pvalue),],
                     NUC_24[which.max(NUC_24$negative_log10_pvalue),]
)


# CB stats (for top cluster)
# Create a table with the total number and 
CB_10 <- df_rank_10_curated[df_rank_10_curated$term_name == "Cajal body curated",]
CB_10$negative_log10_pvalue <- -log10(CB_10$adjusted_p_value)
CB_10$Number_clusters_with_term <- nrow(CB_10[CB_10$significant == TRUE,])

CB_11 <- df_rank_11_curated[df_rank_11_curated$term_name == "Cajal body curated",]
CB_11$negative_log10_pvalue <- -log10(CB_11$adjusted_p_value)
CB_11$Number_clusters_with_term <- nrow(CB_11[CB_11$significant == TRUE,])

CB_12 <- df_rank_12_curated[df_rank_12_curated$term_name == "Cajal body curated",]
CB_12$negative_log10_pvalue <- -log10(CB_12$adjusted_p_value)
CB_12$Number_clusters_with_term <- nrow(CB_12[CB_12$significant == TRUE,])

CB_13 <- df_rank_13_curated[df_rank_13_curated$term_name == "Cajal body curated",]
CB_13$negative_log10_pvalue <- -log10(CB_13$adjusted_p_value)
CB_13$Number_clusters_with_term <- nrow(CB_13[CB_13$significant == TRUE,])

CB_14 <- df_rank_14_curated[df_rank_14_curated$term_name == "Cajal body curated",]
CB_14$negative_log10_pvalue <- -log10(CB_14$adjusted_p_value)
CB_14$Number_clusters_with_term <- nrow(CB_14[CB_14$significant == TRUE,])

CB_15 <- df_rank_15_curated[df_rank_15_curated$term_name == "Cajal body curated",]
CB_15$negative_log10_pvalue <- -log10(CB_15$adjusted_p_value)
CB_15$Number_clusters_with_term <- nrow(CB_15[CB_15$significant == TRUE,])

CB_16 <- df_rank_16_curated[df_rank_16_curated$term_name == "Cajal body curated",]
CB_16$negative_log10_pvalue <- -log10(CB_16$adjusted_p_value)
CB_16$Number_clusters_with_term <- nrow(CB_16[CB_16$significant == TRUE,])

CB_17 <- df_rank_17_curated[df_rank_17_curated$term_name == "Cajal body curated",]
CB_17$negative_log10_pvalue <- -log10(CB_17$adjusted_p_value)
CB_17$Number_clusters_with_term <- nrow(CB_17[CB_17$significant == TRUE,])

CB_18 <- df_rank_18_curated[df_rank_18_curated$term_name == "Cajal body curated",]
CB_18$negative_log10_pvalue <- -log10(CB_18$adjusted_p_value)
CB_18$Number_clusters_with_term <- nrow(CB_18[CB_18$significant == TRUE,])

CB_19 <- df_rank_19_curated[df_rank_19_curated$term_name == "Cajal body curated",]
CB_19$negative_log10_pvalue <- -log10(CB_19$adjusted_p_value)
CB_19$Number_clusters_with_term <- nrow(CB_19[CB_19$significant == TRUE,])

CB_20 <- df_rank_20_curated[df_rank_20_curated$term_name == "Cajal body curated",]
CB_20$negative_log10_pvalue <- -log10(CB_20$adjusted_p_value)
CB_20$Number_clusters_with_term <- nrow(CB_20[CB_20$significant == TRUE,])

CB_21 <- df_rank_21_curated[df_rank_21_curated$term_name == "Cajal body curated",]
CB_21$negative_log10_pvalue <- -log10(CB_21$adjusted_p_value)
CB_21$Number_clusters_with_term <- nrow(CB_21[CB_21$significant == TRUE,])

CB_22 <- df_rank_22_curated[df_rank_22_curated$term_name == "Cajal body curated",]
CB_22$negative_log10_pvalue <- -log10(CB_22$adjusted_p_value)
CB_22$Number_clusters_with_term <- nrow(CB_22[CB_22$significant == TRUE,])

CB_23 <- df_rank_23_curated[df_rank_23_curated$term_name == "Cajal body curated",]
CB_23$negative_log10_pvalue <- -log10(CB_23$adjusted_p_value)
CB_23$Number_clusters_with_term <- nrow(CB_23[CB_23$significant == TRUE,])

CB_24 <- df_rank_24_curated[df_rank_24_curated$term_name == "Cajal body curated",]
CB_24$negative_log10_pvalue <- -log10(CB_24$adjusted_p_value)
CB_24$Number_clusters_with_term <- nrow(CB_24[CB_24$significant == TRUE,])

top_rank_CB <- rbind(CB_10[which.max(CB_10$negative_log10_pvalue),],
                     CB_11[which.max(CB_11$negative_log10_pvalue),],
                     CB_12[which.max(CB_12$negative_log10_pvalue),],
                     CB_13[which.max(CB_13$negative_log10_pvalue),],
                     CB_14[which.max(CB_14$negative_log10_pvalue),],
                     CB_15[which.max(CB_15$negative_log10_pvalue),],
                     CB_16[which.max(CB_16$negative_log10_pvalue),],
                     CB_17[which.max(CB_17$negative_log10_pvalue),],
                     CB_18[which.max(CB_18$negative_log10_pvalue),],
                     CB_19[which.max(CB_19$negative_log10_pvalue),],
                     CB_20[which.max(CB_20$negative_log10_pvalue),],
                     CB_21[which.max(CB_21$negative_log10_pvalue),],
                     CB_22[which.max(CB_22$negative_log10_pvalue),],
                     CB_23[which.max(CB_23$negative_log10_pvalue),],
                     CB_24[which.max(CB_24$negative_log10_pvalue),]
)


# PML stats (for top cluster)
# Create a table with the total number and 
PML_10 <- df_rank_10_curated[df_rank_10_curated$term_name == "PML body curated",]
PML_10$negative_log10_pvalue <- -log10(PML_10$adjusted_p_value)
PML_10$Number_clusters_with_term <- nrow(PML_10[PML_10$significant == TRUE,])

PML_11 <- df_rank_11_curated[df_rank_11_curated$term_name == "PML body curated",]
PML_11$negative_log10_pvalue <- -log10(PML_11$adjusted_p_value)
PML_11$Number_clusters_with_term <- nrow(PML_11[PML_11$significant == TRUE,])

PML_12 <- df_rank_12_curated[df_rank_12_curated$term_name == "PML body curated",]
PML_12$negative_log10_pvalue <- -log10(PML_12$adjusted_p_value)
PML_12$Number_clusters_with_term <- nrow(PML_12[PML_12$significant == TRUE,])

PML_13 <- df_rank_13_curated[df_rank_13_curated$term_name == "PML body curated",]
PML_13$negative_log10_pvalue <- -log10(PML_13$adjusted_p_value)
PML_13$Number_clusters_with_term <- nrow(PML_13[PML_13$significant == TRUE,])

PML_14 <- df_rank_14_curated[df_rank_14_curated$term_name == "PML body curated",]
PML_14$negative_log10_pvalue <- -log10(PML_14$adjusted_p_value)
PML_14$Number_clusters_with_term <- nrow(PML_14[PML_14$significant == TRUE,])

PML_15 <- df_rank_15_curated[df_rank_15_curated$term_name == "PML body curated",]
PML_15$negative_log10_pvalue <- -log10(PML_15$adjusted_p_value)
PML_15$Number_clusters_with_term <- nrow(PML_15[PML_15$significant == TRUE,])

PML_16 <- df_rank_16_curated[df_rank_16_curated$term_name == "PML body curated",]
PML_16$negative_log10_pvalue <- -log10(PML_16$adjusted_p_value)
PML_16$Number_clusters_with_term <- nrow(PML_16[PML_16$significant == TRUE,])

PML_17 <- df_rank_17_curated[df_rank_17_curated$term_name == "PML body curated",]
PML_17$negative_log10_pvalue <- -log10(PML_17$adjusted_p_value)
PML_17$Number_clusters_with_term <- nrow(PML_17[PML_17$significant == TRUE,])

PML_18 <- df_rank_18_curated[df_rank_18_curated$term_name == "PML body curated",]
PML_18$negative_log10_pvalue <- -log10(PML_18$adjusted_p_value)
PML_18$Number_clusters_with_term <- nrow(PML_18[PML_18$significant == TRUE,])

PML_19 <- df_rank_19_curated[df_rank_19_curated$term_name == "PML body curated",]
PML_19$negative_log10_pvalue <- -log10(PML_19$adjusted_p_value)
PML_19$Number_clusters_with_term <- nrow(PML_19[PML_19$significant == TRUE,])

PML_20 <- df_rank_20_curated[df_rank_20_curated$term_name == "PML body curated",]
PML_20$negative_log10_pvalue <- -log10(PML_20$adjusted_p_value)
PML_20$Number_clusters_with_term <- nrow(PML_20[PML_20$significant == TRUE,])

PML_21 <- df_rank_21_curated[df_rank_21_curated$term_name == "PML body curated",]
PML_21$negative_log10_pvalue <- -log10(PML_21$adjusted_p_value)
PML_21$Number_clusters_with_term <- nrow(PML_21[PML_21$significant == TRUE,])

PML_22 <- df_rank_22_curated[df_rank_22_curated$term_name == "PML body curated",]
PML_22$negative_log10_pvalue <- -log10(PML_22$adjusted_p_value)
PML_22$Number_clusters_with_term <- nrow(PML_22[PML_22$significant == TRUE,])

PML_23 <- df_rank_23_curated[df_rank_23_curated$term_name == "PML body curated",]
PML_23$negative_log10_pvalue <- -log10(PML_23$adjusted_p_value)
PML_23$Number_clusters_with_term <- nrow(PML_23[PML_23$significant == TRUE,])

PML_24 <- df_rank_24_curated[df_rank_24_curated$term_name == "PML body curated",]
PML_24$negative_log10_pvalue <- -log10(PML_24$adjusted_p_value)
PML_24$Number_clusters_with_term <- nrow(PML_24[PML_24$significant == TRUE,])

top_rank_PML <- rbind(PML_10[which.max(PML_10$negative_log10_pvalue),],
                     PML_11[which.max(PML_11$negative_log10_pvalue),],
                     PML_12[which.max(PML_12$negative_log10_pvalue),],
                     PML_13[which.max(PML_13$negative_log10_pvalue),],
                     PML_14[which.max(PML_14$negative_log10_pvalue),],
                     PML_15[which.max(PML_15$negative_log10_pvalue),],
                     PML_16[which.max(PML_16$negative_log10_pvalue),],
                     PML_17[which.max(PML_17$negative_log10_pvalue),],
                     PML_18[which.max(PML_18$negative_log10_pvalue),],
                     PML_19[which.max(PML_19$negative_log10_pvalue),],
                     PML_20[which.max(PML_20$negative_log10_pvalue),],
                     PML_21[which.max(PML_21$negative_log10_pvalue),],
                     PML_22[which.max(PML_22$negative_log10_pvalue),],
                     PML_23[which.max(PML_23$negative_log10_pvalue),],
                     PML_24[which.max(PML_24$negative_log10_pvalue),]
)


# Combine all top rank data into one file
df_combined_curated_top <- rbind(top_rank_PS,
                                 top_rank_NS,
                                 top_rank_NUC,
                                 top_rank_CB,
                                 top_rank_PML
                                 )

write.table(df_combined_curated_top, "3_gProfile_Results_NMF_curated_terms_combined_top10ctls_toprank.txt", sep="\t", row.names=FALSE)

#--------
# Plot PS data
library(ggplot2)
library(ggrepel)

#Negative log10 p-values
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_pValue_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=negative_log10_pvalue, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Negative log10 pvalues) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("-log10(pvalue) (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,125))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Recall
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_Recall_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=recall, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Recall) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("recall (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,0.5))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Precision
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_Precision_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=precision, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Precision) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("precision (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,1))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Number of ranks with term
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_NumberRankswTerm_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=Number_clusters_with_term, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Number ranks with sign term) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("Number of ranks with term (label = intersection)")+
  #xlim(c(0,1))+
  #ylim(c(0,1))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# Query size
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_Querysize_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=query_size, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_label_repel(aes(label=intersection_size), nudge_x = 0, size = 2, max.overlaps=100)+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Query Size) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("Query size (label = intersection)")+
  #xlim(c(0,1))+
  ylim(c(0,250))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F1scores
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_F1score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=F1score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:F1-score)")+
  xlab("Total Number of Ranks")+
  ylab("F1-score (2*Prec*Recall)/(Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F0.5 scores
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_F0.5score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=F0.5score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:F0.5-score)")+
  xlab("Total Number of Ranks")+
  ylab("F0.5-score ((1+0.5^2)*Prec*Recall)/(0.5^2*Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

# F2scores
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_F2score_top10ctls.pdf", width = 8, height = 6)
ggplot(data=df_combined_curated_top, aes(x=Total_ranks, y=F2score, group=1)) +
  geom_line()+
  #scale_color_manual(values=c("blue", "red"))+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:F2-score)")+
  xlab("Total Number of Ranks")+
  ylab("F2-score ((1+2^2)*Prec*Recall)/(2^2*Prec+Recall) (label = intersection) Top 10 ctl subtraction")+
  #xlim(c(0,1))+
  ylim(c(0,0.6))+
  theme_bw()+
  facet_wrap(~term_id,scales = "free")
#coord_fixed()
dev.off()

#F1score combined
F1scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Total_ranks)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Total_ranks == i,]
  df_fscores$average_F1score <- mean(df_fscores$F1score)
  avgFscore <- df_fscores[1,c(1,20)]
  F1scores <- rbind(avgFscore, F1scores)
}

# F1score average
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_average_F1score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F1scores, aes(x=Total_ranks, y=average_F1score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("Average F1score (2*Prec*Recall)/(Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

#F0.5score combined
F0.5scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Total_ranks)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Total_ranks == i,]
  df_fscores$average_F0.5score <- mean(df_fscores$F0.5score)
  avgFscore <- df_fscores[1,c(1,20)]
  F0.5scores <- rbind(avgFscore, F0.5scores)
}

# F0.5score average
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_average_F0.5score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F0.5scores, aes(x=Total_ranks, y=average_F0.5score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("Average F2score F0.5-score ((1+0.5^2)*Prec*Recall)/(0.5^2*Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

#F2score combined
F2scores <- data.frame()
num_ranks <- unique(df_combined_curated_top$Total_ranks)

for (i in num_ranks) {
  df_fscores <- df_combined_curated_top[df_combined_curated_top$Total_ranks == i,]
  df_fscores$average_F2score <- mean(df_fscores$F2score)
  avgFscore <- df_fscores[1,c(1,20)]
  F2scores <- rbind(avgFscore, F2scores)
}

# F2score average
pdf("Rplot_NMF_NumRank_testing_gProfiler_TopPvalueRank_average_F2score_top10ctls.pdf", width = 7, height = 6)
ggplot(data=F2scores, aes(x=Total_ranks, y=average_F2score, group=1)) +
  geom_line()+
  geom_point()+
  labs(title="NMF rank testing (gProfiler:Average Fscore[PS,NS,NUC,CB,PML combined]) Top 10 ctl subtraction")+
  xlab("Total Number of Ranks")+
  ylab("Average F2score F2-score ((1+2^2)*Prec*Recall)/(2^2*Prec+Recall)")+
  #xlim(c(0,1))+
  ylim(c(0,0.4))+
  theme_bw()
dev.off()

Fscores_avg_combined <- cbind(F0.5scores,F1scores,F2scores)
Fscores_avg_combined <- Fscores_avg_combined[,c(1,2,4,6)]
write.table(Fscores_avg_combined, "4_gProfile_Results_NMF_rank_test_avgFscores.txt", sep="\t", row.names=FALSE)
