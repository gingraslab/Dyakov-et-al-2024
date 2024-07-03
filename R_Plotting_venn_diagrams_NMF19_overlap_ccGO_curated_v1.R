# Plotting venn diagrams to show overlap between NMF/PPCOR results and GOcc results

## Part 1: Input all files and calculate intersection between different groups
## Part 2: Plot venn diagram to show overlap

library(eulerr)

#------
## Part 1: Input all files and calculate intersection between different groups
NMF_matrix <- read.delim("0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", header=TRUE)
GOcc_curated <- read.delim("GOcc_uniq_proteins_AmiGO_20240304.txt", sep="\t", header=TRUE)
Boris_curated <- read.delim("Boris_curated_lists_wIF_2.csv", sep=",", header=TRUE)

# NMF matrix filtered
PS_NMF_matrix <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X0" | 
                            NMF_matrix$Secondary_cluster_annotation == "X0" |
                            NMF_matrix$Tertiary_cluster_annotation == "X0"),] 

NS_NMF_matrix <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X12" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X12" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X12"),] 

NUC_NMF_matrix <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X4" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X4" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X4"),] 

CB_NMF_matrix <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X10" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X10" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X10"),] 

PML_NMF_matrix <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X14" | 
                                    NMF_matrix$Secondary_cluster_annotation == "X14" |
                                    NMF_matrix$Tertiary_cluster_annotation == "X14"),] 

# Plot venn diagrams (PS)
PS_GOcc <- unique(na.omit(GOcc_curated$PS_GOcc_uniq))
PS_NMF <- PS_NMF_matrix$PreyGene
PS_curated <- na.omit(Boris_curated$PS_Boris_Curated[nzchar(Boris_curated$PS_Boris_Curated)])

PS_input <- list("PS_GOcc" = PS_GOcc, 
             "PS_NMF" = PS_NMF, 
             "PS_curated" = PS_curated)

pdf("Rplot_venn_diagram_PS_NMF19_overlap_GOcc_curated.pdf", width = 5, height = 5)
plot(euler(PS_input), quantities = TRUE)
dev.off()

# Plot venn diagrams (NS)
NS_GOcc <- unique(na.omit(GOcc_curated$NS_GOcc_uniq))
NS_NMF <- NS_NMF_matrix$PreyGene
NS_curated <- na.omit(Boris_curated$NS_Boris_Curated[nzchar(Boris_curated$NS_Boris_Curated)])

NS_input <- list("NS_GOcc" = NS_GOcc, 
                 "NS_NMF" = NS_NMF, 
                 "NS_curated" = NS_curated)

pdf("Rplot_venn_diagram_NS_NMF19_overlap_GOcc_curated.pdf", width = 5, height = 5)
plot(euler(NS_input), quantities = TRUE)
dev.off()

# Plot venn diagrams (NUC)
NUC_GOcc <- unique(na.omit(GOcc_curated$NUC_GOcc_uniq))
NUC_NMF <- NUC_NMF_matrix$PreyGene
NUC_curated <- na.omit(Boris_curated$Nucleolus_Boris_Curated[nzchar(Boris_curated$Nucleolus_Boris_Curated)])

NUC_input <- list("NUC_GOcc" = NUC_GOcc, 
                 "NUC_NMF" = NUC_NMF, 
                 "NUC_curated" = NUC_curated)

pdf("Rplot_venn_diagram_NUC_NMF19_overlap_GOcc_curated.pdf", width = 5, height = 5)
plot(euler(NUC_input), quantities = TRUE)
dev.off()

# Plot venn diagrams (CB)
CB_GOcc <- unique(na.omit(GOcc_curated$CB_GOcc_uniq))
CB_NMF <- CB_NMF_matrix$PreyGene
CB_curated <- na.omit(Boris_curated$CB_Boris_Curated[nzchar(Boris_curated$CB_Boris_Curated)])

CB_input <- list("CB_GOcc" = CB_GOcc, 
                 "CB_NMF" = CB_NMF, 
                 "CB_curated" = CB_curated)

pdf("Rplot_venn_diagram_CB_NMF19_overlap_GOcc_curated.pdf", width = 5, height = 5)
plot(euler(CB_input), quantities = TRUE)
dev.off()

# Plot venn diagrams (PML)
PML_GOcc <- unique(na.omit(GOcc_curated$PML_GOcc_uniq))
PML_NMF <- PML_NMF_matrix$PreyGene
PML_curated <- na.omit(Boris_curated$PML_Boris_Curated[nzchar(Boris_curated$PML_Boris_Curated)])

PML_input <- list("PML_GOcc" = PML_GOcc, 
                 "PML_NMF" = PML_NMF, 
                 "PML_curated" = PML_curated)

pdf("Rplot_venn_diagram_PML_NMF19_overlap_GOcc_curated.pdf", width = 5, height = 5)
plot(euler(PML_input), quantities = TRUE)
dev.off()



