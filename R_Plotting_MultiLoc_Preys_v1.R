# Plot heatmaps showing co-localizing preys between ranks

## Part1: Read in annotated basis matrix and count number of overlaps and create matrix
## Part2: Plot heatmaps showing the number of overlap between each rank

#--------
## Part1: Read in annotated basis matrix and count number of overlaps and create matrix
NMF_matrix <- read.delim("0.5_nmf_prey_matrix_input_wCol_ratios_19ranks_0.05_primary(percent)_0.70_multi.csv", sep=",", header=TRUE)

# Matrix
NMF_matrix_X0 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X0"),]
NMF_matrix_X1 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X1"),]
NMF_matrix_X2 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X2"),]
NMF_matrix_X3 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X3"),]
NMF_matrix_X4 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X4"),]
NMF_matrix_X5 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X5"),]
NMF_matrix_X6 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X6"),]
NMF_matrix_X7 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X7"),]
NMF_matrix_X8 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X8"),]
NMF_matrix_X9 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X9"),]
NMF_matrix_X10 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X10"),]
NMF_matrix_X11 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X11"),]
NMF_matrix_X12 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X12"),]
NMF_matrix_X13 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X13"),]
NMF_matrix_X14 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X14"),]
NMF_matrix_X15 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X15"),]
NMF_matrix_X16 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X16"),]
NMF_matrix_X17 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X17"),]
NMF_matrix_X18 <- NMF_matrix[which(NMF_matrix$Primary_cluster_annotation == "X18"),]

#--------
# Create vectors with 
col_names <- c("X0_Paraspeckle","X1_Transcription_reg","X2_TF-Paraspeckle","X3_ER-NuclearMembrane","X4_Nucleolus","X5_Crytoplasmic_RNP_Granules","X6_Splicing","X7_Cytoskeleton-centrosome","X8_DNA-damage","X9_Transcription_reg","X10_Cajal_Body","X11_RNA_Export_Nucleoplasm","X12_Nuclear_Speckle","X13_Nuclear_Pore_Envelope","X14_PML_Body","X15_ER_Nucleolus","X16_Translation","X17_Polycomb","X18_Other")
row_names <- c("X0_Paraspeckle","X1_Transcription_reg","X2_TF-Paraspeckle","X3_ER-NuclearMembrane","X4_Nucleolus","X5_Crytoplasmic_RNP_Granules","X6_Splicing","X7_Cytoskeleton-centrosome","X8_DNA-damage","X9_Transcription_reg","X10_Cajal_Body","X11_RNA_Export_Nucleoplasm","X12_Nuclear_Speckle","X13_Nuclear_Pore_Envelope","X14_PML_Body","X15_ER_Nucleolus","X16_Translation","X17_Polycomb","X18_Other")


X0 <- c(nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X0" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X0"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X1" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X1"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X2" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X2"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X3" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X3"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X4" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X4"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X5" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X5"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X6" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X6"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X7" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X7"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X8" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X8"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X9" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X9"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X10" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X10"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X11" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X11"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X12" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X12"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X13" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X13"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X14" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X14"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X15" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X15"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X16" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X16"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X17" | 
                                       NMF_matrix_X0$Tertiary_cluster_annotation == "X17"),]),
            nrow(NMF_matrix_X0[which(NMF_matrix_X0$Secondary_cluster_annotation == "X18" | 
                           NMF_matrix_X0$Tertiary_cluster_annotation == "X18"),])
            )

X1 <- c(nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X0" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X0"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X1" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X1"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X2" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X2"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X3" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X3"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X4" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X4"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X5" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X5"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X6" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X6"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X7" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X7"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X8" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X8"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X9" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X9"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X10" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X10"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X11" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X11"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X12" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X12"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X13" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X13"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X14" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X14"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X15" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X15"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X16" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X16"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X17" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X17"),]),
            nrow(NMF_matrix_X1[which(NMF_matrix_X1$Secondary_cluster_annotation == "X18" | 
                                       NMF_matrix_X1$Tertiary_cluster_annotation == "X18"),])
)

X2 <- c(nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X2[which(NMF_matrix_X2$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X2$Tertiary_cluster_annotation == "X18"),])
)

X3 <- c(nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X3[which(NMF_matrix_X3$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X3$Tertiary_cluster_annotation == "X18"),])
)

X4 <- c(nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X4[which(NMF_matrix_X4$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X4$Tertiary_cluster_annotation == "X18"),])
)

X5 <- c(nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X5[which(NMF_matrix_X5$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X5$Tertiary_cluster_annotation == "X18"),])
)

X6 <- c(nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X6[which(NMF_matrix_X6$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X6$Tertiary_cluster_annotation == "X18"),])
)

X7 <- c(nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X7[which(NMF_matrix_X7$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X7$Tertiary_cluster_annotation == "X18"),])
)

X8 <- c(nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X8[which(NMF_matrix_X8$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X8$Tertiary_cluster_annotation == "X18"),])
)

X9 <- c(nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X9[which(NMF_matrix_X9$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X9$Tertiary_cluster_annotation == "X18"),])
)

X10 <- c(nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X10$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X10[which(NMF_matrix_X10$Secondary_cluster_annotation == "X18" | 
                                    NMF_matrix_X10$Tertiary_cluster_annotation == "X18"),])
)

X11 <- c(nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X11$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X11[which(NMF_matrix_X11$Secondary_cluster_annotation == "X18" | 
                                    NMF_matrix_X11$Tertiary_cluster_annotation == "X18"),])
)

X12 <- c(nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X12$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X12[which(NMF_matrix_X12$Secondary_cluster_annotation == "X18" | 
                                    NMF_matrix_X12$Tertiary_cluster_annotation == "X18"),])
)

X13 <- c(nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X13[which(NMF_matrix_X13$Secondary_cluster_annotation == "X18" | 
                                   NMF_matrix_X13$Tertiary_cluster_annotation == "X18"),])
)

X14 <- c(nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X0" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X0"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X1" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X1"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X2" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X2"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X3" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X3"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X4" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X4"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X5" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X5"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X6" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X6"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X7" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X7"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X8" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X8"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X9" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X9"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X10" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X10"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X11" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X11"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X12" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X12"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X13" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X13"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X14" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X14"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X15" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X15"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X16" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X16"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X17" | 
                                   NMF_matrix_X14$Tertiary_cluster_annotation == "X17"),]),
        nrow(NMF_matrix_X14[which(NMF_matrix_X14$Secondary_cluster_annotation == "X18" | 
                                    NMF_matrix_X14$Tertiary_cluster_annotation == "X18"),])
)

X15 <- c(nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X0" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X0"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X1" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X1"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X2" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X2"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X3" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X3"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X4" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X4"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X5" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X5"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X6" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X6"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X7" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X7"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X8" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X8"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X9" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X9"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X10" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X10"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X11" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X11"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X12" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X12"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X13" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X13"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X14" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X14"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X15" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X15"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X16" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X16"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X17" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X17"),]),
         nrow(NMF_matrix_X15[which(NMF_matrix_X15$Secondary_cluster_annotation == "X18" | 
                                     NMF_matrix_X15$Tertiary_cluster_annotation == "X18"),])
)

X16 <- c(nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X0" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X0"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X1" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X1"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X2" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X2"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X3" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X3"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X4" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X4"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X5" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X5"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X6" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X6"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X7" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X7"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X8" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X8"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X9" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X9"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X10" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X10"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X11" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X11"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X12" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X12"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X13" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X13"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X14" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X14"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X15" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X15"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X16" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X16"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X17" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X17"),]),
         nrow(NMF_matrix_X16[which(NMF_matrix_X16$Secondary_cluster_annotation == "X18" | 
                                     NMF_matrix_X16$Tertiary_cluster_annotation == "X18"),])
)

X17 <- c(nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X0" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X0"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X1" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X1"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X2" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X2"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X3" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X3"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X4" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X4"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X5" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X5"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X6" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X6"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X7" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X7"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X8" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X8"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X9" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X9"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X10" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X10"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X11" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X11"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X12" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X12"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X13" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X13"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X14" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X14"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X15" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X15"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X16" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X16"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X17" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X17"),]),
         nrow(NMF_matrix_X17[which(NMF_matrix_X17$Secondary_cluster_annotation == "X18" | 
                                     NMF_matrix_X17$Tertiary_cluster_annotation == "X18"),])
)

X18 <- c(nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X0" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X0"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X1" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X1"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X2" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X2"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X3" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X3"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X4" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X4"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X5" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X5"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X6" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X6"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X7" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X7"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X8" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X8"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X9" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X9"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X10" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X10"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X11" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X11"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X12" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X12"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X13" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X13"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X14" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X14"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X15" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X15"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X16" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X16"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X17" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X17"),]),
         nrow(NMF_matrix_X18[which(NMF_matrix_X18$Secondary_cluster_annotation == "X18" | 
                                     NMF_matrix_X18$Tertiary_cluster_annotation == "X18"),])
)

multiloc_overlap_matrix <- rbind(X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18)
colnames(multiloc_overlap_matrix) <- col_names
rownames(multiloc_overlap_matrix) <- row_names

write.table(multiloc_overlap_matrix, "1_NMF19_counts_multi_loc_preys_col_ratios_primary_0.05_multi_0.70.txt", sep="\t", row.names=TRUE)

#-----------
#Plot heatmap (NOTE: For showing row names need to alter the cEX number size!!! (or may loose rows))
library("gplots")
library("RColorBrewer")

cols <- colorRampPalette(c("white", "blue"))(1000)
#hmcols <- rev(redblue(100))

pdf("Rplot_Heatmap_NMF19_multiloc_number_col_ratio_0.05_primary_0.70_multi.pdf", width = 9, height = 8)
heatmap_plot <- heatmap.2(multiloc_overlap_matrix, cellnote=ifelse(multiloc_overlap_matrix==0, NA, multiloc_overlap_matrix), notecol = "black", dendrogram='none', Colv = FALSE, Rowv = FALSE, scale=NULL, trace='none', col=cols, density.info="none", main = "Multi-loc preys between ranks (multi-score >= 0.50)", cexRow = 1, cexCol =  1)
dev.off()

