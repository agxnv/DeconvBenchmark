#Figure_1b: 
rm(list=ls())
library(immunedeconv)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(MIXTURE)


##Extract estimated proportions for each signature

##Kidney

K_sig <- read.csv("Figure_1/Data/MS/K.csv", row.names = 1, check.names=FALSE)

K_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/Kidney.CIBERSORT.txt") %>% arrange(Mixture)

K_EPIC <- EPIC::EPIC(bulk = K_sig,
                     reference = list(refProfiles = K_sig, sigGenes = rownames(K_sig)),
                     withOtherCells = FALSE)
K_EPIC <- as.data.frame(K_EPIC$cellFractions)
K_EPIC$X1 <- rownames(K_EPIC)
K_EPIC <- K_EPIC %>% arrange(X1)

K_MIXTURE <- MIXTURE(expressionMatrix = K_sig,
                signatureMatrix = K_sig, 
                iter = 0L,
                functionMixture = nu.svm.robust.RFE,
                useCores = 1,
                verbose = TRUE, 
                nullDist = "PopulationBased",  
                fileSave = "") 
K_MIXTURE <- as.data.frame(K_MIXTURE$Subjects$MIXprop)
K_MIXTURE$X1 <- rownames(K_MIXTURE)
K_MIXTURE <- K_MIXTURE %>% arrange(X1)


rownames(K_sig) <- toupper(rownames(K_sig))
write.table(data.frame("X1"=rownames(K_sig),K_sig), file = "Figure_1/Data/MS/K_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
K_quanTIseq <-deconvolute_quantiseq.default(mix.mat = K_sig, 
                                           arrays = FALSE, 
                                           signame = "Figure_1/Data/MS/K_sig", 
                                           tumor = FALSE, 
                                           mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                           rmgenes = "unassigned")[, -c(10)]
K_quanTIseq <- K_quanTIseq %>% arrange(Sample)


K_CIBERSORT_proportions <- data.matrix(K_CIBERSORT[,2:9])
K_CIBERSORT_proportions <-  K_CIBERSORT_proportions[,order(colnames(K_CIBERSORT_proportions))]
K_EPIC_proportions <- data.matrix(K_EPIC[,1:8])
K_EPIC_proportions <- K_EPIC_proportions[,order(colnames(K_EPIC_proportions))]
K_MIXTURE_proportions <- data.matrix(K_MIXTURE[,1:8])
K_MIXTURE_proportions <- K_MIXTURE_proportions[,order(colnames(K_MIXTURE_proportions))]
K_quanTIseq_proportions <- data.matrix(K_quanTIseq[,2:9])
K_quanTIseq_proportions <- K_quanTIseq_proportions[,order(colnames(K_quanTIseq_proportions))]
K_quanTIseq_proportions[K_quanTIseq_proportions == 0] <- NA
K_CIBERSORT_proportions[K_CIBERSORT_proportions == 0] <- NA
K_EPIC_proportions[K_EPIC_proportions == 0] <- NA
K_MIXTURE_proportions[K_MIXTURE_proportions == 0] <- NA

##Liver

L_sig <- read.csv("Figure_1/Data/MS/L.csv", row.names = 1, check.names=FALSE)

L_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/Liver.CIBERSORT.txt") %>% arrange(Mixture)

L_EPIC <- EPIC::EPIC(bulk = L_sig,
                     reference = list(refProfiles = L_sig, sigGenes = rownames(L_sig)),
                     withOtherCells = FALSE)
L_EPIC <- as.data.frame(L_EPIC$cellFractions)
L_EPIC$X1 <- rownames(L_EPIC)
L_EPIC <- L_EPIC %>% arrange(X1)

L_MIXTURE <- MIXTURE(expressionMatrix = L_sig,
                     signatureMatrix = L_sig, 
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE, 
                     nullDist = "PopulationBased",  
                     fileSave = "") 
L_MIXTURE <- as.data.frame(L_MIXTURE$Subjects$MIXprop)
L_MIXTURE$X1 <- rownames(L_MIXTURE)
L_MIXTURE <- L_MIXTURE %>% arrange(X1)


rownames(L_sig) <- toupper(rownames(L_sig))
write.table(data.frame("X1"=rownames(L_sig),L_sig), file = "Figure_1/Data/MS/L_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
L_quanTIseq <-deconvolute_quantiseq.default(mix.mat = L_sig, 
                                            arrays = FALSE, 
                                            signame = "Figure_1/Data/MS/L_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
L_quanTIseq <- L_quanTIseq %>% arrange(Sample)


L_CIBERSORT_proportions <- data.matrix(L_CIBERSORT[,2:9])
L_CIBERSORT_proportions <-  L_CIBERSORT_proportions[,order(colnames(L_CIBERSORT_proportions))]
L_EPIC_proportions <- data.matrix(L_EPIC[,1:8])
L_EPIC_proportions <- L_EPIC_proportions[,order(colnames(L_EPIC_proportions))]
L_MIXTURE_proportions <- data.matrix(L_MIXTURE[,1:8])
L_MIXTURE_proportions <- L_MIXTURE_proportions[,order(colnames(L_MIXTURE_proportions))]
L_quanTIseq_proportions <- data.matrix(L_quanTIseq[,2:9])
L_quanTIseq_proportions <- L_quanTIseq_proportions[,order(colnames(L_quanTIseq_proportions))]
L_quanTIseq_proportions[L_quanTIseq_proportions == 0] <- NA
L_CIBERSORT_proportions[L_CIBERSORT_proportions == 0] <- NA
L_EPIC_proportions[L_EPIC_proportions == 0] <- NA
L_MIXTURE_proportions[L_MIXTURE_proportions == 0] <- NA

#Mammary gland

MG_sig <- read.csv("Figure_1/Data/MS/MG.csv", row.names = 1, check.names=FALSE)

MG_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/Mgland.CIBERSORT.txt") %>% arrange(Mixture)

MG_EPIC <- EPIC::EPIC(bulk = MG_sig,
                      reference = list(refProfiles = MG_sig, sigGenes = rownames(MG_sig)),
                      withOtherCells = FALSE)
MG_EPIC <- as.data.frame(MG_EPIC$cellFractions)
MG_EPIC$X1 <- rownames(MG_EPIC)
MG_EPIC <- MG_EPIC %>% arrange(X1)

MG_MIXTURE <- MIXTURE(expressionMatrix = MG_sig,
                      signatureMatrix = MG_sig, 
                      iter = 0L,
                      functionMixture = nu.svm.robust.RFE,
                      useCores = 1,
                      verbose = TRUE, 
                      nullDist = "PopulationBased",  
                      fileSave = "") 
MG_MIXTURE <- as.data.frame(MG_MIXTURE$Subjects$MIXprop)
MG_MIXTURE$X1 <- rownames(MG_MIXTURE)
MG_MIXTURE <- MG_MIXTURE %>% arrange(X1)


rownames(MG_sig) <- toupper(rownames(MG_sig))
write.table(data.frame("X1"=rownames(MG_sig),MG_sig), file = "Figure_1/Data/MS/MG_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
MG_quanTIseq <-deconvolute_quantiseq.default(mix.mat = MG_sig, 
                                             arrays = FALSE, 
                                             signame = "Figure_1/Data/MS/MG_sig", 
                                             tumor = FALSE, 
                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                             rmgenes = "unassigned")[, -c(10)]
MG_quanTIseq <- MG_quanTIseq %>% arrange(Sample)


MG_CIBERSORT_proportions <- data.matrix(MG_CIBERSORT[,2:9])
MG_CIBERSORT_proportions <-  MG_CIBERSORT_proportions[,order(colnames(MG_CIBERSORT_proportions))]
MG_EPIC_proportions <- data.matrix(MG_EPIC[,1:8])
MG_EPIC_proportions <- MG_EPIC_proportions[,order(colnames(MG_EPIC_proportions))]
MG_MIXTURE_proportions <- data.matrix(MG_MIXTURE[,1:8])
MG_MIXTURE_proportions <- MG_MIXTURE_proportions[,order(colnames(MG_MIXTURE_proportions))]
MG_quanTIseq_proportions <- data.matrix(MG_quanTIseq[,2:9])
MG_quanTIseq_proportions <- MG_quanTIseq_proportions[,order(colnames(MG_quanTIseq_proportions))]
MG_quanTIseq_proportions[MG_quanTIseq_proportions == 0] <- NA
MG_CIBERSORT_proportions[MG_CIBERSORT_proportions == 0] <- NA
MG_EPIC_proportions[MG_EPIC_proportions == 0] <- NA
MG_MIXTURE_proportions[MG_MIXTURE_proportions == 0] <- NA

##Muscle

M_sig <- read.csv("Figure_1/Data/MS/M.csv", row.names = 1, check.names=FALSE)

M_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/Muscle.CIBERSORT.txt") %>% arrange(Mixture)

M_EPIC <- EPIC::EPIC(bulk = M_sig,
                     reference = list(refProfiles = M_sig, sigGenes = rownames(M_sig)),
                     withOtherCells = FALSE)
M_EPIC <- as.data.frame(M_EPIC$cellFractions)
M_EPIC$X1 <- rownames(M_EPIC)
M_EPIC <- M_EPIC %>% arrange(X1)

M_MIXTURE <- MIXTURE(expressionMatrix = M_sig,
                     signatureMatrix = M_sig, 
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE, 
                     nullDist = "PopulationBased",  
                     fileSave = "") 
M_MIXTURE <- as.data.frame(M_MIXTURE$Subjects$MIXprop)
M_MIXTURE$X1 <- rownames(M_MIXTURE)
M_MIXTURE <- M_MIXTURE %>% arrange(X1)


rownames(M_sig) <- toupper(rownames(M_sig))
write.table(data.frame("X1"=rownames(M_sig),M_sig), file = "Figure_1/Data/MS/M_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
M_quanTIseq <-deconvolute_quantiseq.default(mix.mat = M_sig, 
                                            arrays = FALSE, 
                                            signame = "Figure_1/Data/MS/M_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
M_quanTIseq <- M_quanTIseq %>% arrange(Sample)


M_CIBERSORT_proportions <- data.matrix(M_CIBERSORT[,2:9])
M_CIBERSORT_proportions <-  M_CIBERSORT_proportions[,order(colnames(M_CIBERSORT_proportions))]
M_EPIC_proportions <- data.matrix(M_EPIC[,1:8])
M_EPIC_proportions <- M_EPIC_proportions[,order(colnames(M_EPIC_proportions))]
M_MIXTURE_proportions <- data.matrix(M_MIXTURE[,1:8])
M_MIXTURE_proportions <- M_MIXTURE_proportions[,order(colnames(M_MIXTURE_proportions))]
M_quanTIseq_proportions <- data.matrix(M_quanTIseq[,2:9])
M_quanTIseq_proportions <- M_quanTIseq_proportions[,order(colnames(M_quanTIseq_proportions))]
M_quanTIseq_proportions[M_quanTIseq_proportions == 0] <- NA
M_CIBERSORT_proportions[M_CIBERSORT_proportions == 0] <- NA
M_EPIC_proportions[M_EPIC_proportions == 0] <- NA
M_MIXTURE_proportions[M_MIXTURE_proportions == 0] <- NA

#Small intestine

SI_sig <- read.csv("Figure_1/Data/MS/SI.csv", row.names = 1, check.names=FALSE)

SI_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/SInt.CIBERSORT.txt") %>% arrange(Mixture)

SI_EPIC <- EPIC::EPIC(bulk = SI_sig,
                      reference = list(refProfiles = SI_sig, sigGenes = rownames(SI_sig)),
                      withOtherCells = FALSE)
SI_EPIC <- as.data.frame(SI_EPIC$cellFractions)
SI_EPIC$X1 <- rownames(SI_EPIC)
SI_EPIC <- SI_EPIC %>% arrange(X1)

SI_MIXTURE <- MIXTURE(expressionMatrix = SI_sig,
                      signatureMatrix = SI_sig, 
                      iter = 0L,
                      functionMixture = nu.svm.robust.RFE,
                      useCores = 1,
                      verbose = TRUE, 
                      nullDist = "PopulationBased",  
                      fileSave = "") 
SI_MIXTURE <- as.data.frame(SI_MIXTURE$Subjects$MIXprop)
SI_MIXTURE$X1 <- rownames(SI_MIXTURE)
SI_MIXTURE <- SI_MIXTURE %>% arrange(X1)


rownames(SI_sig) <- toupper(rownames(SI_sig))
write.table(data.frame("X1"=rownames(SI_sig),SI_sig), file = "Figure_1/Data/MS/SI_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
SI_quanTIseq <-deconvolute_quantiseq.default(mix.mat = SI_sig, 
                                             arrays = FALSE, 
                                             signame = "Figure_1/Data/MS/SI_sig", 
                                             tumor = FALSE, 
                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                             rmgenes = "unassigned")[, -c(10)]
SI_quanTIseq <- SI_quanTIseq %>% arrange(Sample)


SI_CIBERSORT_proportions <- data.matrix(SI_CIBERSORT[,2:9])
SI_CIBERSORT_proportions <-  SI_CIBERSORT_proportions[,order(colnames(SI_CIBERSORT_proportions))]
SI_EPIC_proportions <- data.matrix(SI_EPIC[,1:8])
SI_EPIC_proportions <- SI_EPIC_proportions[,order(colnames(SI_EPIC_proportions))]
SI_MIXTURE_proportions <- data.matrix(SI_MIXTURE[,1:8])
SI_MIXTURE_proportions <- SI_MIXTURE_proportions[,order(colnames(SI_MIXTURE_proportions))]
SI_quanTIseq_proportions <- data.matrix(SI_quanTIseq[,2:9])
SI_quanTIseq_proportions <- SI_quanTIseq_proportions[,order(colnames(SI_quanTIseq_proportions))]
SI_quanTIseq_proportions[SI_quanTIseq_proportions == 0] <- NA
SI_CIBERSORT_proportions[SI_CIBERSORT_proportions == 0] <- NA
SI_EPIC_proportions[SI_EPIC_proportions == 0] <- NA
SI_MIXTURE_proportions[SI_MIXTURE_proportions == 0] <- NA

#Spleen

S_sig <- read.csv("Figure_1/Data/MS/S.csv", row.names = 1, check.names=FALSE)

S_CIBERSORT <- read_table2("Figure_1/Data/CIBERSORT_estimations/Spleen.CIBERSORT.txt") %>% arrange(Mixture)

S_EPIC <- EPIC::EPIC(bulk = S_sig,
                     reference = list(refProfiles = S_sig, sigGenes = rownames(S_sig)),
                     withOtherCells = FALSE)
S_EPIC <- as.data.frame(S_EPIC$cellFractions)
S_EPIC$X1 <- rownames(S_EPIC)
S_EPIC <- S_EPIC %>% arrange(X1)

S_MIXTURE <- MIXTURE(expressionMatrix = S_sig,
                     signatureMatrix = S_sig, 
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE, 
                     nullDist = "PopulationBased",  
                     fileSave = "") 
S_MIXTURE <- as.data.frame(S_MIXTURE$Subjects$MIXprop)
S_MIXTURE$X1 <- rownames(S_MIXTURE)
S_MIXTURE <- S_MIXTURE %>% arrange(X1)


rownames(S_sig) <- toupper(rownames(S_sig))
write.table(data.frame("X1"=rownames(S_sig),S_sig), file = "Figure_1/Data/S_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
S_quanTIseq <-deconvolute_quantiseq.default(mix.mat = S_sig, 
                                            arrays = FALSE, 
                                            signame = "Figure_1/Data/S_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
S_quanTIseq <- S_quanTIseq %>% arrange(Sample)


S_CIBERSORT_proportions <- data.matrix(S_CIBERSORT[,2:9])
S_CIBERSORT_proportions <-  S_CIBERSORT_proportions[,order(colnames(S_CIBERSORT_proportions))]
S_EPIC_proportions <- data.matrix(S_EPIC[,1:8])
S_EPIC_proportions <- S_EPIC_proportions[,order(colnames(S_EPIC_proportions))]
S_MIXTURE_proportions <- data.matrix(S_MIXTURE[,1:8])
S_MIXTURE_proportions <- S_MIXTURE_proportions[,order(colnames(S_MIXTURE_proportions))]
S_quanTIseq_proportions <- data.matrix(S_quanTIseq[,2:9])
S_quanTIseq_proportions <- S_quanTIseq_proportions[,order(colnames(S_quanTIseq_proportions))]
S_quanTIseq_proportions[S_quanTIseq_proportions == 0] <- NA
S_CIBERSORT_proportions[S_CIBERSORT_proportions == 0] <- NA
S_EPIC_proportions[S_EPIC_proportions == 0] <- NA
S_MIXTURE_proportions[S_MIXTURE_proportions == 0] <- NA

#Renaming cell types for graphic 
cell.types.names <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")

colnames(K_quanTIseq_proportions) <- rownames(K_quanTIseq_proportions) <- cell.types.names
colnames(K_CIBERSORT_proportions) <- rownames(K_CIBERSORT_proportions) <- cell.types.names
colnames(K_EPIC_proportions) <- rownames(K_EPIC_proportions) <- cell.types.names
colnames(K_MIXTURE_proportions) <- rownames(K_MIXTURE_proportions) <- cell.types.names

colnames(L_quanTIseq_proportions) <- rownames(L_quanTIseq_proportions) <- cell.types.names
colnames(L_CIBERSORT_proportions) <- rownames(L_CIBERSORT_proportions) <- cell.types.names
colnames(L_EPIC_proportions) <- rownames(L_EPIC_proportions) <- cell.types.names
colnames(L_MIXTURE_proportions) <- rownames(L_MIXTURE_proportions) <- cell.types.names

colnames(MG_quanTIseq_proportions) <- rownames(MG_quanTIseq_proportions) <- cell.types.names
colnames(MG_CIBERSORT_proportions) <- rownames(MG_CIBERSORT_proportions) <- cell.types.names
colnames(MG_EPIC_proportions) <- rownames(MG_EPIC_proportions) <- cell.types.names
colnames(MG_MIXTURE_proportions) <- rownames(MG_MIXTURE_proportions) <- cell.types.names

colnames(M_quanTIseq_proportions) <- rownames(M_quanTIseq_proportions) <- cell.types.names
colnames(M_CIBERSORT_proportions) <- rownames(M_CIBERSORT_proportions) <- cell.types.names
colnames(M_EPIC_proportions) <- rownames(M_EPIC_proportions) <- cell.types.names
colnames(M_MIXTURE_proportions) <- rownames(M_MIXTURE_proportions) <- cell.types.names

colnames(SI_quanTIseq_proportions) <- rownames(SI_quanTIseq_proportions) <- cell.types.names
colnames(SI_CIBERSORT_proportions) <- rownames(SI_CIBERSORT_proportions) <- cell.types.names
colnames(SI_EPIC_proportions) <- rownames(SI_EPIC_proportions) <- cell.types.names
colnames(SI_MIXTURE_proportions) <- rownames(SI_MIXTURE_proportions) <- cell.types.names

colnames(S_quanTIseq_proportions) <- rownames(S_quanTIseq_proportions) <- cell.types.names
colnames(S_CIBERSORT_proportions) <- rownames(S_CIBERSORT_proportions) <- cell.types.names
colnames(S_EPIC_proportions) <- rownames(S_EPIC_proportions) <- cell.types.names
colnames(S_MIXTURE_proportions) <- rownames(S_MIXTURE_proportions) <- cell.types.names

#Heatmaps
Kidney_heatmap <- 
  Heatmap(K_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Kidney",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(K_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(K_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(K_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)


Liver_heatmap <- 
  Heatmap(L_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Liver",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(L_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(L_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(L_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

MGland_heatmap <- 
  Heatmap(MG_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Mammary gland",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(MG_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MG_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(MG_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

Muscle_heatmap <- 
  Heatmap(M_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Muscle",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(M_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(M_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(M_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

SInt_heatmap <- 
  Heatmap(SI_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Small intestine",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(SI_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(SI_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(SI_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

Spleen_heatmap <- 
  Heatmap(S_CIBERSORT_proportions, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORT",name = "Estimated \nCoefficients",row_title = "Spleen",row_title_gp = gpar(fontsize=20), column_title_gp = gpar(fontsize=20),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = TRUE) +
  Heatmap(S_EPIC_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(S_MIXTURE_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(S_quanTIseq_proportions, cluster_rows = FALSE, column_title_gp = gpar(fontsize=20), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

#Save plot in PDF 12x6
pdf("Spleen_ST.pdf",width=12,height=6)
print(Spleen_heatmap)
dev.off()


#Me falta: Leyenda NAs


