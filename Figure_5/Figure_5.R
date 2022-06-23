library(MIXTURE)
library(readr)
library(immunedeconv)
library(tidyverse)
library(ggpubr)
library(reshape2)

load("Figure_5/Data/simulated_immune.RData")

k_sig <- read.csv("Figure_5/Data/MS/K.csv", row.names = 1, check.names=FALSE)

k_CIBERSORT_inm_sum <- read_csv("Figure_5/Data/CIBERSORT_estimations/k_CIBERSORT_inm_sum_prop.csv")
k_CIBERSORT_inm_sum <- as.data.frame(k_CIBERSORT_inm_sum)
rownames(k_CIBERSORT_inm_sum) <- k_CIBERSORT_inm_sum[,1]
k_CIBERSORT_inm_sum <- k_CIBERSORT_inm_sum[-c((nrow(k_CIBERSORT_inm_sum)-2):((nrow(k_CIBERSORT_inm_sum)))), -c(1,10,11,12,13)]
k_CIBERSORT_inm_sum$Method <- "CIBERSORT"
k_CIBERSORT_inm_sum$Signature <- "Kidney"
k_CIBERSORT_inm_sum$Celltype <- rownames(k_CIBERSORT_inm_sum) 

k_EPIC_inm_sum <- EPIC::EPIC(bulk = k_counts_inm_sum,
                     reference = list(refProfiles = k_sig, sigGenes = rownames(k_sig)),
                     withOtherCells = FALSE)
k_EPIC_inm_sum <- as.data.frame(k_EPIC_inm_sum$cellFractions)
k_EPIC_inm_sum$Method <- "EPIC"
k_EPIC_inm_sum$Signature <- "Kidney"
k_EPIC_inm_sum$Celltype <- rownames(k_EPIC_inm_sum) 

k_MIXTURE_inm_sum <- MIXTURE(expressionMatrix = k_counts_inm_sum,      
                             signatureMatrix = k_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
k_MIXTURE_inm_sum <- as.data.frame(k_MIXTURE_inm_sum$Subjects$MIXprop)
k_MIXTURE_inm_sum$Method <- "MIXTURE"
k_MIXTURE_inm_sum$Signature <- "Kidney"
k_MIXTURE_inm_sum$Celltype <- rownames(k_MIXTURE_inm_sum) 

rownames(k_sig) <- toupper(rownames(k_sig))
write.table(data.frame("X1"=rownames(k_sig),k_sig), file = "Figure_5/Data/MS/K_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(k_counts_inm_sum) <- toupper(rownames(k_counts_inm_sum))
k_quanTIseq_inm_sum <-deconvolute_quantiseq.default(mix.mat = k_counts_inm_sum, 
                                            arrays = FALSE, 
                                            signame = "Figure_5/Data/MS/K_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
k_quanTIseq_inm_sum[,1] <- NULL
k_quanTIseq_inm_sum$Method <- "quanTIseq"
k_quanTIseq_inm_sum$Signature <- "Kidney"
k_quanTIseq_inm_sum$Celltype <- rownames(k_quanTIseq_inm_sum)
names(k_quanTIseq_inm_sum) <- gsub("[.]"," ", names(k_quanTIseq_inm_sum))

k_inm_est <- rbind(k_CIBERSORT_inm_sum,k_EPIC_inm_sum,k_MIXTURE_inm_sum,k_quanTIseq_inm_sum)

l_sig <- read.csv("Figure_5/Data/MS/L.csv", row.names = 1, check.names=FALSE)

l_CIBERSORT_inm_sum_1 <- read_csv("Figure_5/Data/CIBERSORT_estimations/l_CIBERSORT_inm_sum_prop_1.csv")
l_CIBERSORT_inm_sum_1 <- as.data.frame(l_CIBERSORT_inm_sum_1)
rownames(l_CIBERSORT_inm_sum_1) <- l_CIBERSORT_inm_sum_1[,1]
l_CIBERSORT_inm_sum_1 <- l_CIBERSORT_inm_sum_1[-c((nrow(l_CIBERSORT_inm_sum_1)-4):((nrow(l_CIBERSORT_inm_sum_1)))), -c(1,10,11,12,13)]
l_CIBERSORT_inm_sum_1$Method <- "CIBERSORT"
l_CIBERSORT_inm_sum_1$Signature <- "Liver"
l_CIBERSORT_inm_sum_1$Celltype <- rownames(l_CIBERSORT_inm_sum_1) 

l_EPIC_inm_sum_1 <- EPIC::EPIC(bulk = l_counts_inm_sum_1,
                               reference = list(refProfiles = l_sig, sigGenes = rownames(l_sig)),
                               withOtherCells = FALSE)
l_EPIC_inm_sum_1 <- as.data.frame(l_EPIC_inm_sum_1$cellFractions)
l_EPIC_inm_sum_1$Method <- "EPIC"
l_EPIC_inm_sum_1$Signature <- "Liver"
l_EPIC_inm_sum_1$Celltype <- rownames(l_EPIC_inm_sum_1) 

l_MIXTURE_inm_sum_1 <- MIXTURE(expressionMatrix = l_counts_inm_sum_1,      
                               signatureMatrix = l_sig,       
                               iter = 0L,
                               functionMixture = nu.svm.robust.RFE,
                               useCores = 1,
                               verbose = TRUE,
                               nullDist = "PopulationBased") 
l_MIXTURE_inm_sum_1 <- as.data.frame(l_MIXTURE_inm_sum_1$Subjects$MIXprop)
l_MIXTURE_inm_sum_1$Method <- "MIXTURE"
l_MIXTURE_inm_sum_1$Signature <- "Liver"
l_MIXTURE_inm_sum_1$Celltype <- rownames(l_MIXTURE_inm_sum_1) 

rownames(l_sig) <- toupper(rownames(l_sig))
write.table(data.frame("X1"=rownames(l_sig),l_sig), file = "Figure_5/Data/MS/L_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(l_counts_inm_sum_1) <- toupper(rownames(l_counts_inm_sum_1))
l_quanTIseq_inm_sum_1 <-deconvolute_quantiseq.default(mix.mat = l_counts_inm_sum_1, 
                                                      arrays = FALSE, 
                                                      signame = "Figure_5/Data/MS/L_sig", 
                                                      tumor = FALSE, 
                                                      mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                      rmgenes = "unassigned")[, -c(10)]
l_quanTIseq_inm_sum_1[,1] <- NULL
l_quanTIseq_inm_sum_1$Method <- "quanTIseq"
l_quanTIseq_inm_sum_1$Signature <- "Liver"
l_quanTIseq_inm_sum_1$Celltype <- rownames(l_quanTIseq_inm_sum_1)
names(l_quanTIseq_inm_sum_1) <- gsub("[.]"," ", names(l_quanTIseq_inm_sum_1))

l_sig <- read.csv("Figure_5/Data/MS/L.csv", row.names = 1, check.names=FALSE)

l_CIBERSORT_inm_sum_2 <- read_csv("Figure_5/Data/CIBERSORT_estimations/l_CIBERSORT_inm_sum_prop_2.csv")
l_CIBERSORT_inm_sum_2 <- as.data.frame(l_CIBERSORT_inm_sum_2)
rownames(l_CIBERSORT_inm_sum_2) <- l_CIBERSORT_inm_sum_2[,1]
l_CIBERSORT_inm_sum_2 <- l_CIBERSORT_inm_sum_2[-c((nrow(l_CIBERSORT_inm_sum_2)-4):((nrow(l_CIBERSORT_inm_sum_2)))), -c(1,10,11,12,13)]
l_CIBERSORT_inm_sum_2$Method <- "CIBERSORT"
l_CIBERSORT_inm_sum_2$Signature <- "Liver"
l_CIBERSORT_inm_sum_2$Celltype <- rownames(l_CIBERSORT_inm_sum_2) 

l_EPIC_inm_sum_2 <- EPIC::EPIC(bulk = l_counts_inm_sum_2,
                               reference = list(refProfiles = l_sig, sigGenes = rownames(l_sig)),
                               withOtherCells = FALSE)
l_EPIC_inm_sum_2 <- as.data.frame(l_EPIC_inm_sum_2$cellFractions)
l_EPIC_inm_sum_2$Method <- "EPIC"
l_EPIC_inm_sum_2$Signature <- "Liver"
l_EPIC_inm_sum_2$Celltype <- rownames(l_EPIC_inm_sum_2) 

l_MIXTURE_inm_sum_2 <- MIXTURE(expressionMatrix = l_counts_inm_sum_2,      
                               signatureMatrix = l_sig,       
                               iter = 0L,
                               functionMixture = nu.svm.robust.RFE,
                               useCores = 1,
                               verbose = TRUE,
                               nullDist = "PopulationBased") 
l_MIXTURE_inm_sum_2 <- as.data.frame(l_MIXTURE_inm_sum_2$Subjects$MIXprop)
l_MIXTURE_inm_sum_2$Method <- "MIXTURE"
l_MIXTURE_inm_sum_2$Signature <- "Liver"
l_MIXTURE_inm_sum_2$Celltype <- rownames(l_MIXTURE_inm_sum_2) 

rownames(l_counts_inm_sum_2) <- toupper(rownames(l_counts_inm_sum_2))
l_quanTIseq_inm_sum_2 <-deconvolute_quantiseq.default(mix.mat = l_counts_inm_sum_2, 
                                                      arrays = FALSE, 
                                                      signame = "Figure_5/Data/MS/L_sig", 
                                                      tumor = FALSE, 
                                                      mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                      rmgenes = "unassigned")[, -c(10)]
l_quanTIseq_inm_sum_2[,1] <- NULL
l_quanTIseq_inm_sum_2$Method <- "quanTIseq"
l_quanTIseq_inm_sum_2$Signature <- "Liver"
l_quanTIseq_inm_sum_2$Celltype <- rownames(l_quanTIseq_inm_sum_2)
names(l_quanTIseq_inm_sum_2) <- gsub("[.]"," ", names(l_quanTIseq_inm_sum_2))

l_inm_est <- rbind(l_CIBERSORT_inm_sum_1,l_CIBERSORT_inm_sum_2,
                   l_EPIC_inm_sum_1,l_EPIC_inm_sum_2,
                   l_MIXTURE_inm_sum_1,l_MIXTURE_inm_sum_2,
                   l_quanTIseq_inm_sum_1,l_quanTIseq_inm_sum_2)

mg_sig <- read.csv("Figure_5/Data/MS/MG.csv", row.names = 1, check.names=FALSE)

mg_CIBERSORT_inm_sum <- read_csv("Figure_5/Data/CIBERSORT_estimations/mg_CIBERSORT_inm_sum_prop.csv")
mg_CIBERSORT_inm_sum <- as.data.frame(mg_CIBERSORT_inm_sum)
rownames(mg_CIBERSORT_inm_sum) <- mg_CIBERSORT_inm_sum[,1]
mg_CIBERSORT_inm_sum <- mg_CIBERSORT_inm_sum[-c((nrow(mg_CIBERSORT_inm_sum)-2):((nrow(mg_CIBERSORT_inm_sum)))), -c(1,10,11,12,13)]
mg_CIBERSORT_inm_sum$Method <- "CIBERSORT"
mg_CIBERSORT_inm_sum$Signature <- "Mammary gland"
mg_CIBERSORT_inm_sum$Celltype <- rownames(mg_CIBERSORT_inm_sum) 

mg_EPIC_inm_sum <- EPIC::EPIC(bulk = mg_counts_inm_sum,
                             reference = list(refProfiles = mg_sig, sigGenes = rownames(mg_sig)),
                             withOtherCells = FALSE)
mg_EPIC_inm_sum <- as.data.frame(mg_EPIC_inm_sum$cellFractions)
mg_EPIC_inm_sum$Method <- "EPIC"
mg_EPIC_inm_sum$Signature <- "Mammary gland"
mg_EPIC_inm_sum$Celltype <- rownames(mg_EPIC_inm_sum) 

mg_MIXTURE_inm_sum <- MIXTURE(expressionMatrix = mg_counts_inm_sum,      
                             signatureMatrix = mg_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
mg_MIXTURE_inm_sum <- as.data.frame(mg_MIXTURE_inm_sum$Subjects$MIXprop)
mg_MIXTURE_inm_sum$Method <- "MIXTURE"
mg_MIXTURE_inm_sum$Signature <- "Mammary gland"
mg_MIXTURE_inm_sum$Celltype <- rownames(mg_MIXTURE_inm_sum) 

rownames(mg_sig) <- toupper(rownames(mg_sig))
write.table(data.frame("X1"=rownames(mg_sig),mg_sig), file = "Figure_5/Data/MS/MG_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(mg_counts_inm_sum) <- toupper(rownames(mg_counts_inm_sum))
mg_quanTIseq_inm_sum <-deconvolute_quantiseq.default(mix.mat = mg_counts_inm_sum, 
                                                    arrays = FALSE, 
                                                    signame = "Figure_5/Data/MS/MG_sig", 
                                                    tumor = FALSE, 
                                                    mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                    rmgenes = "unassigned")[, -c(10)]
mg_quanTIseq_inm_sum[,1] <- NULL
mg_quanTIseq_inm_sum$Method <- "quanTIseq"
mg_quanTIseq_inm_sum$Signature <- "Mammary gland"
mg_quanTIseq_inm_sum$Celltype <- rownames(mg_quanTIseq_inm_sum)
names(mg_quanTIseq_inm_sum) <- gsub("[.]"," ", names(mg_quanTIseq_inm_sum))

mg_inm_est <- rbind(mg_CIBERSORT_inm_sum,mg_EPIC_inm_sum,mg_MIXTURE_inm_sum,mg_quanTIseq_inm_sum)

m_sig <- read.csv("Figure_5/Data/MS/M.csv", row.names = 1, check.names=FALSE)

m_CIBERSORT_inm_sum <- read_csv("Figure_5/Data/CIBERSORT_estimations/m_CIBERSORT_inm_sum_prop.csv")
m_CIBERSORT_inm_sum <- as.data.frame(m_CIBERSORT_inm_sum)
rownames(m_CIBERSORT_inm_sum) <- m_CIBERSORT_inm_sum[,1]
m_CIBERSORT_inm_sum <- m_CIBERSORT_inm_sum[-c((nrow(m_CIBERSORT_inm_sum)-4):((nrow(m_CIBERSORT_inm_sum)))), -c(1,10,11,12,13)]
m_CIBERSORT_inm_sum$Method <- "CIBERSORT"
m_CIBERSORT_inm_sum$Signature <- "Muscle"
m_CIBERSORT_inm_sum$Celltype <- rownames(m_CIBERSORT_inm_sum) 

m_EPIC_inm_sum <- EPIC::EPIC(bulk = m_counts_inm_sum,
                             reference = list(refProfiles = m_sig, sigGenes = rownames(m_sig)),
                             withOtherCells = FALSE)
m_EPIC_inm_sum <- as.data.frame(m_EPIC_inm_sum$cellFractions)
m_EPIC_inm_sum$Method <- "EPIC"
m_EPIC_inm_sum$Signature <- "Muscle"
m_EPIC_inm_sum$Celltype <- rownames(m_EPIC_inm_sum) 

m_MIXTURE_inm_sum <- MIXTURE(expressionMatrix = m_counts_inm_sum,      
                             signatureMatrix = m_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
m_MIXTURE_inm_sum <- as.data.frame(m_MIXTURE_inm_sum$Subjects$MIXprop)
m_MIXTURE_inm_sum$Method <- "MIXTURE"
m_MIXTURE_inm_sum$Signature <- "Muscle"
m_MIXTURE_inm_sum$Celltype <- rownames(m_MIXTURE_inm_sum) 

rownames(m_sig) <- toupper(rownames(m_sig))
write.table(data.frame("X1"=rownames(m_sig),m_sig), file = "Figure_5/Data/MS/m_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(m_counts_inm_sum) <- toupper(rownames(m_counts_inm_sum))
m_quanTIseq_inm_sum <-deconvolute_quantiseq.default(mix.mat = m_counts_inm_sum, 
                                                    arrays = FALSE, 
                                                    signame = "Figure_5/Data/MS/m_sig", 
                                                    tumor = FALSE, 
                                                    mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                    rmgenes = "unassigned")[, -c(10)]
m_quanTIseq_inm_sum[,1] <- NULL
m_quanTIseq_inm_sum$Method <- "quanTIseq"
m_quanTIseq_inm_sum$Signature <- "Muscle"
m_quanTIseq_inm_sum$Celltype <- rownames(m_quanTIseq_inm_sum)
names(m_quanTIseq_inm_sum) <- gsub("[.]"," ", names(m_quanTIseq_inm_sum))

m_inm_est <- rbind(m_CIBERSORT_inm_sum,m_EPIC_inm_sum,m_MIXTURE_inm_sum,m_quanTIseq_inm_sum)

si_sig <- read.csv("Figure_5/Data/MS/SI.csv", row.names = 1, check.names=FALSE)

si_CIBERSORT_inm_sum_1 <- read_csv("Figure_5/Data/CIBERSORT_estimations/si_CIBERSORT_inm_sum_prop_1.csv")
si_CIBERSORT_inm_sum_1 <- as.data.frame(si_CIBERSORT_inm_sum_1)
rownames(si_CIBERSORT_inm_sum_1) <- si_CIBERSORT_inm_sum_1[,1]
si_CIBERSORT_inm_sum_1 <- si_CIBERSORT_inm_sum_1[-c((nrow(si_CIBERSORT_inm_sum_1)-5):((nrow(si_CIBERSORT_inm_sum_1)))), -c(1,10,11,12,13)]
si_CIBERSORT_inm_sum_1$Method <- "CIBERSORT"
si_CIBERSORT_inm_sum_1$Signature <- "Small intestine"
si_CIBERSORT_inm_sum_1$Celltype <- rownames(si_CIBERSORT_inm_sum_1) 

si_EPIC_inm_sum_1 <- EPIC::EPIC(bulk = si_counts_inm_sum_1,
                                reference = list(refProfiles = si_sig, sigGenes = rownames(si_sig)),
                                withOtherCells = FALSE)
si_EPIC_inm_sum_1 <- as.data.frame(si_EPIC_inm_sum_1$cellFractions)
si_EPIC_inm_sum_1$Method <- "EPIC"
si_EPIC_inm_sum_1$Signature <- "Small intestine"
si_EPIC_inm_sum_1$Celltype <- rownames(si_EPIC_inm_sum_1) 

si_MIXTURE_inm_sum_1 <- MIXTURE(expressionMatrix = si_counts_inm_sum_1,      
                                signatureMatrix = si_sig,       
                                iter = 0L,
                                functionMixture = nu.svm.robust.RFE,
                                useCores = 1,
                                verbose = TRUE,
                                nullDist = "PopulationBased") 
si_MIXTURE_inm_sum_1 <- as.data.frame(si_MIXTURE_inm_sum_1$Subjects$MIXprop)
si_MIXTURE_inm_sum_1$Method <- "MIXTURE"
si_MIXTURE_inm_sum_1$Signature <- "Small intestine"
si_MIXTURE_inm_sum_1$Celltype <- rownames(si_MIXTURE_inm_sum_1) 

rownames(si_sig) <- toupper(rownames(si_sig))
write.table(data.frame("X1"=rownames(si_sig),si_sig), file = "Figure_5/Data/MS/SI_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(si_counts_inm_sum_1) <- toupper(rownames(si_counts_inm_sum_1))
si_quanTIseq_inm_sum_1 <-deconvolute_quantiseq.default(mix.mat = si_counts_inm_sum_1, 
                                                       arrays = FALSE, 
                                                       signame = "Figure_5/Data/MS/SI_sig", 
                                                       tumor = FALSE, 
                                                       mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                       rmgenes = "unassigned")[, -c(10)]
si_quanTIseq_inm_sum_1[,1] <- NULL
si_quanTIseq_inm_sum_1$Method <- "quanTIseq"
si_quanTIseq_inm_sum_1$Signature <- "Small intestine"
si_quanTIseq_inm_sum_1$Celltype <- rownames(si_quanTIseq_inm_sum_1)
names(si_quanTIseq_inm_sum_1) <- gsub("[.]"," ", names(si_quanTIseq_inm_sum_1))

si_sig <- read.csv("Figure_5/Data/MS/SI.csv", row.names = 1, check.names=FALSE)

si_CIBERSORT_inm_sum_2 <- read_csv("Figure_5/Data/CIBERSORT_estimations/si_CIBERSORT_inm_sum_prop_2.csv")
si_CIBERSORT_inm_sum_2 <- as.data.frame(si_CIBERSORT_inm_sum_2)
rownames(si_CIBERSORT_inm_sum_2) <- si_CIBERSORT_inm_sum_2[,1]
si_CIBERSORT_inm_sum_2 <- si_CIBERSORT_inm_sum_2[-c((nrow(si_CIBERSORT_inm_sum_2)-4):((nrow(si_CIBERSORT_inm_sum_2)))), -c(1,10,11,12,13)]
si_CIBERSORT_inm_sum_2$Method <- "CIBERSORT"
si_CIBERSORT_inm_sum_2$Signature <- "Small intestine"
si_CIBERSORT_inm_sum_2$Celltype <- rownames(si_CIBERSORT_inm_sum_2) 

si_EPIC_inm_sum_2 <- EPIC::EPIC(bulk = si_counts_inm_sum_2,
                                reference = list(refProfiles = si_sig, sigGenes = rownames(si_sig)),
                                withOtherCells = FALSE)
si_EPIC_inm_sum_2 <- as.data.frame(si_EPIC_inm_sum_2$cellFractions)
si_EPIC_inm_sum_2$Method <- "EPIC"
si_EPIC_inm_sum_2$Signature <- "Small intestine"
si_EPIC_inm_sum_2$Celltype <- rownames(si_EPIC_inm_sum_2) 

si_MIXTURE_inm_sum_2 <- MIXTURE(expressionMatrix = si_counts_inm_sum_2,      
                                signatureMatrix = si_sig,       
                                iter = 0L,
                                functionMixture = nu.svm.robust.RFE,
                                useCores = 1,
                                verbose = TRUE,
                                nullDist = "PopulationBased") 
si_MIXTURE_inm_sum_2 <- as.data.frame(si_MIXTURE_inm_sum_2$Subjects$MIXprop)
si_MIXTURE_inm_sum_2$Method <- "MIXTURE"
si_MIXTURE_inm_sum_2$Signature <- "Small intestine"
si_MIXTURE_inm_sum_2$Celltype <- rownames(si_MIXTURE_inm_sum_2) 

rownames(si_counts_inm_sum_2) <- toupper(rownames(si_counts_inm_sum_2))
si_quanTIseq_inm_sum_2 <-deconvolute_quantiseq.default(mix.mat = si_counts_inm_sum_2, 
                                                       arrays = FALSE, 
                                                       signame = "Figure_5/Data/MS/SI_sig", 
                                                       tumor = FALSE, 
                                                       mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                       rmgenes = "unassigned")[, -c(10)]
si_quanTIseq_inm_sum_2[,1] <- NULL
si_quanTIseq_inm_sum_2$Method <- "quanTIseq"
si_quanTIseq_inm_sum_2$Signature <- "Small intestine"
si_quanTIseq_inm_sum_2$Celltype <- rownames(si_quanTIseq_inm_sum_2)
names(si_quanTIseq_inm_sum_2) <- gsub("[.]"," ", names(si_quanTIseq_inm_sum_2))

si_inm_est <- rbind(si_CIBERSORT_inm_sum_1,si_CIBERSORT_inm_sum_2,
                    si_EPIC_inm_sum_1,si_EPIC_inm_sum_2,
                    si_MIXTURE_inm_sum_1,si_MIXTURE_inm_sum_2,
                    si_quanTIseq_inm_sum_1,si_quanTIseq_inm_sum_2)

s_sig <- read.csv("Figure_5/Data/MS/S.csv", row.names = 1, check.names=FALSE)

s_CIBERSORT_inm_sum <- read_csv("Figure_5/Data/CIBERSORT_estimations/s_CIBERSORT_inm_sum_prop.csv")
s_CIBERSORT_inm_sum <- as.data.frame(s_CIBERSORT_inm_sum)
rownames(s_CIBERSORT_inm_sum) <- s_CIBERSORT_inm_sum[,1]
s_CIBERSORT_inm_sum <- s_CIBERSORT_inm_sum[-c((nrow(s_CIBERSORT_inm_sum)-2):((nrow(s_CIBERSORT_inm_sum)))), -c(1,10,11,12,13)]
s_CIBERSORT_inm_sum$Method <- "CIBERSORT"
s_CIBERSORT_inm_sum$Signature <- "Spleen"
s_CIBERSORT_inm_sum$Celltype <- rownames(s_CIBERSORT_inm_sum) 

s_EPIC_inm_sum <- EPIC::EPIC(bulk = s_counts_inm_sum,
                             reference = list(refProfiles = s_sig, sigGenes = rownames(s_sig)),
                             withOtherCells = FALSE)
s_EPIC_inm_sum <- as.data.frame(s_EPIC_inm_sum$cellFractions)
s_EPIC_inm_sum$Method <- "EPIC"
s_EPIC_inm_sum$Signature <- "Spleen"
s_EPIC_inm_sum$Celltype <- rownames(s_EPIC_inm_sum) 

s_MIXTURE_inm_sum <- MIXTURE(expressionMatrix = s_counts_inm_sum,      
                             signatureMatrix = s_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
s_MIXTURE_inm_sum <- as.data.frame(s_MIXTURE_inm_sum$Subjects$MIXprop)
s_MIXTURE_inm_sum$Method <- "MIXTURE"
s_MIXTURE_inm_sum$Signature <- "Spleen"
s_MIXTURE_inm_sum$Celltype <- rownames(s_MIXTURE_inm_sum) 

rownames(s_sig) <- toupper(rownames(s_sig))
write.table(data.frame("X1"=rownames(s_sig),s_sig), file = "Figure_5/Data/MS/s_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(s_counts_inm_sum) <- toupper(rownames(s_counts_inm_sum))
s_quanTIseq_inm_sum <-deconvolute_quantiseq.default(mix.mat = s_counts_inm_sum, 
                                                    arrays = FALSE, 
                                                    signame = "Figure_5/Data/MS/s_sig", 
                                                    tumor = FALSE, 
                                                    mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                                    rmgenes = "unassigned")[, -c(10)]
s_quanTIseq_inm_sum[,1] <- NULL
s_quanTIseq_inm_sum$Method <- "quanTIseq"
s_quanTIseq_inm_sum$Signature <- "Spleen"
s_quanTIseq_inm_sum$Celltype <- rownames(s_quanTIseq_inm_sum)
names(s_quanTIseq_inm_sum) <- gsub("[.]"," ", names(s_quanTIseq_inm_sum))

s_inm_est <- rbind(s_CIBERSORT_inm_sum,s_EPIC_inm_sum,s_MIXTURE_inm_sum,s_quanTIseq_inm_sum)

inm_est <- rbind(k_inm_est, l_inm_est, mg_inm_est, m_inm_est, si_inm_est, s_inm_est)
rownames(inm_est) <- NULL

#Table 1

inm_est_melt <- melt(inm_est)
inm_est_melt$Celltype <- gsub("Neutrophil progenitor", "Neutrophil", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Denditric cell", "Dendritic cell", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Marginal zone B cell", "B cell", inm_est_melt$Celltype)

Table1 <- inm_est_melt %>% mutate(result = case_when(Celltype == variable & value > 0 ~ "TP", 
                                                           Celltype != variable & value > 0 ~ "FP",
                                                           Celltype == variable & value == 0 ~ "FN",
                                                           Celltype != variable & value == 0 ~ "TN")) %>% 
  group_by(result,Method) %>% 
  dplyr::summarize(count=n()) %>% 
  ungroup() %>% 
  complete(result, nesting(Method), fill = list(count = 0)) %>% 
  dcast(formula = Method~result,fun.aggregate = sum,value.var = "count")

Table1 <- transform(Table1, Se = TP / (TP + FN))
Table1 <- transform(Table1, Sp = TN / (TN + FP))
Table1 <- transform(Table1, PPV = TP / (TP + FP))
Table1 <- transform(Table1, NPV = TN / (TN + FN))
Table1 <- transform(Table1, F1 = (2 * TP) / ((2 * TP) + FP + FN))
Table1 <- transform(Table1, DOP = sqrt(((Se - 1)^2)+((Sp - 1)^2)+((PPV - 1)^2)+((NPV - 1)^2)))

#Figures

inm_est_melt <- melt(inm_est)
inm_est_melt$Celltype <- gsub("Neutrophil progenitor", "Neutrophil", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Denditric cell", "Dendritic cell", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Marginal zone B cell", "B cell", inm_est_melt$Celltype)
inm_est_melt$value <- as.double(inm_est_melt$value)

inm_est_TP <- inm_est_melt %>% filter(Celltype == variable, value > 0)
inm_est_TP$value <- as.numeric(inm_est_TP$value)
inm_est_TP$Celltype <- as.factor(inm_est_TP$Celltype)

inm_est_TP$CT <- paste0(inm_est_TP$Signature,"-",inm_est_TP$Celltype)

n <- count(inm_est_TP, Method)
n$n <- paste0("n = ", n$n)

Figure_5a <-  ggplot(inm_est_TP, aes(factor(Method), value)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for TP") + scale_y_continuous(labels = scales::percent)

Sup_Figure_3 <- ggplot(inm_est_TP, aes(factor(Signature), value)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0, alpha = 0.5, height = 0)+
  labs(x = "Signature",
       y = "Predicted percentages for TP") +
  facet_wrap(~Method)

inm_est_FP <- inm_est_melt %>% filter(Celltype != variable, value > 0)
inm_est_FP$value <- as.numeric(inm_est_FP$value)
inm_est_FP$Celltype <- as.factor(inm_est_FP$Celltype)
inm_est_FP$variable <- as.factor(inm_est_FP$variable)

n1 <- count(inm_est_FP, Method)
n1$n <- paste0("n = ", n1$n)

Figure_5b <-  ggplot(inm_est_FP, aes(factor(Method), value)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n1, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for FP") + scale_y_continuous(labels = scales::percent)

inm_est_FP$Signature <- as.factor(inm_est_FP$Signature)

Sup_Figure_4 <- ggplot(inm_est_FP, aes(x = Signature, y =value)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  labs(x = "Signature",
       y = "Predicted percentage for FP") +
  facet_wrap(~ Method) + scale_y_continuous(labels = scales::percent)

TP_CM <- subset(inm_est_TP, Method %in% c("CIBERSORT","MIXTURE"))
n2 <- count(TP_CM, Method)
n2$n <- paste0("n = ", n2$n)

inm_est_melt <- melt(inm_est)
inm_est_melt$Celltype <- gsub("Neutrophil progenitor", "Neutrophil", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Denditric cell", "Dendritic cell", inm_est_melt$Celltype)
inm_est_melt$Celltype <- gsub("Marginal zone B cell", "B cell", inm_est_melt$Celltype)
inm_est_melt$value <- as.double(inm_est_melt$value)

inm_est_TP <- inm_est_melt %>% filter(Celltype == variable)
inm_est_TP$value <- as.numeric(inm_est_TP$value)
inm_est_TP$Celltype <- as.factor(inm_est_TP$Celltype)
TP_CM <- subset(inm_est_TP, Method %in% c("CIBERSORT","MIXTURE"))

Sup_Figure_2 <- ggpaired(TP_CM, x = "Method", y = "value",
                         line.color = "gray", line.size = 0.4,
                         palette = "jco", xlab ="Method", ylab = "Predicted percentage for TP")+
  stat_compare_means(paired = TRUE) + scale_y_continuous(labels = scales::percent)

Figure_5 <- ggarrange(Figure_5a, Figure_5b,
                      labels = c("(a)", "(b)"),
                      ncol = 1, nrow = 2)
