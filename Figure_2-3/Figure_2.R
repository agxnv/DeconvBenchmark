library(readr)
library(immunedeconv)
library(MIXTURE)
library(reshape2)
library(tidyverse)
library(ggpubr)


#Code for generating Simulated mixtures:
#coef <- as.integer((ncol(SM)) * 0.7)
#Simulation <- SimulatedMixtures(SM, coef, 500, noisy = TRUE)

load('Figure_2-3/Data/MixSims_BetaSims.RData')

K_sig <- read.csv("Figure_2-3/Data/MS/K.csv", row.names = 1, check.names=FALSE)

K_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Kidney_CIBERSORT.txt")[,-c(10:12)]
names(K_CIBERSORT)[1] <- "ID"
names(K_CIBERSORT) <- gsub("[.]"," ", names(K_CIBERSORT))
K_CIBERSORT <- melt(K_CIBERSORT)
K_CIBERSORT <- left_join(K_CIBERSORT,
                         K_betasim,
                         by = c("ID","variable"))
K_CIBERSORT[,c(1,2)] <- NULL
names(K_CIBERSORT)[1] <- "betahat"
K_CIBERSORT$difs <- K_CIBERSORT$betahat - K_CIBERSORT$betasim
K_CIBERSORT$Method <- "CIBERSORT"

K_EPIC <- EPIC::EPIC(bulk = K_mixsim,
                     reference = list(refProfiles = K_sig, sigGenes = rownames(K_sig)),
                     withOtherCells = FALSE)
K_EPIC <- as.data.frame(K_EPIC$cellFractions)
K_EPIC$ID <- rownames(K_EPIC)
K_EPIC <- melt(K_EPIC, id = "ID")
K_EPIC <- left_join(K_EPIC,
                    K_betasim,
                    by = c("ID","variable"))
K_EPIC[,c(1,2)] <- NULL
names(K_EPIC)[1] <- "betahat"
K_EPIC$difs <- K_EPIC$betahat - K_EPIC$betasim
K_EPIC$Method <- "EPIC"

K_MIXTURE <- MIXTURE(expressionMatrix = K_mixsim,      
                     signatureMatrix = K_sig,       
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE,
                     nullDist = "PopulationBased") 
K_MIXTURE <- as.data.frame(K_MIXTURE$Subjects$MIXprop)
K_MIXTURE$ID <- rownames(K_MIXTURE)
K_MIXTURE <- melt(K_MIXTURE, id = "ID")
K_MIXTURE <- left_join(K_MIXTURE,
                       K_betasim,
                       by = c("ID","variable"))
K_MIXTURE[,c(1,2)] <- NULL
names(K_MIXTURE)[1] <- "betahat"
K_MIXTURE$difs <- K_MIXTURE$betahat - K_MIXTURE$betasim
K_MIXTURE$Method <- "MIXTURE"

rownames(K_sig) <- toupper(rownames(K_sig))
write.table(data.frame("X1"=rownames(K_sig),K_sig), file = "Figure_2-3/Data/MS/K_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(K_mixsim) <- toupper(rownames(K_mixsim))
K_quanTIseq <-deconvolute_quantiseq.default(mix.mat = K_mixsim, 
                                            arrays = FALSE, 
                                            signame = "Figure_2-3/Data/MS/K_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
names(K_quanTIseq)[1] <- "ID"
names(K_quanTIseq) <- gsub("[.]"," ", names(K_quanTIseq))
K_quanTIseq <- melt(K_quanTIseq, id = "ID")
K_quanTIseq <- left_join(K_quanTIseq,
                         K_betasim,
                         by = c("ID","variable"))
K_quanTIseq[,c(1,2)] <- NULL
names(K_quanTIseq)[1] <- "betahat"
K_quanTIseq$difs <- K_quanTIseq$betahat - K_quanTIseq$betasim
K_quanTIseq$Method <- "quanTIseq"

K_final <- rbind(K_CIBERSORT,K_EPIC,K_MIXTURE,K_quanTIseq)

L_sig <- read.csv("Figure_2-3/Data/MS/L.csv", row.names = 1, check.names=FALSE)

L_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Liver_CIBERSORT.txt")[,-c(10:12)]
names(L_CIBERSORT)[1] <- "ID"
names(L_CIBERSORT) <- gsub("[.]"," ", names(L_CIBERSORT))
L_CIBERSORT <- melt(L_CIBERSORT)
L_CIBERSORT <- left_join(L_CIBERSORT,
                         L_betasim,
                         by = c("ID","variable"))
L_CIBERSORT[,c(1,2)] <- NULL
names(L_CIBERSORT)[1] <- "betahat"
L_CIBERSORT$difs <- L_CIBERSORT$betahat - L_CIBERSORT$betasim
L_CIBERSORT$Method <- "CIBERSORT"

L_EPIC <- EPIC::EPIC(bulk = L_mixsim,
                     reference = list(refProfiles = L_sig, sigGenes = rownames(L_sig)),
                     withOtherCells = FALSE)
L_EPIC <- as.data.frame(L_EPIC$cellFractions)
L_EPIC$ID <- rownames(L_EPIC)
L_EPIC <- melt(L_EPIC, id = "ID")
L_EPIC <- left_join(L_EPIC,
                    L_betasim,
                    by = c("ID","variable"))
L_EPIC[,c(1,2)] <- NULL
names(L_EPIC)[1] <- "betahat"
L_EPIC$difs <- L_EPIC$betahat - L_EPIC$betasim
L_EPIC$Method <- "EPIC"

L_MIXTURE <- MIXTURE(expressionMatrix = L_mixsim,      
                     signatureMatrix = L_sig,       
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE,
                     nullDist = "PopulationBased") 
L_MIXTURE <- as.data.frame(L_MIXTURE$Subjects$MIXprop)
L_MIXTURE$ID <- rownames(L_MIXTURE)
L_MIXTURE <- melt(L_MIXTURE, id = "ID")
L_MIXTURE <- left_join(L_MIXTURE,
                       L_betasim,
                       by = c("ID","variable"))
L_MIXTURE[,c(1,2)] <- NULL
names(L_MIXTURE)[1] <- "betahat"
L_MIXTURE$difs <- L_MIXTURE$betahat - L_MIXTURE$betasim
L_MIXTURE$Method <- "MIXTURE"

rownames(L_sig) <- toupper(rownames(L_sig))
write.table(data.frame("X1"=rownames(L_sig),L_sig), file = "Figure_2-3/Data/MS/L_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(L_mixsim) <- toupper(rownames(L_mixsim))
L_quanTIseq <-deconvolute_quantiseq.default(mix.mat = L_mixsim, 
                                            arrays = FALSE, 
                                            signame = "Figure_2-3/Data/MS/L_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
names(L_quanTIseq)[1] <- "ID"
names(L_quanTIseq) <- gsub("[.]"," ", names(L_quanTIseq))
L_quanTIseq <- melt(L_quanTIseq, id = "ID")
L_quanTIseq <- left_join(L_quanTIseq,
                         L_betasim,
                         by = c("ID","variable"))
L_quanTIseq[,c(1,2)] <- NULL
names(L_quanTIseq)[1] <- "betahat"
L_quanTIseq$difs <- L_quanTIseq$betahat - L_quanTIseq$betasim
L_quanTIseq$Method <- "quanTIseq"

L_final <- rbind(L_CIBERSORT,L_EPIC,L_MIXTURE,L_quanTIseq)

MG_sig <- read.csv("Figure_2-3/Data/MS/MG.csv", row.names = 1, check.names=FALSE)

MG_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/MammaryGland_CIBERSORT.txt")[,-c(10:12)]
names(MG_CIBERSORT)[1] <- "ID"
names(MG_CIBERSORT) <- gsub("[.]"," ", names(MG_CIBERSORT))
MG_CIBERSORT <- melt(MG_CIBERSORT)
MG_CIBERSORT <- left_join(MG_CIBERSORT,
                          MG_betasim,
                          by = c("ID","variable"))
MG_CIBERSORT[,c(1,2)] <- NULL
names(MG_CIBERSORT)[1] <- "betahat"
MG_CIBERSORT$difs <- MG_CIBERSORT$betahat - MG_CIBERSORT$betasim
MG_CIBERSORT$Method <- "CIBERSORT"

MG_EPIC <- EPIC::EPIC(bulk = MG_mixsim,
                      reference = list(refProfiles = MG_sig, sigGenes = rownames(MG_sig)),
                      withOtherCells = FALSE)
MG_EPIC <- as.data.frame(MG_EPIC$cellFractions)
MG_EPIC$ID <- rownames(MG_EPIC)
MG_EPIC <- melt(MG_EPIC, id = "ID")
MG_EPIC <- left_join(MG_EPIC,
                     MG_betasim,
                     by = c("ID","variable"))
MG_EPIC[,c(1,2)] <- NULL
names(MG_EPIC)[1] <- "betahat"
MG_EPIC$difs <- MG_EPIC$betahat - MG_EPIC$betasim
MG_EPIC$Method <- "EPIC"

MG_MIXTURE <- MIXTURE(expressionMatrix = MG_mixsim,      
                      signatureMatrix = MG_sig,       
                      iter = 0L,
                      functionMixture = nu.svm.robust.RFE,
                      useCores = 1,
                      verbose = TRUE,
                      nullDist = "PopulationBased") 
MG_MIXTURE <- as.data.frame(MG_MIXTURE$Subjects$MIXprop)
MG_MIXTURE$ID <- rownames(MG_MIXTURE)
MG_MIXTURE <- melt(MG_MIXTURE, id = "ID")
MG_MIXTURE <- left_join(MG_MIXTURE,
                        MG_betasim,
                        by = c("ID","variable"))
MG_MIXTURE[,c(1,2)] <- NULL
names(MG_MIXTURE)[1] <- "betahat"
MG_MIXTURE$difs <- MG_MIXTURE$betahat - MG_MIXTURE$betasim
MG_MIXTURE$Method <- "MIXTURE"

rownames(MG_sig) <- toupper(rownames(MG_sig))
write.table(data.frame("X1"=rownames(MG_sig),MG_sig), file = "Figure_2-3/Data/MS/MG_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(MG_mixsim) <- toupper(rownames(MG_mixsim))
MG_quanTIseq <-deconvolute_quantiseq.default(mix.mat = MG_mixsim, 
                                             arrays = FALSE, 
                                             signame = "Figure_2-3/Data/MS/MG_sig", 
                                             tumor = FALSE, 
                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                             rmgenes = "unassigned")[, -c(10)]
names(MG_quanTIseq)[1] <- "ID"
names(MG_quanTIseq) <- gsub("[.]"," ", names(MG_quanTIseq))
MG_quanTIseq <- melt(MG_quanTIseq, id = "ID")
MG_quanTIseq <- left_join(MG_quanTIseq,
                          MG_betasim,
                          by = c("ID","variable"))
MG_quanTIseq[,c(1,2)] <- NULL
names(MG_quanTIseq)[1] <- "betahat"
MG_quanTIseq$difs <- MG_quanTIseq$betahat - MG_quanTIseq$betasim
MG_quanTIseq$Method <- "quanTIseq"

MG_final <- rbind(MG_CIBERSORT,MG_EPIC,MG_MIXTURE,MG_quanTIseq)

M_sig <- read.csv("Figure_2-3/Data/MS/M.csv", row.names = 1, check.names=FALSE)

M_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Muscle_CIBERSORT.txt")[,-c(10:12)]
names(M_CIBERSORT)[1] <- "ID"
names(M_CIBERSORT) <- gsub("[.]"," ", names(M_CIBERSORT))
M_CIBERSORT <- melt(M_CIBERSORT)
M_CIBERSORT <- left_join(M_CIBERSORT,
                         M_betasim,
                         by = c("ID","variable"))
M_CIBERSORT[,c(1,2)] <- NULL
names(M_CIBERSORT)[1] <- "betahat"
M_CIBERSORT$difs <- M_CIBERSORT$betahat - M_CIBERSORT$betasim
M_CIBERSORT$Method <- "CIBERSORT"

M_EPIC <- EPIC::EPIC(bulk = M_mixsim,
                     reference = list(refProfiles = M_sig, sigGenes = rownames(M_sig)),
                     withOtherCells = FALSE)
M_EPIC <- as.data.frame(M_EPIC$cellFractions)
M_EPIC$ID <- rownames(M_EPIC)
M_EPIC <- melt(M_EPIC, id = "ID")
M_EPIC <- left_join(M_EPIC,
                    M_betasim,
                    by = c("ID","variable"))
M_EPIC[,c(1,2)] <- NULL
names(M_EPIC)[1] <- "betahat"
M_EPIC$difs <- M_EPIC$betahat - M_EPIC$betasim
M_EPIC$Method <- "EPIC"

M_MIXTURE <- MIXTURE(expressionMatrix = M_mixsim,      
                     signatureMatrix = M_sig,       
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE,
                     nullDist = "PopulationBased") 
M_MIXTURE <- as.data.frame(M_MIXTURE$Subjects$MIXprop)
M_MIXTURE$ID <- rownames(M_MIXTURE)
M_MIXTURE <- melt(M_MIXTURE, id = "ID")
M_MIXTURE <- left_join(M_MIXTURE,
                       M_betasim,
                       by = c("ID","variable"))
M_MIXTURE[,c(1,2)] <- NULL
names(M_MIXTURE)[1] <- "betahat"
M_MIXTURE$difs <- M_MIXTURE$betahat - M_MIXTURE$betasim
M_MIXTURE$Method <- "MIXTURE"

rownames(M_sig) <- toupper(rownames(M_sig))
write.table(data.frame("X1"=rownames(M_sig),M_sig), file = "Figure_2-3/Data/MS/M_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(M_mixsim) <- toupper(rownames(M_mixsim))
M_quanTIseq <-deconvolute_quantiseq.default(mix.mat = M_mixsim, 
                                            arrays = FALSE, 
                                            signame = "Figure_2-3/Data/MS/M_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
names(M_quanTIseq)[1] <- "ID"
names(M_quanTIseq) <- gsub("[.]"," ", names(M_quanTIseq))
M_quanTIseq <- melt(M_quanTIseq, id = "ID")
M_quanTIseq <- left_join(M_quanTIseq,
                         M_betasim,
                         by = c("ID","variable"))
M_quanTIseq[,c(1,2)] <- NULL
names(M_quanTIseq)[1] <- "betahat"
M_quanTIseq$difs <- M_quanTIseq$betahat - M_quanTIseq$betasim
M_quanTIseq$Method <- "quanTIseq"

M_final <- rbind(M_CIBERSORT,M_EPIC,M_MIXTURE,M_quanTIseq)

SI_sig <- read.csv("Figure_2-3/Data/MS/SI.csv", row.names = 1, check.names=FALSE)

SI_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/SmallIntestine_CIBERSORT.txt")[,-c(10:12)]
names(SI_CIBERSORT)[1] <- "ID"
names(SI_CIBERSORT) <- gsub("[.]"," ", names(SI_CIBERSORT))
SI_CIBERSORT <- melt(SI_CIBERSORT)
SI_CIBERSORT <- left_join(SI_CIBERSORT,
                          SI_betasim,
                          by = c("ID","variable"))
SI_CIBERSORT[,c(1,2)] <- NULL
names(SI_CIBERSORT)[1] <- "betahat"
SI_CIBERSORT$difs <- SI_CIBERSORT$betahat - SI_CIBERSORT$betasim
SI_CIBERSORT$Method <- "CIBERSORT"

SI_EPIC <- EPIC::EPIC(bulk = SI_mixsim,
                      reference = list(refProfiles = SI_sig, sigGenes = rownames(SI_sig)),
                      withOtherCells = FALSE)
SI_EPIC <- as.data.frame(SI_EPIC$cellFractions)
SI_EPIC$ID <- rownames(SI_EPIC)
SI_EPIC <- melt(SI_EPIC, id = "ID")
SI_EPIC <- left_join(SI_EPIC,
                     SI_betasim,
                     by = c("ID","variable"))
SI_EPIC[,c(1,2)] <- NULL
names(SI_EPIC)[1] <- "betahat"
SI_EPIC$difs <- SI_EPIC$betahat - SI_EPIC$betasim
SI_EPIC$Method <- "EPIC"

SI_MIXTURE <- MIXTURE(expressionMatrix = SI_mixsim,      
                      signatureMatrix = SI_sig,       
                      iter = 0L,
                      functionMixture = nu.svm.robust.RFE,
                      useCores = 1,
                      verbose = TRUE,
                      nullDist = "PopulationBased") 
SI_MIXTURE <- as.data.frame(SI_MIXTURE$Subjects$MIXprop)
SI_MIXTURE$ID <- rownames(SI_MIXTURE)
SI_MIXTURE <- melt(SI_MIXTURE, id = "ID")
SI_MIXTURE <- left_join(SI_MIXTURE,
                        SI_betasim,
                        by = c("ID","variable"))
SI_MIXTURE[,c(1,2)] <- NULL
names(SI_MIXTURE)[1] <- "betahat"
SI_MIXTURE$difs <- SI_MIXTURE$betahat - SI_MIXTURE$betasim
SI_MIXTURE$Method <- "MIXTURE"

rownames(SI_sig) <- toupper(rownames(SI_sig))
write.table(data.frame("X1"=rownames(SI_sig),SI_sig), file = "Figure_2-3/Data/MS/SI_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(SI_mixsim) <- toupper(rownames(SI_mixsim))
SI_quanTIseq <-deconvolute_quantiseq.default(mix.mat = SI_mixsim, 
                                             arrays = FALSE, 
                                             signame = "Figure_2-3/Data/MS/SI_sig", 
                                             tumor = FALSE, 
                                             mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                             rmgenes = "unassigned")[, -c(10)]
names(SI_quanTIseq)[1] <- "ID"
names(SI_quanTIseq) <- gsub("[.]"," ", names(SI_quanTIseq))
SI_quanTIseq <- melt(SI_quanTIseq, id = "ID")
SI_quanTIseq <- left_join(SI_quanTIseq,
                          SI_betasim,
                          by = c("ID","variable"))
SI_quanTIseq[,c(1,2)] <- NULL
names(SI_quanTIseq)[1] <- "betahat"
SI_quanTIseq$difs <- SI_quanTIseq$betahat - SI_quanTIseq$betasim
SI_quanTIseq$Method <- "quanTIseq"

SI_final <- rbind(SI_CIBERSORT,SI_EPIC,SI_MIXTURE,SI_quanTIseq)

S_sig <- read.csv("Figure_2-3/Data/MS/S.csv", row.names = 1, check.names=FALSE)

S_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Spleen_CIBERSORT.txt")[,-c(10:12)]
names(S_CIBERSORT)[1] <- "ID"
names(S_CIBERSORT) <- gsub("[.]"," ", names(S_CIBERSORT))
S_CIBERSORT <- melt(S_CIBERSORT)
S_CIBERSORT <- left_join(S_CIBERSORT,
                         S_betasim,
                         by = c("ID","variable"))
S_CIBERSORT[,c(1,2)] <- NULL
names(S_CIBERSORT)[1] <- "betahat"
S_CIBERSORT$difs <- S_CIBERSORT$betahat - S_CIBERSORT$betasim
S_CIBERSORT$Method <- "CIBERSORT"

S_EPIC <- EPIC::EPIC(bulk = S_mixsim,
                     reference = list(refProfiles = S_sig, sigGenes = rownames(S_sig)),
                     withOtherCells = FALSE)
S_EPIC <- as.data.frame(S_EPIC$cellFractions)
S_EPIC$ID <- rownames(S_EPIC)
S_EPIC <- melt(S_EPIC, id = "ID")
S_EPIC <- left_join(S_EPIC,
                    S_betasim,
                    by = c("ID","variable"))
S_EPIC[,c(1,2)] <- NULL
names(S_EPIC)[1] <- "betahat"
S_EPIC$difs <- S_EPIC$betahat - S_EPIC$betasim
S_EPIC$Method <- "EPIC"

S_MIXTURE <- MIXTURE(expressionMatrix = S_mixsim,      
                     signatureMatrix = S_sig,       
                     iter = 0L,
                     functionMixture = nu.svm.robust.RFE,
                     useCores = 1,
                     verbose = TRUE,
                     nullDist = "PopulationBased") 
S_MIXTURE <- as.data.frame(S_MIXTURE$Subjects$MIXprop)
S_MIXTURE$ID <- rownames(S_MIXTURE)
S_MIXTURE <- melt(S_MIXTURE, id = "ID")
S_MIXTURE <- left_join(S_MIXTURE,
                       S_betasim,
                       by = c("ID","variable"))
S_MIXTURE[,c(1,2)] <- NULL
names(S_MIXTURE)[1] <- "betahat"
S_MIXTURE$difs <- S_MIXTURE$betahat - S_MIXTURE$betasim
S_MIXTURE$Method <- "MIXTURE"

rownames(S_sig) <- toupper(rownames(S_sig))
write.table(data.frame("X1"=rownames(S_sig),S_sig), file = "Figure_2-3/Data/MS/S_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)
rownames(S_mixsim) <- toupper(rownames(S_mixsim))
S_quanTIseq <-deconvolute_quantiseq.default(mix.mat = S_mixsim, 
                                            arrays = FALSE, 
                                            signame = "Figure_2-3/Data/MS/S_sig", 
                                            tumor = FALSE, 
                                            mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                            rmgenes = "unassigned")[, -c(10)]
names(S_quanTIseq)[1] <- "ID"
names(S_quanTIseq) <- gsub("[.]"," ", names(S_quanTIseq))
S_quanTIseq <- melt(S_quanTIseq, id = "ID")
S_quanTIseq <- left_join(S_quanTIseq,
                         S_betasim,
                         by = c("ID","variable"))
S_quanTIseq[,c(1,2)] <- NULL
names(S_quanTIseq)[1] <- "betahat"
S_quanTIseq$difs <- S_quanTIseq$betahat - S_quanTIseq$betasim
S_quanTIseq$Method <- "quanTIseq"

S_final <- rbind(S_CIBERSORT,S_EPIC,S_MIXTURE,S_quanTIseq)

source('Figure_2-3/Data/BlandAltman.R')
K_BA <- BA.plot(K_final,title = "Kidney",graph = "Kidney")
L_BA <- BA.plot(L_final,title = "Liver",graph = "Liver")
MG_BA <- BA.plot(MG_final,title = "Mammary gland",graph = "Mammary gland")
M_BA <- BA.plot(M_final,title = "Muscle",graph = "Muscle")
SI_BA <- BA.plot(SI_final,title = "Small intestine",graph = "Small intestine")
S_BA <- BA.plot(S_final,title = "Spleen",graph = "Spleen")

ggsave(filename = "Figure_2.pdf",plot=grid.arrange(K_BA, L_BA, MG_BA, M_BA, SI_BA,S_BA,
                                                   ncol = 2, nrow = 3),device = "pdf",width = unit(24,"cm"),height = unit(20,"cm"))

#Sup_Table2-3

K_final$Signature <- "Kidney"
L_final$Signature <- "Liver"
MG_final$Signature <- "Mammary gland"
M_final$Signature <- "Muscle"
SI_final$Signature <- "Small intestine"
S_final$Signature <- "Spleen"

Sup_Table2 <- rbind(K_final,L_final,MG_final,M_final,SI_final,S_final) %>%
  mutate(result = case_when(betahat > 0 & betasim > 0 ~ "TP", 
                            betahat > 0 & betasim == 0 ~ "FP",
                            betahat == 0 & betasim == 0 ~ "TN",
                            betahat == 0 & betasim > 0 ~ "FN"
                            )) %>% 
  group_by(result,Signature,Method) %>% 
  dplyr::summarize(count=n()) %>% 
  ungroup() %>% 
  complete(result, nesting(Signature,Method), fill = list(count = 0))

Sup_Table2$MethodSig <- paste0(Sup_Table2$Method,"-",Sup_Table2$Signature)
Sup_Table2 <- dcast(data = Sup_Table2,formula = MethodSig~result,fun.aggregate = sum,value.var = "count")

Sup_Table2 <- transform(Sup_Table2, Se = TP / (TP + FN))
Sup_Table2 <- transform(Sup_Table2, Sp = TN / (TN + FP))
Sup_Table2 <- transform(Sup_Table2, PPV = TP / (TP + FP))
Sup_Table2 <- transform(Sup_Table2, NPV = TN / (TN + FN))
Sup_Table2 <- transform(Sup_Table2, F1 = (2 * TP) / ((2 * TP) + FP + FN))
Sup_Table2 <- transform(Sup_Table2, DOP = sqrt(((Se - 1)^2)+((Sp - 1)^2)+((PPV - 1)^2)+((NPV - 1)^2)))
Sup_Table2 <- transform(Sup_Table2, ER = (FP + FN)/(FP + FN + TP + TN))

#define quantiles of interest
q = c(.25, .5, .75)

Sup_Table3 <- rbind(K_final,L_final,MG_final,M_final,SI_final,S_final) %>%
  group_by(Method,Signature) %>%
  dplyr::summarize(quant25 = quantile(difs, probs = q[1]), 
            quant50 = quantile(difs, probs = q[2]),
            quant75 = quantile(difs, probs = q[3]),
            min = min(difs),
            max = max(difs),
            mean = mean(difs),
            sd = sd(difs),
            range = max(difs) - min(difs),
            correlation = cor(betahat,betasim)) %>% 
  ungroup()
