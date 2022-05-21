library(tidyverse)
library(immunedeconv)
library(MIXTURE)
library(ggplot2)
library(readr)
library(reshape2)

load('Figure_2-3/Data/MixSims_BetaSims.RData')

K_sig <- read.csv("Figure_2-3/Data/MS/K.csv", row.names = 1, check.names=FALSE)
K_betasim <- subset(K_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(K_betasim)[2] <- "simcelltypes"

K_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Kidney_CIBERSORT.txt")[,-c(10:12)]
names(K_CIBERSORT)[1] <- "ID"
names(K_CIBERSORT) <- gsub("[.]"," ", names(K_CIBERSORT))
K_CIBERSORT <- melt(K_CIBERSORT)
K_CIBERSORT <- subset(K_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(K_CIBERSORT)[2] <- "estcelltypes"
K_CIBERSORT$Signature <- "Kidney"
K_CIBERSORT$Method <- "CIBERSORT"

K_EPIC <- EPIC::EPIC(bulk = K_mixsim,
                     reference = list(refProfiles = K_sig, sigGenes = rownames(K_sig)),
                     withOtherCells = FALSE)
K_EPIC <- as.data.frame(K_EPIC$cellFractions)
K_EPIC$ID <- rownames(K_EPIC)
K_EPIC <- melt(K_EPIC, id = "ID")
K_EPIC <- subset(K_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(K_EPIC)[2] <- "estcelltypes"
K_EPIC$Signature <- "Kidney"
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
K_MIXTURE <- subset(K_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(K_MIXTURE)[2] <- "estcelltypes"
K_MIXTURE$Signature <- "Kidney"
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
K_quanTIseq <- subset(K_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(K_quanTIseq)[2] <- "estcelltypes"
K_quanTIseq$Signature <- "Kidney"
K_quanTIseq$Method <- "quanTIseq"

K_final <- rbind(K_CIBERSORT,K_EPIC,K_MIXTURE,K_quanTIseq)
K_final <- left_join(K_final,
                     K_betasim,
                     by = "ID")
K_final$ID <- NULL

L_sig <- read.csv("Figure_2-3/Data/MS/L.csv", row.names = 1, check.names=FALSE)
L_betasim <- subset(L_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(L_betasim)[2] <- "simcelltypes"

L_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Liver_CIBERSORT.txt")[,-c(10:12)]
names(L_CIBERSORT)[1] <- "ID"
names(L_CIBERSORT) <- gsub("[.]"," ", names(L_CIBERSORT))
L_CIBERSORT <- melt(L_CIBERSORT)
L_CIBERSORT <- subset(L_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(L_CIBERSORT)[2] <- "estcelltypes"
L_CIBERSORT$Signature <- "Liver"
L_CIBERSORT$Method <- "CIBERSORT"

L_EPIC <- EPIC::EPIC(bulk = L_mixsim,
                     reference = list(refProfiles = L_sig, sigGenes = rownames(L_sig)),
                     withOtherCells = FALSE)
L_EPIC <- as.data.frame(L_EPIC$cellFractions)
L_EPIC$ID <- rownames(L_EPIC)
L_EPIC <- melt(L_EPIC, id = "ID")
L_EPIC <- subset(L_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(L_EPIC)[2] <- "estcelltypes"
L_EPIC$Signature <- "Liver"
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
L_MIXTURE <- subset(L_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(L_MIXTURE)[2] <- "estcelltypes"
L_MIXTURE$Signature <- "Liver"
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
L_quanTIseq <- subset(L_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(L_quanTIseq)[2] <- "estcelltypes"
L_quanTIseq$Signature <- "Liver"
L_quanTIseq$Method <- "quanTIseq"

L_final <- rbind(L_CIBERSORT,L_EPIC,L_MIXTURE,L_quanTIseq)
L_final <- left_join(L_final,
                     L_betasim,
                     by = "ID")
L_final$ID <- NULL

MG_sig <- read.csv("Figure_2-3/Data/MS/MG.csv", row.names = 1, check.names=FALSE)
MG_betasim <- subset(MG_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(MG_betasim)[2] <- "simcelltypes"

MG_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/MammaryGland_CIBERSORT.txt")[,-c(10:12)]
names(MG_CIBERSORT)[1] <- "ID"
names(MG_CIBERSORT) <- gsub("[.]"," ", names(MG_CIBERSORT))
MG_CIBERSORT <- melt(MG_CIBERSORT)
MG_CIBERSORT <- subset(MG_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(MG_CIBERSORT)[2] <- "estcelltypes"
MG_CIBERSORT$Signature <- "Mammary gland"
MG_CIBERSORT$Method <- "CIBERSORT"

MG_EPIC <- EPIC::EPIC(bulk = MG_mixsim,
                     reference = list(refProfiles = MG_sig, sigGenes = rownames(MG_sig)),
                     withOtherCells = FALSE)
MG_EPIC <- as.data.frame(MG_EPIC$cellFractions)
MG_EPIC$ID <- rownames(MG_EPIC)
MG_EPIC <- melt(MG_EPIC, id = "ID")
MG_EPIC <- subset(MG_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(MG_EPIC)[2] <- "estcelltypes"
MG_EPIC$Signature <- "Mammary gland"
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
MG_MIXTURE <- subset(MG_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(MG_MIXTURE)[2] <- "estcelltypes"
MG_MIXTURE$Signature <- "Mammary gland"
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
MG_quanTIseq <- subset(MG_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(MG_quanTIseq)[2] <- "estcelltypes"
MG_quanTIseq$Signature <- "Mammary gland"
MG_quanTIseq$Method <- "quanTIseq"

MG_final <- rbind(MG_CIBERSORT,MG_EPIC,MG_MIXTURE,MG_quanTIseq)
MG_final <- left_join(MG_final,
                     MG_betasim,
                     by = "ID")
MG_final$ID <- NULL

M_sig <- read.csv("Figure_2-3/Data/MS/M.csv", row.names = 1, check.names=FALSE)
M_betasim <- subset(M_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(M_betasim)[2] <- "simcelltypes"

M_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Muscle_CIBERSORT.txt")[,-c(10:12)]
names(M_CIBERSORT)[1] <- "ID"
names(M_CIBERSORT) <- gsub("[.]"," ", names(M_CIBERSORT))
M_CIBERSORT <- melt(M_CIBERSORT)
M_CIBERSORT <- subset(M_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(M_CIBERSORT)[2] <- "estcelltypes"
M_CIBERSORT$Signature <- "Muscle"
M_CIBERSORT$Method <- "CIBERSORT"

M_EPIC <- EPIC::EPIC(bulk = M_mixsim,
                     reference = list(refProfiles = M_sig, sigGenes = rownames(M_sig)),
                     withOtherCells = FALSE)
M_EPIC <- as.data.frame(M_EPIC$cellFractions)
M_EPIC$ID <- rownames(M_EPIC)
M_EPIC <- melt(M_EPIC, id = "ID")
M_EPIC <- subset(M_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(M_EPIC)[2] <- "estcelltypes"
M_EPIC$Signature <- "Muscle"
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
M_MIXTURE <- subset(M_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(M_MIXTURE)[2] <- "estcelltypes"
M_MIXTURE$Signature <- "Muscle"
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
M_quanTIseq <- subset(M_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(M_quanTIseq)[2] <- "estcelltypes"
M_quanTIseq$Signature <- "Muscle"
M_quanTIseq$Method <- "quanTIseq"

M_final <- rbind(M_CIBERSORT,M_EPIC,M_MIXTURE,M_quanTIseq)
M_final <- left_join(M_final,
                     M_betasim,
                     by = "ID")
M_final$ID <- NULL

SI_sig <- read.csv("Figure_2-3/Data/MS/SI.csv", row.names = 1, check.names=FALSE)
SI_betasim <- subset(SI_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(SI_betasim)[2] <- "simcelltypes"

SI_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/SmallIntestine_CIBERSORT.txt")[,-c(10:12)]
names(SI_CIBERSORT)[1] <- "ID"
names(SI_CIBERSORT) <- gsub("[.]"," ", names(SI_CIBERSORT))
SI_CIBERSORT <- melt(SI_CIBERSORT)
SI_CIBERSORT <- subset(SI_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(SI_CIBERSORT)[2] <- "estcelltypes"
SI_CIBERSORT$Signature <- "Small intestine"
SI_CIBERSORT$Method <- "CIBERSORT"

SI_EPIC <- EPIC::EPIC(bulk = SI_mixsim,
                     reference = list(refProfiles = SI_sig, sigGenes = rownames(SI_sig)),
                     withOtherCells = FALSE)
SI_EPIC <- as.data.frame(SI_EPIC$cellFractions)
SI_EPIC$ID <- rownames(SI_EPIC)
SI_EPIC <- melt(SI_EPIC, id = "ID")
SI_EPIC <- subset(SI_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(SI_EPIC)[2] <- "estcelltypes"
SI_EPIC$Signature <- "Small intestine"
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
SI_MIXTURE <- subset(SI_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(SI_MIXTURE)[2] <- "estcelltypes"
SI_MIXTURE$Signature <- "Small intestine"
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
SI_quanTIseq <- subset(SI_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(SI_quanTIseq)[2] <- "estcelltypes"
SI_quanTIseq$Signature <- "Small intestine"
SI_quanTIseq$Method <- "quanTIseq"

SI_final <- rbind(SI_CIBERSORT,SI_EPIC,SI_MIXTURE,SI_quanTIseq)
SI_final <- left_join(SI_final,
                     SI_betasim,
                     by = "ID")
SI_final$ID <- NULL

S_sig <- read.csv("Figure_2-3/Data/MS/S.csv", row.names = 1, check.names=FALSE)
S_betasim <- subset(S_betasim, betasim > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(S_betasim)[2] <- "simcelltypes"

S_CIBERSORT <- read_table("Figure_2-3/Data/CIBERSORT_estimations/Spleen_CIBERSORT.txt")[,-c(10:12)]
names(S_CIBERSORT)[1] <- "ID"
names(S_CIBERSORT) <- gsub("[.]"," ", names(S_CIBERSORT))
S_CIBERSORT <- melt(S_CIBERSORT)
S_CIBERSORT <- subset(S_CIBERSORT, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(S_CIBERSORT)[2] <- "estcelltypes"
S_CIBERSORT$Signature <- "Spleen"
S_CIBERSORT$Method <- "CIBERSORT"

S_EPIC <- EPIC::EPIC(bulk = S_mixsim,
                     reference = list(refProfiles = S_sig, sigGenes = rownames(S_sig)),
                     withOtherCells = FALSE)
S_EPIC <- as.data.frame(S_EPIC$cellFractions)
S_EPIC$ID <- rownames(S_EPIC)
S_EPIC <- melt(S_EPIC, id = "ID")
S_EPIC <- subset(S_EPIC, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(S_EPIC)[2] <- "estcelltypes"
S_EPIC$Signature <- "Spleen"
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
S_MIXTURE <- subset(S_MIXTURE, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(S_MIXTURE)[2] <- "estcelltypes"
S_MIXTURE$Signature <- "Spleen"
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
S_quanTIseq <- subset(S_quanTIseq, value > 0 ) %>% 
  group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  ungroup()
names(S_quanTIseq)[2] <- "estcelltypes"
S_quanTIseq$Signature <- "Spleen"
S_quanTIseq$Method <- "quanTIseq"

S_final <- rbind(S_CIBERSORT,S_EPIC,S_MIXTURE,S_quanTIseq)
S_final <- left_join(S_final,
                     S_betasim,
                     by = "ID")
S_final$ID <- NULL

counts_FINAL <- rbind(K_final, L_final, MG_final, M_final, SI_final, S_final)

#Figure

counts_FINAL$simcelltypes <- as.factor(counts_FINAL$simcelltypes)

Figure_3 <- ggplot(counts_FINAL, aes(x=simcelltypes, y=estcelltypes)) +
  geom_boxplot(aes(colour = Method),show.legend=TRUE,outlier.shape = NA, fatten = NULL) +
  geom_jitter(width = 0.2, alpha = 0.1, height = 0) +
  labs( x = "True number of coefficients",y = "Estimated number of coefficients", fill = "Method") + 
  scale_y_continuous(name ="Estimated number of coefficients", breaks = c(0:22),labels=as.character(c(0:22))) +
  theme(legend.key = element_blank(),plot.title = element_text(hjust=0.5,size = 16),axis.line = element_line(colour = "black"),
        panel.background = element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 14),
        legend.position = "bottom") + 
  facet_grid(Method ~ Signature) 
