library(MIXTURE)
library(readr)
source('Figure_4/Data/quanTIseq.EPIC.abs.R')
load("Figure_4/Data/simulated_nonimmune.RData")

k_sig <- read.csv("Figure_4/Data/MS/K.csv", row.names = 1, check.names=FALSE)

k_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/K_CIBERSORT_par_sum.csv")
k_CIBERSORT_par_sum <- as.data.frame(k_CIBERSORT_par_sum)
rownames(k_CIBERSORT_par_sum) <- k_CIBERSORT_par_sum[,1]
k_CIBERSORT_par_sum <- k_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
k_CIBERSORT_par_sum$Signature <- "Kidney"
k_CIBERSORT_par_sum$Method <- "CIBERSORT"

k_EPIC_par_sum <- MyEPIC(bulk = k_counts_par_sum,
                         reference =  list(refProfiles=k_sig, sigGenes = rownames(k_sig)),
                         withOtherCells = FALSE)
k_EPIC_par_sum <- k_EPIC_par_sum$absCellFractions
k_EPIC_par_sum <- as.data.frame(k_EPIC_par_sum)
k_EPIC_par_sum$Signature <- "Kidney"
k_EPIC_par_sum$Method <- "EPIC"

k_MIXTURE_par_sum <- MIXTURE(expressionMatrix = k_counts_par_sum,      
                             signatureMatrix = k_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
k_MIXTURE_par_sum <- as.data.frame(k_MIXTURE_par_sum$Subjects$MIXabs)
k_MIXTURE_par_sum$Signature <- "Kidney"
k_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(k_counts_par_sum) <- toupper(rownames(k_counts_par_sum))
rownames(k_sig) <- toupper(rownames(k_sig))
write.table(data.frame("X1"=rownames(k_sig),k_sig), file = "Figure_4/Data/MS/k_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

k_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = k_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/k_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
k_quanTIseq_par_sum <- k_quanTIseq_par_sum$Abs[,-ncol(k_quanTIseq_par_sum$Abs)]
k_quanTIseq_par_sum <- as.data.frame(t(k_quanTIseq_par_sum))
k_quanTIseq_par_sum$Signature <- "Kidney"
k_quanTIseq_par_sum$Method <- "quanTIseq"
names(k_quanTIseq_par_sum) <- gsub("[.]"," ", names(k_quanTIseq_par_sum))

k_null_est <- rbind(k_CIBERSORT_par_sum,k_EPIC_par_sum,k_MIXTURE_par_sum,k_quanTIseq_par_sum)

l_sig <- read.csv("Figure_4/Data/MS/L.csv", row.names = 1, check.names=FALSE)

l_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/l_CIBERSORT_par_sum.csv")
l_CIBERSORT_par_sum <- as.data.frame(l_CIBERSORT_par_sum)
rownames(l_CIBERSORT_par_sum) <- l_CIBERSORT_par_sum[,1]
l_CIBERSORT_par_sum <- l_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
l_CIBERSORT_par_sum$Signature <- "Liver"
l_CIBERSORT_par_sum$Method <- "CIBERSORT"

l_EPIC_par_sum <- MyEPIC(bulk = l_counts_par_sum,
                         reference =  list(refProfiles=l_sig, sigGenes = rownames(l_sig)),
                         withOtherCells = FALSE)
l_EPIC_par_sum <- l_EPIC_par_sum$absCellFractions
l_EPIC_par_sum <- as.data.frame(l_EPIC_par_sum)
l_EPIC_par_sum$Signature <- "Liver"
l_EPIC_par_sum$Method <- "EPIC"

l_MIXTURE_par_sum <- MIXTURE(expressionMatrix = l_counts_par_sum,      
                             signatureMatrix = l_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
l_MIXTURE_par_sum <- as.data.frame(l_MIXTURE_par_sum$Subjects$MIXabs)
l_MIXTURE_par_sum$Signature <- "Liver"
l_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(l_counts_par_sum) <- toupper(rownames(l_counts_par_sum))
rownames(l_sig) <- toupper(rownames(l_sig))
write.table(data.frame("X1"=rownames(l_sig),l_sig), file = "Figure_4/Data/MS/l_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

l_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = l_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/l_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
l_quanTIseq_par_sum <- l_quanTIseq_par_sum$Abs[,-ncol(l_quanTIseq_par_sum$Abs)]
l_quanTIseq_par_sum <- as.data.frame(t(l_quanTIseq_par_sum))
l_quanTIseq_par_sum$Signature <- "Liver"
l_quanTIseq_par_sum$Method <- "quanTIseq"
names(l_quanTIseq_par_sum) <- gsub("[.]"," ", names(l_quanTIseq_par_sum))

l_null_est <- rbind(l_CIBERSORT_par_sum,l_EPIC_par_sum,l_MIXTURE_par_sum,l_quanTIseq_par_sum)

mg_sig <- read.csv("Figure_4/Data/MS/MG.csv", row.names = 1, check.names=FALSE)

mg_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/mg_CIBERSORT_par_sum.csv")
mg_CIBERSORT_par_sum <- as.data.frame(mg_CIBERSORT_par_sum)
rownames(mg_CIBERSORT_par_sum) <- mg_CIBERSORT_par_sum[,1]
mg_CIBERSORT_par_sum <- mg_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
mg_CIBERSORT_par_sum$Signature <- "Mammary gland"
mg_CIBERSORT_par_sum$Method <- "CIBERSORT"

mg_EPIC_par_sum <- MyEPIC(bulk = mg_counts_par_sum,
                         reference =  list(refProfiles=mg_sig, sigGenes = rownames(mg_sig)),
                         withOtherCells = FALSE)
mg_EPIC_par_sum <- mg_EPIC_par_sum$absCellFractions
mg_EPIC_par_sum <- as.data.frame(mg_EPIC_par_sum)
mg_EPIC_par_sum$Signature <- "Mammary gland"
mg_EPIC_par_sum$Method <- "EPIC"

mg_MIXTURE_par_sum <- MIXTURE(expressionMatrix = mg_counts_par_sum,      
                             signatureMatrix = mg_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
mg_MIXTURE_par_sum <- as.data.frame(mg_MIXTURE_par_sum$Subjects$MIXabs)
mg_MIXTURE_par_sum$Signature <- "Mammary gland"
mg_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(mg_counts_par_sum) <- toupper(rownames(mg_counts_par_sum))
rownames(mg_sig) <- toupper(rownames(mg_sig))
write.table(data.frame("X1"=rownames(mg_sig),mg_sig), file = "Figure_4/Data/MS/mg_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

mg_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = mg_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/mg_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
mg_quanTIseq_par_sum <- mg_quanTIseq_par_sum$Abs[,-ncol(mg_quanTIseq_par_sum$Abs)]
mg_quanTIseq_par_sum <- as.data.frame(t(mg_quanTIseq_par_sum))
mg_quanTIseq_par_sum$Signature <- "Mammary gland"
mg_quanTIseq_par_sum$Method <- "quanTIseq"
names(mg_quanTIseq_par_sum) <- gsub("[.]"," ", names(mg_quanTIseq_par_sum))

mg_null_est <- rbind(mg_CIBERSORT_par_sum,mg_EPIC_par_sum,mg_MIXTURE_par_sum,mg_quanTIseq_par_sum)

m_sig <- read.csv("Figure_4/Data/MS/M.csv", row.names = 1, check.names=FALSE)

m_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/m_CIBERSORT_par_sum.csv")
m_CIBERSORT_par_sum <- as.data.frame(m_CIBERSORT_par_sum)
rownames(m_CIBERSORT_par_sum) <- m_CIBERSORT_par_sum[,1]
m_CIBERSORT_par_sum <- m_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
m_CIBERSORT_par_sum$Signature <- "Muscle"
m_CIBERSORT_par_sum$Method <- "CIBERSORT"

m_EPIC_par_sum <- MyEPIC(bulk = m_counts_par_sum,
                         reference =  list(refProfiles=m_sig, sigGenes = rownames(m_sig)),
                         withOtherCells = FALSE)
m_EPIC_par_sum <- m_EPIC_par_sum$absCellFractions
m_EPIC_par_sum <- as.data.frame(m_EPIC_par_sum)
m_EPIC_par_sum$Signature <- "Muscle"
m_EPIC_par_sum$Method <- "EPIC"

m_MIXTURE_par_sum <- MIXTURE(expressionMatrix = m_counts_par_sum,      
                             signatureMatrix = m_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
m_MIXTURE_par_sum <- as.data.frame(m_MIXTURE_par_sum$Subjects$MIXabs)
m_MIXTURE_par_sum$Signature <- "Muscle"
m_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(m_counts_par_sum) <- toupper(rownames(m_counts_par_sum))
rownames(m_sig) <- toupper(rownames(m_sig))
write.table(data.frame("X1"=rownames(m_sig),m_sig), file = "Figure_4/Data/MS/m_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

m_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = m_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/m_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
m_quanTIseq_par_sum <- m_quanTIseq_par_sum$Abs[,-ncol(m_quanTIseq_par_sum$Abs)]
m_quanTIseq_par_sum <- as.data.frame(t(m_quanTIseq_par_sum))
m_quanTIseq_par_sum$Signature <- "Muscle"
m_quanTIseq_par_sum$Method <- "quanTIseq"
names(m_quanTIseq_par_sum) <- gsub("[.]"," ", names(m_quanTIseq_par_sum))

m_null_est <- rbind(m_CIBERSORT_par_sum,m_EPIC_par_sum,m_MIXTURE_par_sum,m_quanTIseq_par_sum)

si_sig <- read.csv("Figure_4/Data/MS/SI.csv", row.names = 1, check.names=FALSE)

si_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/si_CIBERSORT_par_sum.csv")
si_CIBERSORT_par_sum <- as.data.frame(si_CIBERSORT_par_sum)
rownames(si_CIBERSORT_par_sum) <- si_CIBERSORT_par_sum[,1]
si_CIBERSORT_par_sum <- si_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
si_CIBERSORT_par_sum$Signature <- "Small intestine"
si_CIBERSORT_par_sum$Method <- "CIBERSORT"

si_EPIC_par_sum <- MyEPIC(bulk = si_counts_par_sum,
                         reference =  list(refProfiles=si_sig, sigGenes = rownames(si_sig)),
                         withOtherCells = FALSE)
si_EPIC_par_sum <- si_EPIC_par_sum$absCellFractions
si_EPIC_par_sum <- as.data.frame(si_EPIC_par_sum)
si_EPIC_par_sum$Signature <- "Small intestine"
si_EPIC_par_sum$Method <- "EPIC"

si_MIXTURE_par_sum <- MIXTURE(expressionMatrix = si_counts_par_sum,      
                             signatureMatrix = si_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
si_MIXTURE_par_sum <- as.data.frame(si_MIXTURE_par_sum$Subjects$MIXabs)
si_MIXTURE_par_sum$Signature <- "Small intestine"
si_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(si_counts_par_sum) <- toupper(rownames(si_counts_par_sum))
rownames(si_sig) <- toupper(rownames(si_sig))
write.table(data.frame("X1"=rownames(si_sig),si_sig), file = "Figure_4/Data/MS/si_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

si_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = si_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/si_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
si_quanTIseq_par_sum <- si_quanTIseq_par_sum$Abs[,-ncol(si_quanTIseq_par_sum$Abs)]
si_quanTIseq_par_sum <- as.data.frame(t(si_quanTIseq_par_sum))
si_quanTIseq_par_sum$Signature <- "Small intestine"
si_quanTIseq_par_sum$Method <- "quanTIseq"
names(si_quanTIseq_par_sum) <- gsub("[.]"," ", names(si_quanTIseq_par_sum))

si_null_est <- rbind(si_CIBERSORT_par_sum,si_EPIC_par_sum,si_MIXTURE_par_sum,si_quanTIseq_par_sum)

s_sig <- read.csv("Figure_4/Data/MS/S.csv", row.names = 1, check.names=FALSE)

s_CIBERSORT_par_sum <- read_csv("Figure_4/Data/CIBERSORT_estimations/s_CIBERSORT_par_sum.csv")
s_CIBERSORT_par_sum <- as.data.frame(s_CIBERSORT_par_sum)
rownames(s_CIBERSORT_par_sum) <- s_CIBERSORT_par_sum[,1]
s_CIBERSORT_par_sum <- s_CIBERSORT_par_sum[-c(2:10), -c(1,10,11,12,13)]
s_CIBERSORT_par_sum$Signature <- "Spleen"
s_CIBERSORT_par_sum$Method <- "CIBERSORT"

s_EPIC_par_sum <- MyEPIC(bulk = s_counts_par_sum,
                         reference =  list(refProfiles=s_sig, sigGenes = rownames(s_sig)),
                         withOtherCells = FALSE)
s_EPIC_par_sum <- s_EPIC_par_sum$absCellFractions
s_EPIC_par_sum <- as.data.frame(s_EPIC_par_sum)
s_EPIC_par_sum$Signature <- "Spleen"
s_EPIC_par_sum$Method <- "EPIC"

s_MIXTURE_par_sum <- MIXTURE(expressionMatrix = s_counts_par_sum,      
                             signatureMatrix = s_sig,       
                             iter = 0L,
                             functionMixture = nu.svm.robust.RFE,
                             useCores = 1,
                             verbose = TRUE,
                             nullDist = "PopulationBased") 
s_MIXTURE_par_sum <- as.data.frame(s_MIXTURE_par_sum$Subjects$MIXabs)
s_MIXTURE_par_sum$Signature <- "Spleen"
s_MIXTURE_par_sum$Method <- "MIXTURE"

rownames(s_counts_par_sum) <- toupper(rownames(s_counts_par_sum))
rownames(s_sig) <- toupper(rownames(s_sig))
write.table(data.frame("X1"=rownames(s_sig),s_sig), file = "Figure_4/Data/MS/s_sig_signature.txt", sep = c("\t"), row.names = FALSE, quote = FALSE)

s_quanTIseq_par_sum <-deconvolute_quantiseq.default2(mix.mat = s_counts_par_sum,
                                                     signame = "Figure_4/Data/MS/s_sig",
                                                     arrays = FALSE,
                                                     tumor = FALSE,
                                                     mRNAscale = FALSE, method = "lsei", btotalcells = FALSE, rmgenes = "unassigned")
s_quanTIseq_par_sum <- s_quanTIseq_par_sum$Abs[,-ncol(s_quanTIseq_par_sum$Abs)]
s_quanTIseq_par_sum <- as.data.frame(t(s_quanTIseq_par_sum))
s_quanTIseq_par_sum$Signature <- "Spleen"
s_quanTIseq_par_sum$Method <- "quanTIseq"
names(s_quanTIseq_par_sum) <- gsub("[.]"," ", names(s_quanTIseq_par_sum))

s_null_est <- rbind(s_CIBERSORT_par_sum,s_EPIC_par_sum,s_MIXTURE_par_sum,s_quanTIseq_par_sum)

null_est <-rbind(k_null_est,l_null_est,mg_null_est,m_null_est,si_null_est,s_null_est)
rownames(null_est) <- NULL

#Figure

library(reshape2)
library(tidyverse)
library(ggpubr)

null_est_melt <- melt(null_est)

NCT_null_test <- null_est_melt %>% filter(value != 0) %>% group_by(Method, Signature) %>% summarize(count=n()) %>% ungroup() 

NCT_null_test$Method <- as.factor(NCT_null_test$Method)
NCT_null_test$variable <- as.factor(NCT_null_test$variable)

my_comparisons <- list( c("CIBERSORT", "quanTIseq"), c("EPIC", "quanTIseq"), c("MIXTURE", "quanTIseq") )
Figure_4 <- ggboxplot(NCT_null_test, x = "Method", y = "count", palette = "jco", xlab ="Method", ylab = "NCT")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 10) +
  scale_y_continuous(
    breaks = get_breaks(by = 1, from = 0, to = 8)
  )

Sup_Figure_1 <- ggpaired(subset(NCT_null_test, Method %in% c("MIXTURE", "CIBERSORT")), x = "Method", y = "count",
         line.color = "gray", line.size = 0.4,
         palette = "jco", xlab ="Method", ylab = "NCT")+
  stat_compare_means(paired = TRUE) 
