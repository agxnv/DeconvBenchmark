library(ADAPTS)
library(WGCNA)
library(DeconRNASeq)
library(reshape2)
library(tidyverse)
library(scales)
library(MIXTURE)
library(immunedeconv)
library(reshape2)
library(ggpubr)
library(readr)


sigs <- list(BM = read.csv("Figure_5/Data/MS/BM.csv", row.names = 1, check.names=FALSE),
             K = read.csv("Figure_5/Data/MS/K.csv", row.names = 1, check.names=FALSE),
             L = read.csv("Figure_5/Data/MS/L.csv", row.names = 1, check.names=FALSE),
             LU = read.csv("Figure_5/Data/MS/LU.csv", row.names = 1, check.names=FALSE),
             MG = read.csv("Figure_5/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

load("Figure_5/Data/simulated_immune.RData")


inm_test <- list( BM = bm_counts_inm_tpm,
                  K = k_counts_inm_tpm,
                  L = l_counts_inm_tpm,
                  LU = lu_counts_inm_tpm,
                  MG = mg_counts_inm_tpm)

CIBERSORTx_inmtest_list <- list(BM = read_csv("Figure_5/Data/CIBERSORTx_estimations/bm_inm_CIBERSORTx.csv")[,c(1:11)],
  K = read_csv("Figure_5/Data/CIBERSORTx_estimations/k_inm_CIBERSORTx.csv")[,c(1:9)],
                                L = read_csv("Figure_5/Data/CIBERSORTx_estimations/l_inm_CIBERSORTx.csv")[,c(1:9)],
                                LU = read_csv("Figure_5/Data/CIBERSORTx_estimations/lu_inm_CIBERSORTx.csv")[,c(1:10)],
                                MG = read_csv("Figure_5/Data/CIBERSORTx_estimations/mg_inm_CIBERSORTx.csv")[,c(1:9)])

CIBERSORTx_inmtest_list <- lapply(CIBERSORTx_inmtest_list, function(x) setNames(melt(x, id.vars = "Mixture"), nm = c("ID","variable","beta_hat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "CIBERSORTx")
CIBERSORTx_inmtest_list$variable <- gsub("\\."," ", CIBERSORTx_inmtest_list$variable)
CIBERSORTx_inmtest_list$ID <- gsub("\\."," ", CIBERSORTx_inmtest_list$ID) 

DCQ_inmtest_list <- mapply(FUN = estCellPercent.DCQ ,refExpr = sigs, geneExpr = inm_test,SIMPLIFY=FALSE)
DCQ_inmtest_list <- lapply(DCQ_inmtest_list, function(x) setNames(melt(t((x/100))), nm = c("ID","variable","beta_hat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "DCQ")
DCQ_inmtest_list <- DCQ_inmtest_list[DCQ_inmtest_list$variable != "others",]

DeconRNASeq_inmtest_list <- mapply(FUN = estCellPercent.DeconRNASeq,refExpr = sigs, geneExpr = inm_test,SIMPLIFY=FALSE)
DeconRNASeq_inmtest_list <- lapply(DeconRNASeq_inmtest_list, function(x) setNames(melt(t(x/100)), nm = c("ID","variable","beta_hat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "DeconRNASeq")
DeconRNASeq_inmtest_list <- DeconRNASeq_inmtest_list[DeconRNASeq_inmtest_list$variable != "others",]

EPIC_inmtest_list <- mapply(FUN = EPIC::EPIC, bulk = inm_test, reference = lapply(sigs, function(x) list(refProfiles = x, sigGenes = rownames(x))), withOtherCells = FALSE, 
                             SIMPLIFY=FALSE)

EPIC_inmtest_list <- lapply(EPIC_inmtest_list, function(x) setNames(melt(x$cellFractions), nm = c("ID","variable","beta_hat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "EPIC")

MIXTURE_inmtest_list <- mapply(FUN = MIXTURE, expressionMatrix = inm_test,
                                signatureMatrix = sigs, 
                                SIMPLIFY=FALSE)

MIXTURE_inmtest_list <- lapply(MIXTURE_inmtest_list, function(x) setNames(melt(x$Subjects$MIXprop), nm = c("ID","variable","beta_hat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "MIXTURE")

inm_test <- lapply(inm_test,function(x) {rownames(x) <- toupper(rownames(x));x})

quanTIseq_inmtest_list <- mapply(FUN = deconvolute_quantiseq.default,
                                  mix.mat = inm_test, 
                                  arrays = FALSE, 
                                  signame = paste0("Figure_1/Data/MS/",names(sigs)),
                                  tumor = FALSE, 
                                  mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                                  rmgenes = "unassigned", 
                                  SIMPLIFY=FALSE)

quanTIseq_inmtest_list <- lapply(quanTIseq_inmtest_list, function(x) {x <- as.data.frame(x) %>% select(-c("Other"));
x <- melt(x, id.vars="Sample");
names(x) <- c("ID","variable","beta_hat");
x$variable <- gsub("\\."," ", x$variable);
x}) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "quanTIseq")

inm_est <-rbind(CIBERSORTx_inmtest_list,DCQ_inmtest_list,
                 DeconRNASeq_inmtest_list, EPIC_inmtest_list,
                 MIXTURE_inmtest_list,quanTIseq_inmtest_list)

#Table 1

inm_est_melt <- inm_est

Table1 <- inm_est_melt %>% mutate(result = case_when(ID == variable & beta_hat > 0 ~ "TP", 
                                                     ID != variable & beta_hat > 0 ~ "FP",
                                                     ID == variable & beta_hat == 0 ~ "FN",
                                                     ID != variable & beta_hat == 0 ~ "TN")) %>% 
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
Table1 <- transform(Table1, ER = (FP + FN)/(FP + FN + TP + TN))

#Figures

inm_est_melt <- inm_est


inm_est_TP <- inm_est_melt %>% filter(ID == variable, beta_hat > 0)
inm_est_TP$ID <- as.factor(inm_est_TP$ID)

inm_est_TP$CT <- paste0(inm_est_TP$Signature,"-",inm_est_TP$ID)

n <- inm_est_TP %>% 
  group_by(Method) %>% 
  dplyr::summarize(count = n())
n$n <- paste0("n = ", n$count)

Figure_5a <-  ggplot(inm_est_TP, aes(factor(Method), beta_hat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for TP") + scale_y_continuous(labels = scales::percent)

Sup_Figure_3 <- ggplot(inm_est_TP, aes(factor(Signature), beta_hat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0, alpha = 0.5, height = 0)+
  labs(x = "Signature",
       y = "Predicted percentages for TP") +
  facet_wrap(~Method) + scale_y_continuous(labels = scales::percent)

inm_est_FP <- inm_est_melt %>% filter(ID != variable, beta_hat > 0)

n1 <- inm_est_FP %>% 
  group_by(Method) %>% 
  dplyr::summarize(count = n())
n1$n <- paste0("n = ", n1$count)

Figure_5b <-  ggplot(inm_est_FP, aes(factor(Method), beta_hat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n1, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for FP") + scale_y_continuous(labels = scales::percent)

inm_est_FP$Signature <- as.factor(inm_est_FP$Signature)

Sup_Figure_4 <- ggplot(inm_est_FP, aes(x = Signature, y =beta_hat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  labs(x = "Signature",
       y = "Predicted percentage for FP") +
  facet_wrap(~ Method) + scale_y_continuous(labels = scales::percent)

TP_CM <- subset(inm_est_TP, Method %in% c("CIBERSORTx","MIXTURE"))
n2 <- TP_CM %>% 
  group_by(Method) %>% 
  dplyr::summarize(count = n())
n2$n <- paste0("n = ", n2$count)

inm_est_melt <- inm_est
inm_est_melt$beta_hat <- as.double(inm_est_melt$beta_hat)

inm_est_TP <- inm_est_melt %>% filter(ID == variable)
inm_est_TP$beta_hat <- as.numeric(inm_est_TP$beta_hat)
TP_CM <- subset(inm_est_TP, Method %in% c("CIBERSORTx","MIXTURE"))

Sup_Figure_2 <- ggpaired(TP_CM, x = "Method", y = "beta_hat",
                         line.color = "gray", line.size = 0.4,
                         palette = "jco", xlab ="Method", ylab = "Predicted percentage for TP")+
  stat_compare_means(paired = TRUE) + scale_y_continuous(labels = scales::percent)

Figure_5 <- ggarrange(Figure_5a, Figure_5b,
                      labels = c("A", "B"),
                      ncol = 1, nrow = 2)

