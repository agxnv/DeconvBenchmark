library(ADAPTS)
library(WGCNA)
library(DeconRNASeq)
library(reshape2)
library(tidyverse)
library(scales)
library(MIXTURE)
library(immunedeconv)
library(reshape2)

#Code for generating Simulated mixtures:
#coef <- as.integer((ncol(sigs$K)) * 0.7)
#Simulation <- SimulatedMixtures(SM, coef, 500, noisy = TRUE)

sigs <- list(BM = read.csv("Figure_2-3/Data/MS/BM.csv", row.names = 1, check.names=FALSE),
             K = read.csv("Figure_2-3/Data/MS/K.csv", row.names = 1, check.names=FALSE),
             L = read.csv("Figure_2-3/Data/MS/L.csv", row.names = 1, check.names=FALSE),
             LU = read.csv("Figure_2-3/Data/MS/LU.csv", row.names = 1, check.names=FALSE),
             MG = read.csv("Figure_2-3/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

load('Figure_2-3/Data/simulated_biastest.RData')

CIBERSORTx_mixsims_list <- list(BM = read_table2("Figure_2-3/Data/CIBERSORTx_estimations/BoneMarrow_CIBERSORTx.txt")[,c(1:11)],
                              K = read_table2("Figure_2-3/Data/CIBERSORTx_estimations/Kidney_CIBERSORTx.txt")[,c(1:9)],
                              L = read_table2("Figure_2-3/Data/CIBERSORTx_estimations/Liver_CIBERSORTx.txt")[,c(1:9)],
                              LU = read_table2("Figure_2-3/Data/CIBERSORTx_estimations/Lung_CIBERSORTx.txt")[,c(1:10)],
                              MG = read_table2("Figure_2-3/Data/CIBERSORTx_estimations/MammaryGland_CIBERSORTx.txt")[,c(1:9)])

CIBERSORTx_mixsims_list <- lapply(CIBERSORTx_mixsims_list, function(x) setNames(melt(x), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") 

CIBERSORTx_mixsims_list$variable <- gsub("\\."," ",CIBERSORTx_mixsims_list$variable) 

CIBERSORTx_BA <- left_join(CIBERSORTx_mixsims_list,
          betasims,
          by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "CIBERSORTxx")

DCQ_mixsims_list <- mapply(FUN = estCellPercent.DCQ,refExpr = sigs, geneExpr = mixsims,SIMPLIFY=FALSE)
DCQ_mixsims_list <- lapply(DCQ_mixsims_list, function(x) setNames(melt(t((x/100)[!(row.names(x) == "others"), ])), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature")

DCQ_BA <- left_join(DCQ_mixsims_list,
                    betasims,
                    by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "DCQ")

DeconRNASeq_mixsims_list <- mapply(FUN = estCellPercent.DeconRNASeq,refExpr = sigs, geneExpr = mixsims,SIMPLIFY=FALSE)
DeconRNASeq_mixsims_list <- lapply(DeconRNASeq_mixsims_list, function(x) setNames(melt(t(x/100)), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature")
DeconRNASeq_mixsims_list <- DeconRNASeq_mixsims_list[DeconRNASeq_mixsims_list$variable != "others",]
DeconRNASeq_BA <- left_join(DeconRNASeq_mixsims_list,
                            betasims,
                            by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "DeconRNASeq")

EPIC_mixsims_list <- mapply(FUN = EPIC::EPIC, bulk = mixsims, reference = lapply(sigs, function(x) list(refProfiles = x, sigGenes = rownames(x))), withOtherCells = FALSE, 
                    SIMPLIFY=FALSE)

EPIC_mixsims_list <- lapply(EPIC_mixsims_list, function(x) setNames(melt(x$cellFractions), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") 

EPIC_BA <- left_join(EPIC_mixsims_list,
                            betasims,
                            by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "EPIC")

MIXTURE_mixsims_list <- mapply(FUN = MIXTURE, expressionMatrix = mixsims,
          signatureMatrix = sigs, 
          SIMPLIFY=FALSE)

MIXTURE_mixsims_list <- lapply(MIXTURE_mixsims_list, function(x) setNames(melt(x$Subjects$MIXprop), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") 

MIXTURE_BA <- left_join(MIXTURE_mixsims_list,
                     betasims,
                     by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "MIXTURE")

mixsims <- lapply(mixsims,function(x) {rownames(x) <- toupper(rownames(x));x})

quanTIseq_mixsims_list <- mapply(FUN = deconvolute_quantiseq.default,
                         mix.mat = mixsims, 
                         arrays = FALSE, 
                         signame = paste0("Figure_2-3/Data/MS/",names(sigs)),
                         tumor = FALSE, 
                         mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                         rmgenes = "unassigned", 
                         SIMPLIFY=FALSE)

quanTIseq_mixsims_list <- lapply(quanTIseq_mixsims_list, function(x) {x <- x %>% select(-c("Other"));
                                                              x <- melt(as.data.frame(x));
                                                              names(x) <- c("ID","variable","betahat");
                                                              x$variable <- gsub("\\."," ", x$variable);
                                                              x}) %>% 
  bind_rows(.id = "Signature")

quanTIseq_BA <- left_join(quanTIseq_mixsims_list,
                        betasims,
                        by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "quanTIseq")



source('Figure_2-3/Data/BlandAltman.R')
library(plyr)

BM_BA <- BA.plot(rbind(CIBERSORTx_BA[CIBERSORTx_BA$Signature == "BM",c(2:5)],DCQ_BA[DCQ_BA$Signature == "BM",c(2:5)],DeconRNASeq_BA[DeconRNASeq_BA$Signature == "BM",c(2:5)],EPIC_BA[EPIC_BA$Signature == "BM",c(2:5)],MIXTURE_BA[MIXTURE_BA$Signature == "BM",c(2:5)],quanTIseq_BA[quanTIseq_BA$Signature == "BM",c(2:5)]),title = "Bone marrow",graph = "Bone marrow")
K_BA <- BA.plot(rbind(CIBERSORTx_BA[CIBERSORTx_BA$Signature == "K",c(2:5)],DCQ_BA[DCQ_BA$Signature == "K",c(2:5)],DeconRNASeq_BA[DeconRNASeq_BA$Signature == "K",c(2:5)],EPIC_BA[EPIC_BA$Signature == "K",c(2:5)],MIXTURE_BA[MIXTURE_BA$Signature == "K",c(2:5)],quanTIseq_BA[quanTIseq_BA$Signature == "K",c(2:5)]),title = "Kidney",graph = "Kidney")
L_BA <- BA.plot(rbind(CIBERSORTx_BA[CIBERSORTx_BA$Signature == "L",c(2:5)],DCQ_BA[DCQ_BA$Signature == "L",c(2:5)],DeconRNASeq_BA[DeconRNASeq_BA$Signature == "L",c(2:5)],EPIC_BA[EPIC_BA$Signature == "L",c(2:5)],MIXTURE_BA[MIXTURE_BA$Signature == "L",c(2:5)],quanTIseq_BA[quanTIseq_BA$Signature == "L",c(2:5)]),title = "Liver",graph = "Liver")
LU_BA <- BA.plot(rbind(CIBERSORTx_BA[CIBERSORTx_BA$Signature == "LU",c(2:5)],DCQ_BA[DCQ_BA$Signature == "LU",c(2:5)],DeconRNASeq_BA[DeconRNASeq_BA$Signature == "LU",c(2:5)],EPIC_BA[EPIC_BA$Signature == "LU",c(2:5)],MIXTURE_BA[MIXTURE_BA$Signature == "LU",c(2:5)],quanTIseq_BA[quanTIseq_BA$Signature == "LU",c(2:5)]),title = "Lung",graph = "Lung")
MG_BA <- BA.plot(rbind(CIBERSORTx_BA[CIBERSORTx_BA$Signature == "MG",c(2:5)],DCQ_BA[DCQ_BA$Signature == "MG",c(2:5)],DeconRNASeq_BA[DeconRNASeq_BA$Signature == "MG",c(2:5)],EPIC_BA[EPIC_BA$Signature == "MG",c(2:5)],MIXTURE_BA[MIXTURE_BA$Signature == "MG",c(2:5)],quanTIseq_BA[quanTIseq_BA$Signature == "MG",c(2:5)]),title = "Mammary gland",graph = "Mammary gland")

Sup_Table2 <- rbind(CIBERSORTx_BA,DeconRNASeq_BA,DCQ_BA,EPIC_BA,MIXTURE_BA,quanTIseq_BA) %>%
  mutate(result = case_when(betahat > 0 & betasim > 0 ~ "TP", 
                            betahat > 0 & betasim == 0 ~ "FP",
                            betahat == 0 & betasim == 0 ~ "TN",
                            betahat == 0 & betasim > 0 ~ "FN"
  )) %>% 
  group_by(result,Method) %>% 
  dplyr::summarize(count=n()) %>% 
  ungroup() %>% 
  complete(result, nesting(Method), fill = list(count = 0))

Sup_Table2 <- dcast(data = Sup_Table2,formula = Method~result,fun.aggregate = sum,value.var = "count")

Sup_Table2 <- Sup_Table2 %>% 
  transform(Se = TP / (TP + FN)) %>% 
  transform(Sp = TN / (TN + FP)) %>% 
  transform(PPV = TP / (TP + FP)) %>% 
  transform(NPV = TN / (TN + FN)) %>% 
  transform(F1 = (2 * TP) / ((2 * TP) + FP + FN)) %>%
  transform(DOP = sqrt(((Se - 1)^2)+((Sp - 1)^2)+((PPV - 1)^2)+((NPV - 1)^2))) %>% 
  transform(ER = (FP + FN)/(FP + FN + TP + TN))

#define quantiles of interest
q = c(.25, .5, .75)

Sup_Table3 <- rbind(CIBERSORTx_BA,DeconRNASeq_BA,DCQ_BA,EPIC_BA,MIXTURE_BA,quanTIseq_BA) %>%
  group_by(Method,Signature) %>%
  dplyr::summarize(quant25 = quantile(difs, probs = q[1]), 
                   quant50 = quantile(difs, probs = q[2]),
                   quant75 = quantile(difs, probs = q[3]),
                   min = min(difs),
                   max = max(difs),
                   mean = mean(difs),
                   sd = sd(difs),
                   range = max(difs) - min(difs),
                   correlation = cor(betasim,betahat)) %>% 
  ungroup()

#####

betasims_counts <- betasims[betasims$betasim > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(simcelltypes = n()) %>%
  ungroup() 

CIBERSORTx_counts_list <- CIBERSORTx_mixsims_list[CIBERSORTx_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "CIBERSORTx")

DCQ_counts_list <- DCQ_mixsims_list[DCQ_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "DCQ")

DeconRNASeq_counts_list <- DeconRNASeq_mixsims_list[DeconRNASeq_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "DeconRNASeq")

EPIC_counts_list <- EPIC_mixsims_list[EPIC_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "EPIC")

MIXTURE_counts_list <- MIXTURE_mixsims_list[MIXTURE_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "MIXTURE")

quanTIseq_counts_list <- quanTIseq_mixsims_list[quanTIseq_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "quanTIseq")

counts_a <- rbind(CIBERSORTx_counts_list,DeconRNASeq_counts_list,DCQ_counts_list,EPIC_counts_list,MIXTURE_counts_list,quanTIseq_counts_list)
counts_FINAL <- left_join(as.data.frame(counts_a),
                          as.data.frame(betasims_counts),
                          by = c("ID", "Signature"))[,c(2:5)]
#Figure

counts_total_heatmap <- counts_FINAL %>%  group_by(estcelltypes,simcelltypes, Method) %>% dplyr::summarize( total = n())
counts_totalbeta_heatmap <- counts_FINAL %>%  group_by(simcelltypes, Method) %>% dplyr::summarize( totalbeta = n())

conts_percent_heatmap <- counts_totalbeta_heatmap %>% left_join(counts_total_heatmap,
                                                                by = c("simcelltypes", "Method"))

conts_percent_heatmap <- transform(conts_percent_heatmap, percent = total / totalbeta)

conts_percent_heatmap$total <- NULL
conts_percent_heatmap$totalbeta <- NULL

conts_percent_heatmap <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method) %>%
  complete(estcelltypes = full_seq(1:10, 1), fill = list(percent = 0))


greater7 <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes > 7])) %>% 
  ungroup()
greater7$estcelltypes <- ">7"

lower2 <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes < 2])) %>% 
  ungroup()
lower2$estcelltypes <- "<2"

conts_percent_heatmap_final <- rbind(greater7,lower2,subset(conts_percent_heatmap, estcelltypes >= 2 & estcelltypes <= 7))

conts_percent_heatmap_final <- conts_percent_heatmap_final %>%
  group_by(simcelltypes, Method, estcelltypes) %>%
  dplyr::summarise(percent = sum(percent)/5) %>% 
  ungroup()

conts_percent_heatmap_final <- rbind(greater7,lower2,subset(conts_percent_heatmap, estcelltypes >= 2 & estcelltypes <= 7))

conts_percent_heatmap_final$estcelltypes <- conts_percent_heatmap_final$estcelltypes %>% 
  factor(levels=c("<2","2","3","4","5","6","7",">7"))

conts_percent_heatmap_final$simcelltypes <- conts_percent_heatmap_final$simcelltypes %>% 
  factor(levels=c("2","3","4","5","6","7"))

library(scales)
Figure_3 <- ggplot(conts_percent_heatmap_final, aes(x=simcelltypes, y=estcelltypes, fill=percent, text=percent)) + 
  geom_tile() + 
  geom_text(aes(label = scales::percent(percent,accuracy = 0.1L)), size = 8 / .pt)  +
  facet_grid(~Method) +
  xlab(bquote(K[r])) + ylab("Estimated NCTs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Percentage", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue",labels = percent) +
  theme(legend.position = "none")

counts_total_heatmap <- counts_FINAL %>%  group_by(estcelltypes,simcelltypes, Method, Signature) %>% dplyr::summarize( total = n())
counts_totalbeta_heatmap <- counts_FINAL %>%  group_by(simcelltypes, Method, Signature) %>% dplyr::summarize( totalbeta = n())

conts_percent_heatmap <- counts_totalbeta_heatmap %>% left_join(counts_total_heatmap,
                                                                by = c("simcelltypes", "Method","Signature"))

conts_percent_heatmap <- transform(conts_percent_heatmap, percent = total / totalbeta)

conts_percent_heatmap$total <- NULL
conts_percent_heatmap$totalbeta <- NULL

conts_percent_heatmap <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method, Signature) %>%
  complete(estcelltypes = full_seq(1:10, 1), fill = list(percent = 0))

conts_percent_heatmap$percent <- round(conts_percent_heatmap$percent,3) 

greater7 <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method, Signature) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes > 7])) %>% 
  ungroup()
greater7$estcelltypes <- ">7"

lower2 <- conts_percent_heatmap %>%
  group_by(simcelltypes, Method, Signature) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes < 2])) %>% 
  ungroup()
lower2$estcelltypes <- "<2"

conts_percent_heatmap_final <- rbind(greater7,lower2,subset(conts_percent_heatmap, estcelltypes >= 2 & estcelltypes <= 7))

conts_percent_heatmap_final <- conts_percent_heatmap_final %>%
  group_by(simcelltypes, Method, estcelltypes) %>%
  dplyr::summarise(percent = sum(percent)/5) %>% 
  ungroup()

conts_percent_heatmap_final <- rbind(greater7,lower2,subset(conts_percent_heatmap, estcelltypes >= 2 & estcelltypes <= 7))

conts_percent_heatmap_final$estcelltypes <- conts_percent_heatmap_final$estcelltypes %>% 
  factor(levels=c("<2","2","3","4","5","6","7",">7"))

conts_percent_heatmap_final$simcelltypes <- conts_percent_heatmap_final$simcelltypes %>% 
  factor(levels=c("2","3","4","5","6","7"))

SupFigure_1 <- ggplot(conts_percent_heatmap_final, aes(x=simcelltypes, y=estcelltypes, fill=percent, text=percent)) + 
  geom_tile() + 
  geom_text(aes(label = scales::percent(percent,accuracy = 0.1L)), size = 8 / .pt)  +
  facet_grid(Signature~Method) +
  xlab(bquote(K[r])) + ylab("Estimated NCTs") + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Percentage", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue",labels = percent) +
  theme(legend.position = "none")
