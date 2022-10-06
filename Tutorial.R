#This document its a step-by-step guide to perform the benchmarking algorithm 
#for deconvolution methods and molecular signatures described in the paper:
#"Agreement and cell-type identification metrics for an appropriate evaluation 
#of molecular signature-based deconvolution methods". 
#These tests allow to evaluate the estimation accuracy and classificatory 
#competence of each deconvolution method, detect bias over its results, and 
#serve as a fair comparision framework between other molecular signature based 
#deconvolution methods. They're also useful to evaluate how the performance of 
#a DM changes when paired to different molecular signatures. 

#In the examples presented bellow, we're going to test the EPIC and MIXTURE
#methods, paired to two murine tissue-specific MSs (Chen et al.,
#DOI: 10.1093/bioinformatics/btz672).  

#Step 0: Load the DMs packages and the molecular signatures

library(ADAPTS)
library(WGCNA)
library(DeconRNASeq)
library(reshape2)
library(tidyverse)
library(scales)
library(MIXTURE)
library(immunedeconv)
library(readr)
library(ComplexHeatmap)
library(circlize)

sigs <- list(K = read.csv("Data/MS/K.csv", row.names = 1, check.names=FALSE),
  MG = read.csv("Data/MS/MG.csv", row.names = 1, check.names=FALSE))

# Test 1: Self-test
#For this test, each tested method's input must be the MS that is paired to it, 
#so every cell type GEP in the tested MS is going to act as a separated case.
#The desired output it's that each tested DM detect only 100% of every cell
#type in the MS. 


EPIC_selftest_list <- mapply(FUN = EPIC, 
                             bulk = sigs, 
                             reference = lapply(sigs, function(x)
                             list(refProfiles = x, sigGenes = rownames(x))),
                             withOtherCells = FALSE, 
                             SIMPLIFY=FALSE)

MIXTURE_selftest_list <- mapply(FUN = MIXTURE, 
                                expressionMatrix = sigs,
                                signatureMatrix = sigs, 
                                SIMPLIFY=FALSE)

#After running the deconvolutions, a heatmap with the results can be 
#generated to compare them easily.

EPIC_selftest_list <- lapply(EPIC_selftest_list, function(x) 
  as.data.frame(x$cellFractions) %>% 
  mutate(X1 = rownames(x)) %>% 
  arrange(X1) %>%
  na_if(0)) 

EPIC_selftest_list <- lapply(EPIC_selftest_list, function(x) 
  data.matrix(x[,!(colnames(x) == "X1")]
              [,order(colnames(x[,!(colnames(x) == "X1")]))]))

colnames(EPIC_selftest_list$MG) <- 
  rownames(EPIC_selftest_list$MG) <- 
  c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")

colnames(EPIC_selftest_list$K) <- 
  rownames(EPIC_selftest_list$K) <-
  c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")


MIXTURE_selftest_list <- lapply(MIXTURE_selftest_list, function(x) 
  as.data.frame(x$Subjects$MIXprop) %>% 
                      mutate(X1 = rownames(x)) %>% 
                      arrange(X1) %>% na_if(0)) 

MIXTURE_selftest_list <- lapply(MIXTURE_selftest_list, function(x)
  data.matrix(x[,!(colnames(x) == "X1")][,order(
    colnames(x[,!(colnames(x) == "X1")]))]))

colnames(MIXTURE_selftest_list$MG) <-
  rownames(MIXTURE_selftest_list$MG) <-
  cell.types.names.MG

colnames(MIXTURE_selftest_list$K) <-
  rownames(MIXTURE_selftest_list$K) <-
  cell.types.names.k.l.mg


library(ComplexHeatmap)
library(circlize)

MG_heatmap <- 
  Heatmap(EPIC_selftest_list$MG, cluster_rows = FALSE,
          column_title_gp = gpar(fontsize=15), show_column_names = TRUE,
          show_row_names = FALSE, cluster_columns = FALSE,
          column_title = "EPIC",name = "Mammary gland",
          col = colorRamp2(c(0, 1), c("blue",  "red")),
          show_heatmap_legend = F) +
  Heatmap(MIXTURE_selftest_list$MG, cluster_rows = FALSE,
          column_title_gp = gpar(fontsize=15), show_row_names = FALSE,
          cluster_columns = FALSE, column_title = "MIXTURE", name = "MIXTURE",
          col = colorRamp2(c(0, 1), c("blue",  "red")),
          show_heatmap_legend = FALSE, show_column_names = TRUE) 
  
#To evaluate how the correlation between the profiles of the tested MS's affect
#the falsely detected cell-type coefficients, we can plot them against the 
#inter-celltype profile correlation.

sigs_cor_5 <- list(K = read.csv("Data/MS/K.csv", row.names = 1, check.names=FALSE),
                   L = read.csv("Data/MS/L.csv", row.names = 1, check.names=FALSE),
                   MG = read.csv("Data/MS/MG.csv", row.names = 1, check.names=FALSE))

sigs_cor_5 <- lapply(sigs_cor_5, function(x) {x <- x[,order(colnames(x))];
colnames(x) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T");
x <- as.matrix(round(cor(x),2));
x})

sigs_cor_5$K[lower.tri(sigs_cor_5$K, diag = TRUE)] <- NA
K_melted_table <- melt(sigs_cor_5$K, na.rm = TRUE, value.name = "cor")
K_melted_table$Signature <- "K"

sigs_cor_5$MG[lower.tri(sigs_cor_5$MG, diag = TRUE)] <- NA
MG_melted_table <- melt(sigs_cor_5$MG, na.rm = TRUE, value.name = "cor")
MG_melted_table$Signature <- "MG"

cortotal <- rbind(K_melted_table,MG_melted_table) %>% arrange(Signature)
cortotal[166:330,] <- cortotal[,c(2,1,3,4)]


EPIC_selftest_list_supfig1 <- lapply(EPIC_selftest_list, function(x) melt(x))
EPIC_selftest_list_supfig1 <-  EPIC_selftest_list_supfig1 %>% bind_rows(.id = "Signature")
EPIC_selftest_list_supfig1 <- EPIC_selftest_list_supfig1[EPIC_selftest_list_supfig1$Var1 != EPIC_selftest_list_supfig1$Var2,]
EPIC_selftest_list_supfig1 <- EPIC_selftest_list_supfig1[complete.cases(EPIC_selftest_list_supfig1),]
EPIC_selftest_list_supfig1$Method <- "EPIC"

Supfig1_df <- rbind(EPIC_selftest_list_supfig1)
Supfig1_df<- left_join(Supfig1_df,
                       cortotal,
                       by = c("Signature","Var1","Var2"))

Cor_offdiagonal_analysis <- ggplot(Supfig1_df, aes(x=cor, y=value)) +
  geom_point() +
  facet_grid(~Method) +
  xlab("Correlation between the deconvolved MS profile and the detected cell type") +
  ylab("Estimated coefficient")+ stat_smooth(method="lm", se=FALSE, color="red")


#Test 2: Bias test
#In this test, simulated gene expression profiles are generated with the
#function 'MIXTURE::SimulatedMixtures'. 
#Step 1: Generate simulated samples

Simulation_K <- SimulatedMixtures(sigs$K, 5, 500, noisy = TRUE)
Simulation_MG <- SimulatedMixtures(sigs$MG, 7, 500, noisy = TRUE)

#In this example, 500 samples are generated that have a 
#maximum number of simulated cell-types of 5 and 7, with a noise component. 
#It returns a list with 2 components, a dataframe containing the simulated 
#bulk cell mixtures and another data frame with their corresponding
#beta coefficients.

mixsims <- list(MG = Simulation_MG$M,
                K = Simulation_K$M)

betasims <- list(MG = Simulation_MG$B,
                 K = Simulation_K$B) %>% 
  lapply(function(x) {x <- as.data.frame(x);
  x$ID <- paste0("Sim",seq(1:nrow(x)));
x <- melt(as.data.frame(x), id.vars = "ID", value.name = "betasim");
x}) %>% 
  bind_rows(.id = "Signature")

#Step 2: Deconvolve the simulated samples

EPIC_mixsims_list <- mapply(FUN = EPIC::EPIC, 
                            bulk = mixsims, 
                            reference = lapply(sigs, function(x) 
                              list(refProfiles = x, sigGenes = rownames(x))), 
                            withOtherCells = FALSE, 
                            SIMPLIFY=FALSE)

MIXTURE_mixsims_list <- mapply(FUN = MIXTURE, 
                               expressionMatrix = mixsims,
                               signatureMatrix = sigs, 
                               SIMPLIFY=FALSE)

#Step 3: Bland-Altmann method
#-3a: Calculate the difference between the estimated and simulated coefficients

EPIC_mixsims_list <- lapply(EPIC_mixsims_list, function(x) 
  setNames(melt(x$cellFractions), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") 

EPIC_BA <- left_join(EPIC_mixsims_list,
                     betasims,
                     by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "EPIC")


MIXTURE_mixsims_list <- lapply(MIXTURE_mixsims_list, function(x) 
  setNames(melt(x$Subjects$MIXprop), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") 

MIXTURE_BA <- left_join(MIXTURE_mixsims_list,
                        betasims,
                        by = c("Signature", "ID","variable"))[,c(1,4,5)] %>% 
  mutate(difs = betahat - betasim) %>% 
  mutate(Method = "MIXTURE")

#-3b: Bland-Altmann plot

source('Data/BlandAltman.R')

MG_BA <- BA.plot(rbind(EPIC_BA[EPIC_BA$Signature == "MG",c(2:5)],
                       MIXTURE_BA[MIXTURE_BA$Signature == "MG",c(2:5)]),
                 title = "Bone marrow",
                 graph = "Bone marrow")

K_BA <- BA.plot(rbind(EPIC_BA[EPIC_BA$Signature == "K",c(2:5)],
                      MIXTURE_BA[MIXTURE_BA$Signature == "K",c(2:5)]),
                title = "Kidney",
                graph = "Kidney")

#-3c: Summary statistics

q = c(.25, .5, .75)

BiasTest_table1 <- rbind(EPIC_BA,MIXTURE_BA) %>%
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

#Step 4: Sensitivity, specificity, F1-score, DOP and total error estimation

BiasTest_table2 <- rbind(EPIC_BA,MIXTURE_BA) %>%
  mutate(result = case_when(betahat > 0 & betasim > 0 ~ "TP", 
                            betahat > 0 & betasim == 0 ~ "FP",
                            betahat == 0 & betasim == 0 ~ "TN",
                            betahat == 0 & betasim > 0 ~ "FN")) %>% 
  group_by(result,Method) %>% 
  dplyr::summarize(count=n()) %>% 
  ungroup() %>% 
  complete(result, nesting(Method), fill = list(count = 0))

BiasTest_table2 <- dcast(data = BiasTest_table2,
                         formula = Method~result,
                         fun.aggregate = sum,
                         value.var = "count") %>% 
  transform(Se = TP / (TP + FN)) %>% 
  transform(Sp = TN / (TN + FP)) %>% 
  transform(PPV = TP / (TP + FP)) %>% 
  transform(NPV = TN / (TN + FN)) %>% 
  transform(F1 = (2 * TP) / ((2 * TP) + FP + FN)) %>%
  transform(DOP = 
  sqrt(((Se - 1)^2)+((Sp - 1)^2)+((PPV - 1)^2)+((NPV - 1)^2))) %>% 
  transform(ER = (FP + FN)/(FP + FN + TP + TN))
 
#Step 5: Estimate the number of cell types (NCT) detected by the tested DMs, 
#according to which MS has been paired to them and the simulated NCTs.

betasims_counts <- betasims[betasims$betasim > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(simcelltypes = n()) %>%
  ungroup()

EPIC_counts_list <- EPIC_mixsims_list[EPIC_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "EPIC")

MIXTURE_counts_list <- 
  MIXTURE_mixsims_list[MIXTURE_mixsims_list$betahat > 0,] %>% 
  group_by(ID,Signature) %>% 
  dplyr::summarise(estcelltypes = n()) %>%
  ungroup() %>% 
  mutate(Method = "MIXTURE")

counts_FINAL <- left_join(rbind(EPIC_counts_list,MIXTURE_counts_list),
                          betasims_counts,
                          by = c("ID", "Signature"))[,c(2:5)]

counts_total_heatmap <- counts_FINAL %>%
  group_by(estcelltypes,simcelltypes, Method, Signature) %>%
  dplyr::summarize( total = n())

counts_totalbeta_heatmap <- counts_FINAL %>%
  group_by(simcelltypes, Method, Signature) %>% 
  dplyr::summarize( totalbeta = n())

counts_percent_heatmap <- counts_totalbeta_heatmap %>% 
  left_join(counts_total_heatmap,
            by = c("simcelltypes", "Method","Signature")) %>% 
  transform(percent = total / totalbeta) %>% 
  transform(percent = round(percent,2))%>%
  select(-matches(c("total","totalbeta"))) %>% 
  group_by(simcelltypes, Method, Signature) %>%
  complete(estcelltypes = full_seq(1:10, 1), fill = list(percent = 0)) %>% 
  ungroup()

greater7 <- counts_percent_heatmap %>%
  group_by(simcelltypes, Method, Signature) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes > 7])) %>% 
  ungroup() %>% 
  mutate(estcelltypes = ">7")

lower2 <- counts_percent_heatmap %>%
  group_by(simcelltypes, Method, Signature) %>%
  dplyr::summarise(percent = sum(percent[estcelltypes < 2])) %>% 
  ungroup() %>% 
  mutate(estcelltypes = "<2")

counts_heatmap_by_MS <- rbind(greater7,
                                     lower2,
                                     subset(counts_percent_heatmap,
                                            estcelltypes >= 2 &
                                              estcelltypes <= 7))

counts_heatmap_by_MS$estcelltypes <- counts_heatmap_by_MS$estcelltypes %>% 
  factor(levels=c("<2","2","3","4","5","6","7",">7"))

counts_heatmap_by_MS$simcelltypes <- counts_heatmap_by_MS$simcelltypes %>% 
  factor(levels=c("2","3","4","5","6","7"))

Figure_NCT_heatmap <- 
  ggplot(counts_heatmap_by_MS, 
                      aes(x=simcelltypes,
                          y=estcelltypes,
                          fill=percent,
                          text=percent)) + 
  geom_tile() + 
  geom_text(aes(label = scales::percent(percent,accuracy = 0.1L)),
            size = 8 / .pt)  +
  facet_wrap(Signature~Method,scales="free",ncol = 2) +
  xlab(bquote(K[r])) + 
  ylab("Estimated NCTs") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Percentage", reverse=T)) + 
  scale_x_discrete(label=abbreviate) + 
  scale_fill_continuous(high = "red", low = "light blue",labels = percent) +
  theme(legend.position = "none")

#Test 3: NULL test
#For this test the GEP of a mixture of non-immune-cells is required. In this 
#example we're going to use semi-simulated pooled samples that were generated
#adding the GEP of parenchymatous and stromal cells from murine tissues (cita MCA). 
#It is also required to take the non-constrained-to-one coefficients outputed
#by the tested methods. Since EPIC doesn't present the "Absolute mode" feature,
#this test can't be performed with this method. 

load("Data/simulated_nonimmune.RData")

#Step 1: Load the simulated samples

null_test <- list(K = k_nonimmune,
                 MG = mg_nonimmune)

#Step 2: Run the deconvolutions

MIXTURE_nulltest_list <- mapply(FUN = MIXTURE, expressionMatrix = null_test,
                                signatureMatrix = sigs, 
                                SIMPLIFY=FALSE)

#Step 3: Count the NCTs detected by each method

MIXTURE_nulltest_list <- lapply(MIXTURE_nulltest_list, function(x)
  setNames(melt(x$Subjects$MIXabs), nm = c("ID","variable","beta_abs"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "MIXTURE")
MIXTURE_nulltest_list[is.na(MIXTURE_nulltest_list)] <- 0

NCT_null_test <- MIXTURE_nulltest_list %>% 
  filter(beta_abs > 0) %>% 
  group_by(Method, Signature,ID) %>% 
  dplyr::summarize(count=n()) %>%
  ungroup() %>%
  complete(Method, nesting(ID, Signature), fill = list(count = 0)) %>% 
  mutate(Method = as.factor(Method))

#The results can be plotted in boxplots. 

Nulltest_boxplot <- ggboxplot(NCT_null_test, x = "Method", 
                       y = "count", 
                       palette = "jco",
                       xlab ="Method", 
                       ylab = "NCT")+
  scale_y_continuous(
    breaks = get_breaks(by = 1, from = 0, to = 9)
  )

#Test 4: Unique cell test

#For this test the GEP of a pooled sample of just one type of immune-cells is 
#required. In this  example we're going to use semi-simulated pooled samples
#that were generated adding the GEP of murine kidney and mammary gland cells (cita MCA). 

#Step 1: Load the samples

load("Figure_5/Data/simulated_immune.RData")

inm_test <- list(K = k_counts_inm_tpm,
                 MG = mg_counts_inm_tpm)

#Step 2: Run the deconvolution

EPIC_inmtest_list <- mapply(FUN = EPIC::EPIC, 
                            bulk = inm_test, 
                            reference = lapply(sigs, function(x) 
                              list(refProfiles = x, 
                                   sigGenes = rownames(x))), 
                            withOtherCells = FALSE, 
                            SIMPLIFY=FALSE)

MIXTURE_inmtest_list <- mapply(FUN = MIXTURE,
                               expressionMatrix = inm_test,
                               signatureMatrix = sigs, 
                               SIMPLIFY=FALSE)

#Step 3: Sensitivity, specificity, F1-score, DOP and total error estimation

EPIC_inmtest_list <- lapply(EPIC_inmtest_list, function(x)
  setNames(melt(x$cellFractions), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "EPIC")

MIXTURE_inmtest_list <- lapply(MIXTURE_inmtest_list, function(x)
  setNames(melt(x$Subjects$MIXprop), nm = c("ID","variable","betahat"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "MIXTURE")

inm_est_melt <-rbind(EPIC_inmtest_list,
                MIXTURE_inmtest_list)

Uniquetest_Table <- inm_est_melt %>% 
  mutate(ID = as.character(ID))%>%
  mutate(variable = as.character(variable))%>%
  mutate(
  result = case_when(ID == variable & betahat > 0 ~ "TP", 
                     ID != variable & betahat > 0 ~ "FP",
                     ID == variable & betahat == 0 ~ "FN",
                     ID != variable & betahat == 0 ~ "TN")) %>% 
  group_by(result,Method) %>% 
  dplyr::summarize(count=n()) %>% 
  ungroup() %>% 
  complete(result, nesting(Method), fill = list(count = 0)) %>% 
  dcast(formula = Method~result,fun.aggregate = sum,value.var = "count") %>% 
  transform(Se = TP / (TP + FN)) %>% 
  transform(Sp = TN / (TN + FP)) %>% 
  transform(PPV = TP / (TP + FP)) %>% 
  transform(NPV = TN / (TN + FN)) %>% 
  transform(F1 = (2 * TP) / ((2 * TP) + FP + FN)) %>% 
  transform(DOP = 
              sqrt(((Se - 1)^2)+((Sp - 1)^2)+((PPV - 1)^2)+((NPV - 1)^2))) %>% 
  transform(ER = (FP + FN)/(FP + FN + TP + TN))

#True and false positive estimated coefficients can be compared to compare the
#accuracy of each tested method.

inm_est_TP <- inm_est_melt %>% 
  mutate(ID = as.character(ID))%>%
  mutate(variable = as.character(variable))%>%
  filter(ID == variable, betahat > 0)  %>%
  mutate(CT = paste0(Signature,"-",ID))

n <- inm_est_TP %>% 
  group_by(Method) %>% 
  dplyr::summarize(count = n()) %>% 
  mutate(n = paste0("n = ", count))


Uniquetest_TruePositives_Boxplot <-  ggplot(inm_est_TP, 
                                            aes(factor(Method), 
                                                betahat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for TP") +
  scale_y_continuous(labels = scales::percent)

Uniquetest_TruePositives_bySig_Boxplot <- ggplot(inm_est_TP,
                                                 aes(factor(Signature),
                                                     betahat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0, alpha = 0.5, height = 0)+
  labs(x = "Signature",
       y = "Predicted percentages for TP") +
  facet_wrap(~Method) + scale_y_continuous(labels = scales::percent)

inm_est_FP <- inm_est_melt %>% filter(ID != variable, betahat > 0)

n1 <- inm_est_FP %>% 
  group_by(Method) %>% 
  dplyr::summarize(count = n())
n1$n <- paste0("n = ", n1$count)

Uniquetest_FalsePositives_Boxplot <-  ggplot(inm_est_FP,
                                             aes(factor(Method),
                                                 betahat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.2, alpha = 0.5, height = 0) +
  geom_text(data = n1, aes(y = 1.05, label = n)) +
  labs(x = "Methods",
       y = "Predicted percentages for FP") + 
  scale_y_continuous(labels = scales::percent)

Uniquetest_FalsePositives_bySig_Boxplot <- ggplot(inm_est_FP,
                                                  aes(x = as.factor(Signature),
                                                      y =betahat)) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  labs(x = "Signature",
       y = "Predicted percentage for FP") +
  facet_wrap(~ Method) + scale_y_continuous(labels = scales::percent)

#To evaluate how the correlation between the profiles of the tested MS's affect
#the falsely detected cell-type coefficients, we can plot them against the 
#inter-cell-type profile correlation.

sigs_cor_5 <- list(K = read.csv("Data/MS/K.csv", row.names = 1,
                                check.names=FALSE),
                   MG = read.csv("Data/MS/MG.csv", row.names = 1,
                                 check.names=FALSE))

sigs_cor_5 <- lapply(sigs_cor_5, function(x) {x <- x[,order(colnames(x))];
x <- as.matrix(round(cor(x),2));
x})

sigs_cor_5$K[lower.tri(sigs_cor_5$K, diag = TRUE)] <- NA
K_melted_table <- melt(sigs_cor_5$K, na.rm = TRUE, value.name = "cor")
K_melted_table$Signature <- "K"
sigs_cor_5$MG[lower.tri(sigs_cor_5$MG, diag = TRUE)] <- NA
MG_melted_table <- melt(sigs_cor_5$MG, na.rm = TRUE, value.name = "cor")
MG_melted_table$Signature <- "MG"

cortotal <- rbind(K_melted_table,MG_melted_table) %>% arrange(Signature)
cortotal$names <- paste0(cortotal$Var1,", ",cortotal$Var2)
cortotal$names <- sapply(lapply(strsplit(cortotal$names, split = ",\\s*"),
                                sort), paste, collapse = ", ")
cortotal <- cortotal[,-c(1,2)]

inm_est_FP$names <- paste0(inm_est_FP$ID,", ",inm_est_FP$variable)
inm_est_FP$names <- sapply(lapply(strsplit(inm_est_FP$names, split = ",\\s*"),
                                  sort), paste, collapse = ", ")

inm_est_FP <- left_join(inm_est_FP,
                        cortotal,
                        by = c("names","Signature"))

cor_FP_analysis <- ggplot(inm_est_FP, aes(x=cor, y=beta_hat)) +
  geom_point() +
  facet_grid(~Method) +
  xlab("Correlation between the deconvolved cell type and the detected cell
       type") +
  ylab("Estimated coefficient")+ 
  stat_smooth(method="lm", se=FALSE, color="red")

