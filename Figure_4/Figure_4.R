library(ADAPTS)
library(WGCNA)
library(DeconRNASeq)
library(tidyverse)
library(scales)
library(MIXTURE)
library(immunedeconv)
library(reshape2)
library(ggpubr)
library(readr)


sigs <- list(K = read.csv("Figure_4/Data/MS/K.csv", row.names = 1, check.names=FALSE),
             L = read.csv("Figure_4/Data/MS/L.csv", row.names = 1, check.names=FALSE),
             LU = read.csv("Figure_4/Data/MS/LU.csv", row.names = 1, check.names=FALSE),
             MG = read.csv("Figure_4/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

load("Figure_4/Data/simulated_nonimmune.RData")

null_test <- list(K = k_nonimmune,
                  L = l_nonimmune,
                  LU = lu_nonimmune,
                  MG = mg_nonimmune)

CIBERSORTx_nulltest_list <- list(K = read_csv("Figure_4/Data/CIBERSORTx_estimations/k_par_CIBERSORTx.csv")[,c(1:9)],
                          L = read_csv("Figure_4/Data/CIBERSORTx_estimations/l_par_CIBERSORTx.csv")[,c(1:9)],
                          LU = read_csv("Figure_4/Data/CIBERSORTx_estimations/lu_par_CIBERSORTx.csv")[,c(1:10)],
                          MG = read_csv("Figure_4/Data/CIBERSORTx_estimations/mg_par_CIBERSORTx.csv")[,c(1:9)])
CIBERSORTx_nulltest_list <- lapply(CIBERSORTx_nulltest_list, function(x) setNames(melt(x, id.vars = "Mixture"), nm = c("ID","variable","beta_abs"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "CIBERSORTx")

CIBERSORTx_nulltest_list$ID <- gsub('\\.'," ",CIBERSORTx_nulltest_list$ID)

MIXTURE_nulltest_list <- mapply(FUN = MIXTURE, expressionMatrix = null_test,
                               signatureMatrix = sigs, 
                               SIMPLIFY=FALSE)

MIXTURE_nulltest_list <- lapply(MIXTURE_nulltest_list, function(x) setNames(melt(x$Subjects$MIXabs), nm = c("ID","variable","beta_abs"))) %>% 
  bind_rows(.id = "Signature") %>% 
  mutate(Method = "MIXTURE")
MIXTURE_nulltest_list[is.na(MIXTURE_nulltest_list)] <- 0

null_est <-rbind(CIBERSORTx_nulltest_list,
                 MIXTURE_nulltest_list)

#Figure

null_est_melt <- null_est

NCT_null_test <- as.data.frame(null_est_melt) %>% filter(beta_abs > 0) %>% group_by(Method, Signature,ID) %>% dplyr::summarize(count=n()) %>% ungroup() %>%
  complete(Method, nesting(ID, Signature), fill = list(count = 0))

NCT_null_test$Method <- as.factor(NCT_null_test$Method)

Figure_4 <- ggpaired(NCT_null_test, x = "Method", y = "count",
                      line.color = "gray", line.size = 0.4,
                      palette = "jco", xlab ="Method", ylab = "NCT")+
  stat_compare_means(paired = TRUE, label.y=9) 


