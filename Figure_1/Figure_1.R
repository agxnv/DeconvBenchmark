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

sigs <- list(BM = read.csv("Figure_1/Data/MS/BM.csv", row.names = 1, check.names=FALSE),
             K = read.csv("Figure_1/Data/MS/K.csv", row.names = 1, check.names=FALSE),
             L = read.csv("Figure_1/Data/MS/L.csv", row.names = 1, check.names=FALSE),
             LU = read.csv("Figure_1/Data/MS/LU.csv", row.names = 1, check.names=FALSE),
             MG = read.csv("Figure_1/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

#CIBERSORTx

CIBERSORTx_list <- list(BM = read_table2("Figure_1/Data/CIBERSORTx_estimations/Bmarrow.CIBERSORTx.txt")[,c(1:11)],
                       K = read_table2("Figure_1/Data/CIBERSORTx_estimations/Kidney.CIBERSORTx.txt")[,c(1:9)],
                       L = read_table2("Figure_1/Data/CIBERSORTx_estimations/Liver.CIBERSORTx.txt")[,c(1:9)],
                       LU = read_table2("Figure_1/Data/CIBERSORTx_estimations/Lung.CIBERSORTx.txt")[,c(1:10)],
                       MG = read_table2("Figure_1/Data/CIBERSORTx_estimations/Mgland.CIBERSORTx.txt")[,c(1:9)])

CIBERSORTx_list <- lapply(CIBERSORTx_list, function(x) as.data.frame(x) %>%
                         arrange(Mixture) %>% na_if(0)) 

CIBERSORTx_list <- lapply(CIBERSORTx_list, function(x) data.matrix(x[,!(colnames(x) == "Mixture")][,order(colnames(x[,!(colnames(x) == "Mixture")]))]))

cell.types.names.k.l.mg <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
cell.types.names.lu <- c("B","Bsp","Dnd","Eo","Mcr","Mnc","PMN","NK","T")
cell.types.names.bm <- c("B","Bsp","Dnd","Eo","Mcr","Mst", "Mnc","PMN","NK","T")

colnames(CIBERSORTx_list$BM) <- rownames(CIBERSORTx_list$BM) <- cell.types.names.bm
colnames(CIBERSORTx_list$K) <- rownames(CIBERSORTx_list$K) <- cell.types.names.k.l.mg
colnames(CIBERSORTx_list$L) <- rownames(CIBERSORTx_list$L) <- cell.types.names.k.l.mg
colnames(CIBERSORTx_list$LU) <- rownames(CIBERSORTx_list$LU) <- cell.types.names.lu
colnames(CIBERSORTx_list$MG) <- rownames(CIBERSORTx_list$MG) <- cell.types.names.k.l.mg

DCQ_list <- mapply(FUN = estCellPercent.DCQ,refExpr = sigs, geneExpr = sigs,SIMPLIFY=FALSE)
DCQ_list <- lapply(DCQ_list, function(x) as.data.frame(x[!(row.names(x) == "others"), ]/100) %>% 
                     mutate(X1 = rownames(x[!(row.names(x) == "others"), ])) %>% 
                     arrange(X1) %>% na_if(0)) 

DCQ_list <- lapply(DCQ_list, function(x) data.matrix(x[,!(colnames(x) == "X1")][,order(colnames(x[,!(colnames(x) == "X1")]))]))

colnames(DCQ_list$BM) <- rownames(DCQ_list$BM) <- cell.types.names.bm
colnames(DCQ_list$K) <- rownames(DCQ_list$K) <- cell.types.names.k.l.mg
colnames(DCQ_list$L) <- rownames(DCQ_list$L) <- cell.types.names.k.l.mg
colnames(DCQ_list$LU) <- rownames(DCQ_list$LU) <- cell.types.names.lu
colnames(DCQ_list$MG) <- rownames(DCQ_list$MG) <- cell.types.names.k.l.mg

DeconRNASeq_list <- mapply(FUN = estCellPercent.DeconRNASeq,refExpr = sigs, geneExpr = sigs,SIMPLIFY=FALSE)


DeconRNASeq_list <- lapply(DeconRNASeq_list, function(x) as.data.frame(x[!(row.names(x) == "others"), ]/100) %>% 
                             mutate(X1 = rownames(x[!(row.names(x) == "others"), ])) %>% 
                             arrange(X1) %>% na_if(0)) 

DeconRNASeq_list <- lapply(DeconRNASeq_list, function(x) data.matrix(x[,!(colnames(x) == "X1")][,order(colnames(x[,!(colnames(x) == "X1")]))]))

colnames(DeconRNASeq_list$BM) <- rownames(DeconRNASeq_list$BM) <- cell.types.names.bm
colnames(DeconRNASeq_list$K) <- rownames(DeconRNASeq_list$K) <- cell.types.names.k.l.mg
colnames(DeconRNASeq_list$L) <- rownames(DeconRNASeq_list$L) <- cell.types.names.k.l.mg
colnames(DeconRNASeq_list$LU) <- rownames(DeconRNASeq_list$LU) <- cell.types.names.lu
colnames(DeconRNASeq_list$MG) <- rownames(DeconRNASeq_list$MG) <- cell.types.names.k.l.mg

#EPIC

EPIC_list <- mapply(FUN = EPIC::EPIC, bulk = sigs, reference = lapply(sigs, function(x) list(refProfiles = x, sigGenes = rownames(x))), withOtherCells = FALSE, 
               SIMPLIFY=FALSE)

EPIC_list <- lapply(EPIC_list, function(x) x$cellFractions)

EPIC_list <- lapply(EPIC_list, function(x) as.data.frame(x) %>% 
                             mutate(X1 = rownames(x)) %>% 
                             arrange(X1) %>% na_if(0)) 

EPIC_list <- lapply(EPIC_list, function(x) data.matrix(x[,!(colnames(x) == "X1")][,order(colnames(x[,!(colnames(x) == "X1")]))]))

colnames(EPIC_list$BM) <- rownames(EPIC_list$BM) <- cell.types.names.bm
colnames(EPIC_list$K) <- rownames(EPIC_list$K) <- cell.types.names.k.l.mg
colnames(EPIC_list$L) <- rownames(EPIC_list$L) <- cell.types.names.k.l.mg
colnames(EPIC_list$LU) <- rownames(EPIC_list$LU) <- cell.types.names.lu
colnames(EPIC_list$MG) <- rownames(EPIC_list$MG) <- cell.types.names.k.l.mg

#MIXTURE
MIXTURE_list <- mapply(FUN = MIXTURE, expressionMatrix = sigs,
                       signatureMatrix = sigs, 
                    SIMPLIFY=FALSE)

MIXTURE_list <- lapply(MIXTURE_list, function(x) x$Subjects$MIXprop)

MIXTURE_list <- lapply(MIXTURE_list, function(x) as.data.frame(x) %>% 
                      mutate(X1 = rownames(x)) %>% 
                      arrange(X1) %>% na_if(0)) 

MIXTURE_list <- lapply(MIXTURE_list, function(x) data.matrix(x[,!(colnames(x) == "X1")][,order(colnames(x[,!(colnames(x) == "X1")]))]))

colnames(MIXTURE_list$BM) <- rownames(MIXTURE_list$BM) <- cell.types.names.bm
colnames(MIXTURE_list$K) <- rownames(MIXTURE_list$K) <- cell.types.names.k.l.mg
colnames(MIXTURE_list$L) <- rownames(MIXTURE_list$L) <- cell.types.names.k.l.mg
colnames(MIXTURE_list$LU) <- rownames(MIXTURE_list$LU) <- cell.types.names.lu
colnames(MIXTURE_list$MG) <- rownames(MIXTURE_list$MG) <- cell.types.names.k.l.mg


#QUANTISEQ

rownames(sigs$BM) <- toupper(rownames(sigs$BM))
rownames(sigs$K) <- toupper(rownames(sigs$K))
rownames(sigs$L) <- toupper(rownames(sigs$L))
rownames(sigs$LU) <- toupper(rownames(sigs$LU))
rownames(sigs$MG) <- toupper(rownames(sigs$MG))


quanTIseq_list <- mapply(FUN = deconvolute_quantiseq.default,
                         mix.mat = sigs, 
                         arrays = FALSE, 
                         signame = paste0("Figure_1/Data/MS/",names(sigs)),
                         tumor = FALSE, 
                         mRNAscale = FALSE, method = "lsei", btotalcells = FALSE,
                         rmgenes = "unassigned", 
                    SIMPLIFY=FALSE)

quanTIseq_list <- lapply(quanTIseq_list, function(x) x %>% select(-c("Other")))

quanTIseq_list <- lapply(quanTIseq_list, function(x) as.data.frame(x) %>%  
                      arrange(Sample) %>% na_if(0)) 

quanTIseq_list <- lapply(quanTIseq_list, function(x) data.matrix(x[,!(colnames(x) == "Sample")][,order(colnames(x[,!(colnames(x) == "Sample")]))]))

colnames(quanTIseq_list$BM) <- rownames(quanTIseq_list$BM) <- cell.types.names.bm
colnames(quanTIseq_list$K) <- rownames(quanTIseq_list$K) <- cell.types.names.k.l.mg
colnames(quanTIseq_list$L) <- rownames(quanTIseq_list$L) <- cell.types.names.k.l.mg
colnames(quanTIseq_list$LU) <- rownames(quanTIseq_list$LU) <- cell.types.names.lu
colnames(quanTIseq_list$MG) <- rownames(quanTIseq_list$MG) <- cell.types.names.k.l.mg


sigs_cor_5 <- list(K = read.csv("Figure_1/Data/MS/K.csv", row.names = 1, check.names=FALSE),
                   L = read.csv("Figure_1/Data/MS/L.csv", row.names = 1, check.names=FALSE),
                   MG = read.csv("Figure_1/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

sigs_cor_5 <- lapply(sigs_cor_5, function(x) {x <- x[,order(colnames(x))];
colnames(x) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T");
x <- as.matrix(round(cor(x),2));
x})

lu_cor <- read.csv("Figure_1/Data/MS/LU.csv", row.names = 1, check.names=FALSE)
lu_cor <- lu_cor[,order(colnames(lu_cor))]
colnames(lu_cor) <- c("B","Bsp","Dnd","Eo","Mcr","Mnc","PMN","NK","T")
lu_cor <- round(cor(lu_cor),2)

bm_cor <- read.csv("Figure_1/Data/MS/BM.csv", row.names = 1, check.names=FALSE)
bm_cor <- bm_cor[,order(colnames(bm_cor))]
colnames(bm_cor) <- c("B","Bsp","Dnd","Eo","Mcr","Mst", "Mnc","PMN","NK","T")
bm_cor <- round(cor(bm_cor),2)

#Heatmaps
library(ComplexHeatmap)
library(circlize)

BMarrow_heatmap <- 
  Heatmap(bm_cor, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "Correlation",name = "Estimated \nCoefficients",row_title = "Bone Marrow",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("light blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(CIBERSORTx_list$BM, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORTx",name = "Estimated \nCoefficients",row_title = "Bone Marrow",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(DCQ_list$BM, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DCQ",name = "DCQ",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(DeconRNASeq_list$BM, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DRNASeq",name = "DRNASeq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(EPIC_list$BM, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MIXTURE_list$BM, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(quanTIseq_list$BM, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

Kidney_heatmap <- 
  Heatmap(sigs_cor_5$K, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "Correlation",name = "Estimated \nCoefficients",row_title = "Kidney",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("light blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(CIBERSORTx_list$K, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORTx",name = "Estimated \nCoefficients",row_title = "Kidney",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(DCQ_list$K, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DCQ",name = "DCQ",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(DeconRNASeq_list$K, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DRNASeq",name = "DRNASeq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(EPIC_list$K, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MIXTURE_list$K, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(quanTIseq_list$K, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

Liver_heatmap <- 
  Heatmap(sigs_cor_5$L, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "Correlation",name = "Estimated \nCoefficients",row_title = "Liver",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("light blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(CIBERSORTx_list$L, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORTx",name = "Estimated \nCoefficients",row_title = "Liver",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(DCQ_list$L, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DCQ",name = "DCQ",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(DeconRNASeq_list$L, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DRNASeq",name = "DRNASeq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(EPIC_list$L, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MIXTURE_list$L, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(quanTIseq_list$L, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

Lung_heatmap <- 
  Heatmap(lu_cor, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "Correlation",name = "Estimated \nCoefficients",row_title = "Lung",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("light blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(CIBERSORTx_list$LU, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORTx",name = "Estimated \nCoefficients",row_title = "Lung",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(DCQ_list$LU, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DCQ",name = "DCQ",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(DeconRNASeq_list$LU, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DRNASeq",name = "DRNASeq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(EPIC_list$LU, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MIXTURE_list$LU, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(quanTIseq_list$LU, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

MGland_heatmap <- 
  Heatmap(sigs_cor_5$MG, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "Correlation",name = "Estimated \nCoefficients",row_title = "Mammary gland",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("light blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(CIBERSORTx_list$MG, cluster_rows = FALSE, show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "CIBERSORTx",name = "Estimated \nCoefficients",row_title = "Mammary gland",row_title_gp = gpar(fontsize=15), column_title_gp = gpar(fontsize=15),col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE) +
  Heatmap(DCQ_list$MG, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DCQ",name = "DCQ",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(DeconRNASeq_list$MG, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "DRNASeq",name = "DRNASeq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(EPIC_list$MG, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE,column_title = "EPIC",name = "EPIC",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F) +
  Heatmap(MIXTURE_list$MG, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = FALSE, cluster_columns = FALSE, column_title = "MIXTURE",name = "MIXTURE",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = FALSE, show_column_names = TRUE) +
  Heatmap(quanTIseq_list$MG, cluster_rows = FALSE, column_title_gp = gpar(fontsize=15), show_row_names = TRUE, cluster_columns = FALSE, column_title = "quanTIseq",name = "quanTIseq",col = colorRamp2(c(0, 1), c("blue",  "red")),show_heatmap_legend = F, show_column_names = TRUE)

bm_cor[lower.tri(bm_cor, diag = TRUE)] <- NA
BM_melted_table <- melt(bm_cor, na.rm = TRUE)
BM_melted_table$signature <- "Bone marrow"
sigs_cor_5$K[lower.tri(sigs_cor_5$K, diag = TRUE)] <- NA
K_melted_table <- melt(sigs_cor_5$K, na.rm = TRUE)
K_melted_table$signature <- "Kidney"
sigs_cor_5$L[lower.tri(sigs_cor_5$L, diag = TRUE)] <- NA
L_melted_table <- melt(sigs_cor_5$L, na.rm = TRUE)
L_melted_table$signature <- "Liver"
lu_cor[lower.tri(lu_cor, diag = TRUE)] <- NA
LU_melted_table <- melt(lu_cor, na.rm = TRUE)
LU_melted_table$signature <- "Lung"
sigs_cor_5$MG[lower.tri(sigs_cor_5$MG, diag = TRUE)] <- NA
MG_melted_table <- melt(sigs_cor_5$MG, na.rm = TRUE)
MG_melted_table$signature <- "Mammary gland"

cortotal <- rbind(K_melted_table,L_melted_table,MG_melted_table,BM_melted_table,LU_melted_table) %>% arrange(signature)

q = c(.25, .5, .75)

Sup_Table1 <- cortotal %>%
  group_by(signature) %>%
  dplyr::summarize(quant25 = quantile(value, probs = q[1]), 
            quant50 = quantile(value, probs = q[2]),
            quant75 = quantile(value, probs = q[3]),
            min = min(value),
            max = max(value),
            mean = mean(value),
            sd = sd(value)) %>% 
  ungroup()


#Sup figure 1

DCQ_list_supfig1 <- lapply(DCQ_list, function(x) melt(x))
DCQ_list_supfig1 <-  DCQ_list_supfig1 %>% bind_rows(.id = "Signature")
DCQ_list_supfig1 <- DCQ_list_supfig1[DCQ_list_supfig1$Var1 != DCQ_list_supfig1$Var2,]
DCQ_list_supfig1 <- DCQ_list_supfig1[complete.cases(DCQ_list_supfig1),]
DCQ_list_supfig1$Method <- "DCQ"

CIBERSORTx_list_supfig1 <- lapply(CIBERSORTx_list, function(x) melt(x))
CIBERSORTx_list_supfig1 <-  CIBERSORTx_list_supfig1 %>% bind_rows(.id = "Signature")
CIBERSORTx_list_supfig1 <- CIBERSORTx_list_supfig1[CIBERSORTx_list_supfig1$Var1 != CIBERSORTx_list_supfig1$Var2,]
CIBERSORTx_list_supfig1 <- CIBERSORTx_list_supfig1[complete.cases(CIBERSORTx_list_supfig1),]
CIBERSORTx_list_supfig1$Method <- "CIBERSORTx"

EPIC_list_supfig1 <- lapply(EPIC_list, function(x) melt(x))
EPIC_list_supfig1 <-  EPIC_list_supfig1 %>% bind_rows(.id = "Signature")
EPIC_list_supfig1 <- EPIC_list_supfig1[EPIC_list_supfig1$Var1 != EPIC_list_supfig1$Var2,]
EPIC_list_supfig1 <- EPIC_list_supfig1[complete.cases(EPIC_list_supfig1),]
EPIC_list_supfig1$Method <- "EPIC"

quanTIseq_list_supfig1 <- lapply(quanTIseq_list, function(x) melt(x))
quanTIseq_list_supfig1 <-  quanTIseq_list_supfig1 %>% bind_rows(.id = "Signature")
quanTIseq_list_supfig1 <- quanTIseq_list_supfig1[quanTIseq_list_supfig1$Var1 != quanTIseq_list_supfig1$Var2,]
quanTIseq_list_supfig1 <- quanTIseq_list_supfig1[complete.cases(quanTIseq_list_supfig1),]
quanTIseq_list_supfig1$Method <- "quanTIseq"

sigs_cor_5 <- list(K = read.csv("Figure_1/Data/MS/K.csv", row.names = 1, check.names=FALSE),
                   L = read.csv("Figure_1/Data/MS/L.csv", row.names = 1, check.names=FALSE),
                   MG = read.csv("Figure_1/Data/MS/MG.csv", row.names = 1, check.names=FALSE))

sigs_cor_5 <- lapply(sigs_cor_5, function(x) {x <- x[,order(colnames(x))];
colnames(x) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T");
x <- as.matrix(round(cor(x),2));
x})

lu_cor <- read.csv("Figure_1/Data/MS/LU.csv", row.names = 1, check.names=FALSE)
lu_cor <- lu_cor[,order(colnames(lu_cor))]
colnames(lu_cor) <- c("B","Bsp","Dnd","Eo","Mcr","Mnc","PMN","NK","T")
lu_cor <- round(cor(lu_cor),2)

bm_cor <- read.csv("Figure_1/Data/MS/BM.csv", row.names = 1, check.names=FALSE)
bm_cor <- bm_cor[,order(colnames(bm_cor))]
colnames(bm_cor) <- c("B","Bsp","Dnd","Eo","Mcr","Mst", "Mnc","PMN","NK","T")
bm_cor <- round(cor(bm_cor),2)

bm_cor[lower.tri(bm_cor, diag = TRUE)] <- NA
BM_melted_table <- melt(bm_cor, na.rm = TRUE, value.name = "cor")
BM_melted_table$Signature <- "BM"
sigs_cor_5$K[lower.tri(sigs_cor_5$K, diag = TRUE)] <- NA
K_melted_table <- melt(sigs_cor_5$K, na.rm = TRUE, value.name = "cor")
K_melted_table$Signature <- "K"
sigs_cor_5$L[lower.tri(sigs_cor_5$L, diag = TRUE)] <- NA
L_melted_table <- melt(sigs_cor_5$L, na.rm = TRUE, value.name = "cor")
L_melted_table$Signature <- "L"
lu_cor[lower.tri(lu_cor, diag = TRUE)] <- NA
LU_melted_table <- melt(lu_cor, na.rm = TRUE, value.name = "cor")
LU_melted_table$Signature <- "LU"
sigs_cor_5$MG[lower.tri(sigs_cor_5$MG, diag = TRUE)] <- NA
MG_melted_table <- melt(sigs_cor_5$MG, na.rm = TRUE, value.name = "cor")
MG_melted_table$Signature <- "MG"

cortotal <- rbind(K_melted_table,L_melted_table,MG_melted_table,BM_melted_table,LU_melted_table) %>% arrange(Signature)
cortotal[166:330,] <- cortotal[,c(2,1,3,4)]

Supfig1_df <- rbind(EPIC_list_supfig1,CIBERSORTx_list_supfig1,DCQ_list_supfig1,quanTIseq_list_supfig1)
Supfig1_df<- left_join(Supfig1_df,
                       cortotal,
                       by = c("Signature","Var1","Var2"))

Sup_fig_1 <- ggplot(Supfig1_df, aes(x=cor, y=value)) +
  geom_point() +
  facet_grid(~Method) +
  xlab("Correlation between the deconvolved MS profile and the detected cell type") +
  ylab("Estimated coefficient")+ stat_smooth(method="lm", se=FALSE, color="red")
