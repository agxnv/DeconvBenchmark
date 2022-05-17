library(openxlsx)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

K_sig <- read.csv("Figure_1/Data/MS/K.csv", row.names = 1)
K_sig <- K_sig[,order(colnames(K_sig))]
colnames(K_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
K_cormat <- round(cor(K_sig),2)

K_melted <- melt(K_cormat)
K_melted$Var2 = with(K_melted, factor(Var2, levels = rev(levels(Var2))))

K_corplot <- ggplot(data = K_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Kidney") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

L_sig <- read.csv("Figure_1/Data/MS/L.csv", row.names = 1)
L_sig <- L_sig[,order(colnames(L_sig))]
colnames(L_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
L_cormat <- round(cor(L_sig),2)

L_melted <- melt(L_cormat)
L_melted$Var2 = with(L_melted, factor(Var2, levels = rev(levels(Var2))))

L_corplot <- ggplot(data = L_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Liver") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

MG_sig <- read.csv("Figure_1/Data/MS/MG.csv", row.names = 1)
MG_sig <- MG_sig[,order(colnames(MG_sig))]
colnames(MG_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
MG_cormat <- round(cor(MG_sig),2)

MG_melted <- melt(MG_cormat)
MG_melted$Var2 = with(MG_melted, factor(Var2, levels = rev(levels(Var2))))

MG_corplot <- ggplot(data = MG_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Mammary gland") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

M_sig <- read.csv("Figure_1/Data/MS/M.csv", row.names = 1)
M_sig <- M_sig[,order(colnames(M_sig))]
colnames(M_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
M_cormat <- round(cor(M_sig),2)

M_melted <- melt(M_cormat)
M_melted$Var2 = with(M_melted, factor(Var2, levels = rev(levels(Var2))))

M_corplot <- ggplot(data = M_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Muscle") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

SI_sig <- read.csv("Figure_1/Data/MS/SI.csv", row.names = 1)
SI_sig <- SI_sig[,order(colnames(SI_sig))]
colnames(SI_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
SI_cormat <- round(cor(SI_sig),2)

SI_melted <- melt(SI_cormat)
SI_melted$Var2 = with(SI_melted, factor(Var2, levels = rev(levels(Var2))))

SI_corplot <- ggplot(data = SI_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Small intestine") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

S_sig <- read.csv("Figure_1/Data/MS/S.csv", row.names = 1)
S_sig <- S_sig[,order(colnames(S_sig))]
colnames(S_sig) <- c("B","Bsp","Dnd","Mcr","Mnc","PMN","NK","T")
S_cormat <- round(cor(S_sig),2)

S_melted <- melt(S_cormat)
S_melted$Var2 = with(S_melted, factor(Var2, levels = rev(levels(Var2))))

S_corplot <- ggplot(data = S_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_y_discrete(position = "right", label=abbreviate) + ggtitle("Spleen") +
  xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Correlation\nCoefficient", reverse=T)) + scale_x_discrete(label=abbreviate) + scale_fill_continuous(high = "red", low = "light blue")

############Sup_Table1

K_cormat[lower.tri(K_cormat, diag = TRUE)] <- NA
K_melted_table <- melt(K_cormat, na.rm = TRUE)

L_cormat[lower.tri(L_cormat, diag = TRUE)] <- NA
L_melted_table <- melt(L_cormat, na.rm = TRUE)

MG_cormat[lower.tri(MG_cormat, diag = TRUE)] <- NA
MG_melted_table <- melt(MG_cormat, na.rm = TRUE)

M_cormat[lower.tri(M_cormat, diag = TRUE)] <- NA
M_melted_table <- melt(M_cormat, na.rm = TRUE)

SI_cormat[lower.tri(SI_cormat, diag = TRUE)] <- NA
SI_melted_table <- melt(SI_cormat, na.rm = TRUE)

S_cormat[lower.tri(S_cormat, diag = TRUE)] <- NA
S_melted_table <- melt(S_cormat, na.rm = TRUE)

K_melted_table$signature <- "Kidney"
L_melted_table$signature <- "Liver"
MG_melted_table$signature <- "Mammary gland"
M_melted_table$signature <- "Muscle"
SI_melted_table$signature <- "Small intestine"
S_melted_table$signature <- "Spleen"

cortotal <- rbind(K_melted_table,L_melted_table,MG_melted_table,M_melted_table,PB_melted_table,SI_melted_table,S_melted_table) 

#define quantiles of interest
q = c(.25, .5, .75)

Sup_Table1 <- cortotal %>%
  group_by(signature) %>%
  summarize(quant25 = quantile(value, probs = q[1]), 
            quant50 = quantile(value, probs = q[2]),
            quant75 = quantile(value, probs = q[3]),
            min = min(value),
            max = max(value),
            mean = mean(value),
            sd = sd(value)) %>% 
  ungroup()

