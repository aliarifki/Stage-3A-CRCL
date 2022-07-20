library(openxlsx)
library(dplyr)
library(ggplot2)

#Importation
genes1 = readWorkbook("E:rnaseq_melpredict/TPM_first_line.xlsx", colNames = F)
genes2 = readWorkbook("E:rnaseq_melpredict/TPM_second_line.xlsx", colNames = F)
datas<-read.table("C:/Users/aliar/Desktop/INSA/Stage/Distances/Melpred-cd8.txt", sep = ",", header = TRUE)
id_tumeurs <- na.omit(datas$id_tumeurs)
listegenes = c("CD86", "CD163", "NOS2", "MRC1", "HLA-DRA", "SIGLEC1", "IL10", "TREM2")

#Variables
genes1 = na.omit(genes1[genes1$X1== "gene_id" | genes1$X1 %in% listegenes,])
genes1 = t(genes1)
genes1 = as.data.frame(genes1)
colnames(genes1)=genes1[1,]
genes1 = genes1[!(genes1$TREM2 == 'TREM2'),]

genes2 = na.omit(genes2[genes2$X1== "gene_id" | genes2$X1 %in% listegenes,])
genes2 = t(genes2)
genes2 = as.data.frame(genes2)
colnames(genes2)=genes2[1,]
genes2 = genes2[!(genes2$TREM2 == 'TREM2'),]

data = union(genes1, genes2)

#Ajout de ZEB1
ZEB1 = readWorkbook("C:/Users/aliar/Desktop/INSA/Stage/Scatterplot 3D/reponses_traitements_melpredict.xlsx", "MELPREDICT + PAIR tot")
ZEB1 = subset(ZEB1, select = c(n_bloc, ZEB1.statut))
colnames(ZEB1) = c("gene_id", "ZEB1")
newdata = na.omit(full_join(data, ZEB1, by="gene_id"))
newdata$TREM2 = as.numeric(newdata$TREM2)
newdata$CD163 = as.numeric(newdata$CD163)
newdata$`HLA-DRA` = as.numeric(newdata$`HLA-DRA`)
newdata$CD86 = as.numeric(newdata$CD86)
newdata$SIGLEC1 = as.numeric(newdata$SIGLEC1)
newdata$MRC1 = as.numeric(newdata$MRC1)
newdata$IL10 = as.numeric(newdata$IL10)
newdata$NOS2 = as.numeric(newdata$NOS2)
newdata$ZEB1 = as.factor(newdata$ZEB1)


#TREM2
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_TREM2.png")
statistique = wilcox.test(as.numeric(newdata$TREM2)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=TREM2)) + geom_boxplot() + ggtitle(paste0("Expression de TREM2", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#CD163
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_CD163.png")
statistique = wilcox.test(as.numeric(newdata$CD163)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=CD163)) + geom_boxplot() + ggtitle(paste0("Expression de CD163", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#HLA-DRA
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_HLA-DRA.png")
statistique = wilcox.test(as.numeric(newdata$`HLA-DRA`)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=`HLA-DRA`)) + geom_boxplot() + ggtitle(paste0("Expression de HLA-DRA", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#CD86
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_CD86.png")
statistique = wilcox.test(as.numeric(newdata$CD86)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=CD86)) + geom_boxplot() + ggtitle(paste0("Expression de CD86", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#SIGLEC1
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_SIGLEC1.png")
statistique = wilcox.test(as.numeric(newdata$SIGLEC1)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=SIGLEC1)) + geom_boxplot() + ggtitle(paste0("Expression de SIGLEC1", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#MRC1
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_MRC1.png")
statistique = wilcox.test(as.numeric(newdata$MRC1)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=MRC1)) + geom_boxplot() + ggtitle(paste0("Expression de MRC1", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#IL10
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_IL10.png")
statistique = wilcox.test(as.numeric(newdata$IL10)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=IL10)) + geom_boxplot() + ggtitle(paste0("Expression de IL10", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()

#NOS2
png("C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/ZEB1_NOS2.png")
statistique = wilcox.test(as.numeric(newdata$NOS2)~as.factor(newdata$ZEB1), data=newdata)
ggplot(newdata, aes(x=ZEB1, y=NOS2)) + geom_boxplot() + ggtitle(paste0("Expression de NOS2", "  pvalue = ", signif(as.numeric(statistique["p.value"]), 2)))
dev.off()
