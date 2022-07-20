library(openxlsx)

#Importation
genes1 = readWorkbook("E:rnaseq_melpredict/TPM_first_line.xlsx")
genes2 = readWorkbook("E:rnaseq_melpredict/TPM_second_line.xlsx")
proportions = readWorkbook("C:/Users/aliar/Desktop/INSA/Stage/decompte_phenotypes.xlsx", "Cell prop (tumor and stroma)")
datas<-read.table("C:/Users/aliar/Desktop/INSA/Stage/Distances/Melpred-cd8.txt", sep = ",", header = TRUE)
id_tumeurs <- na.omit(datas$id_tumeurs)
path <- na.omit(datas$path)
marqueurs <- na.omit(datas$marqueurs)
opal <- na.omit(datas$opal)
zone_marqueurs <- na.omit(datas$zone_marqueurs)
dossier_data <- na.omit(datas$dossier_data)
phenotypes <- c("SOX10+", "SOX10+Ki67+","SOX10+ZEB1+","CD8+",
                               "CD8+PD1+","CD4+",
                               "CD4+PD1+","CD4+Ki67+", "Total Cells")

#Variables
genes1 = na.omit(genes1[genes1$GUCA2B== "gene_id" | genes1$GUCA2B == "TCF7" | genes1$GUCA2B == "TOX" | genes1$GUCA2B == "HAVCR2" | genes1$GUCA2B == "LAG3" | genes1$GUCA2B == "CTLA4", ])
row.names(genes1) = c("LAG3", "HAVCR2", "CTLA4", "TCF7", "TOX", "ID")
genes1 = t(genes1)
genes1 = as.data.frame(genes1)

#Sélection 
genes1 = genes1[genes1$ID %in% proportions$ID_tumeur,]
proportions = proportions[proportions$ID_tumeur %in% genes1$ID,]

#Ordre
genes1 = genes1[order(genes1$ID),]
proportions = proportions[order(proportions$ID_tumeur),]

#Liaison
data = cbind(genes1, subset(proportions, select = CD8.PD1.))

#Plot
png(filename = "C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/LAG3.png")
statistiques = cor.test(as.numeric(data$LAG3),as.numeric(data$CD8.PD1.), method = "pearson")
plot(data$LAG3, data$CD8.PD1., xlab = "Expression de LAG3", ylab = "Proportion de CD8+PD1+", main = paste0("R = ", round(statistiques$estimate, 2), "  pvalue = ", signif(statistiques$p.value, 3)))
dev.off()

png(filename = "C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/HAVCR2.png")
statistiques = cor.test(as.numeric(data$HAVCR2),as.numeric(data$CD8.PD1.), method = "pearson")
plot(data$HAVCR2, data$CD8.PD1., xlab = "Expression de HAVCR2", ylab = "Proportion de CD8+PD1+", main = paste0("R = ", round(statistiques$estimate, 2), "  pvalue = ", signif(statistiques$p.value, 3)))
dev.off()

png(filename = "C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/CTLA4.png")
statistiques = cor.test(as.numeric(data$CTLA4),as.numeric(data$CD8.PD1.), method = "pearson")
plot(data$CTLA4, data$CD8.PD1., xlab = "Expression de CTLA4", ylab = "Proportion de CD8+PD1+", main = paste0("R = ", round(statistiques$estimate, 2), "  pvalue = ", signif(statistiques$p.value, 3)))
dev.off()

png(filename = "C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/TOX.png")
statistiques = cor.test(as.numeric(data$TOX),as.numeric(data$CD8.PD1.), method = "pearson")
plot(data$TOX, data$CD8.PD1., xlab = "Expression de TOX", ylab = "Proportion de CD8+PD1+", main = paste0("R = ", round(statistiques$estimate, 2), "  pvalue = ", signif(statistiques$p.value, 3)))
dev.off()

png(filename = "C:/Users/aliar/Desktop/INSA/Stage/Expression de genes/TCF7.png")
statistiques = cor.test(as.numeric(data$TCF7),as.numeric(data$CD8.PD1.), method = "pearson")
plot(data$TCF7, data$CD8.PD1., xlab = "Expression de TCF7", ylab = "Proportion de CD8+PD1+", main = paste0("R = ", round(statistiques$estimate, 2), "  pvalue = ", signif(statistiques$p.value, 3)))
dev.off()

