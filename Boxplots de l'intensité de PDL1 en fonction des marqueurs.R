######
#Librairies
library(phenoptr)
library(phenoptrReports)
library(ggplot2)
library(openxlsx)
library(dplyr)

######
#Noms des dossiers
datas<-read.table("C:/Users/aliar/Desktop/INSA/Stage/Distances/Melpred-cd8.txt", sep = ",", header = TRUE)
id_tumeurs <- na.omit(datas$id_tumeurs)
par(mar = c(10, 5, 2, 2))

######
#Tous les phénotypes
for(i in 1:length(id_tumeurs)){
  a = read_cell_seg_data(paste0("D:melpredict_cd8/", id_tumeurs[i], "/Merge_cell_seg_data.txt"))
  a = a[a$`Tissue Category`=="Tumor", ]
  #a = a[!(a$Phenotype %in% c("CD68+")),]
  PDL1 = a$`Entire Cell Opal 620 Mean`
  phenotype = as.factor(a$Phenotype)
  jpeg(paste0("C:/Users/aliar/Desktop/INSA/Stage/Boxplot valentin/Boxplot_", id_tumeurs[i], ".jpg")) #début de l'enregistrement
  plot(phenotype, PDL1, xlab = "Phénotype", ylab = "Intensité de PDL1", las=2, font.lab=2) #boxplot
  dev.off() #fin de l'enregistrement
}

phenotype

head(a$Phenotype)

getwd()
table(phenotype)

i=5

reponse = readWorkbook("C:/Users/aliar/Desktop/INSA/Stage/Scatterplot 3D/reponses_traitements_melpredict.xlsx", "MELPREDICT + PAIR tot")
reponse<-as.data.frame(cbind(reponse$Reponse.ICP.1.an, reponse$`n_bloc`, reponse$ZEB1.statut))
colnames(reponse)<-c("Rep", "ID", "ZEB1")
reponse$Rep[reponse$Rep %in% c("LR", "LR*")] <- "R"
reponse$Rep[reponse$Rep %in% c("SR", "NR ")] <- "NR"
reponse = reponse[reponse$ID %in% id_tumeurs,]
reponse = reponse[order(reponse$ID),]
id_tumeurs = id_tumeurs[id_tumeurs %in% reponse$ID]

a = as.data.frame(0)
b = as.data.frame(0)
######
#Juste CD68, CD163CD68+/-PDL1, CD68pSTAT1+/-PDL1
for(i in 1:length(id_tumeurs)){
  a = read_cell_seg_data(paste0("D:melpredict_cd8/", id_tumeurs[i], "/Merge_cell_seg_data.txt"))
  a = a[a$`Tissue Category`=="Tumor", ]
  #a = a[!(a$Phenotype %in% c("CD68+")),]
  a = subset(a, select = c(Phenotype, `Entire Cell Opal 620 Mean`))
  a$Groupe = 1
  a = a[(a$Phenotype == "CD68+" |a$Phenotype == "CD68+CD163+" |a$Phenotype == "CD68+CD163+PDL1+" |a$Phenotype == "CD68+pSTAT1+PDL1+" |a$Phenotype == "CD68+pSTAT1+"), ]
  a$ID = id_tumeurs[i]
  b = bind_rows(b, a)
}

data = full_join(b, reponse, by="ID")
data = subset(data, select = -`0`)
data$Groupe[data$Phenotype=="CD68+"] = "CD68+"
data$Groupe[data$Phenotype=="CD68+CD163+"|data$Phenotype=="CD68+CD163+PDL1+"] = "CD68+CD163+"
data$Groupe[data$Phenotype=="CD68+pSTAT1+PDL1+"|data$Phenotype=="CD68+pSTAT1+PDL1+"] = "CD68+pSTAT1+"
data = na.omit(data)
PDL1 = data$`Entire Cell Opal 620 Mean`
data$Groupe = as.factor(data$Groupe)
Rep = as.factor(data$Rep)
ZEB1 = as.factor(data$ZEB1)


#R NR
ggplot(data = data, aes(x=Groupe, y=PDL1, fill=Rep)) +
  geom_boxplot() +
  labs(x = "Phénotype", fill = "Réponse")

#R NR
ggplot(data = data, aes(x=Groupe, y=PDL1, fill=ZEB1)) +
  geom_boxplot() +
  labs(x = "Phénotype", fill = "ZEB1")


png(paste0("C:/Users/aliar/Desktop/INSA/Stage/Boxplot valentin/Boxplot_", id_tumeurs[i], ".png")) #début de l'enregistrement
#ggplot à mettre ici
dev.off() #fin de l'enregistrement