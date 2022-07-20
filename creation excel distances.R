##### Reprise du script genere par phenoptrReports
#sur les calculs de distances des cellules voisine
rm(list = ls())

library(readr)
library(dplyr)
library(remotes)
library(rtree)
library(rlang)
library(phenoptr)
library(phenoptrReports)
library(openxlsx)
library(plotly)
library(scatterplot3d)
library(xlsx)

#####
#chargement des donnees (Melpredict)
datas<-read.table("C:/Users/aliar/Desktop/INSA/Stage/Distances/Melpred-cd8.txt", sep = ",", header = TRUE)
id_tumeurs <- na.omit(datas$id_tumeurs)
phenotypes <- c("SOX10+", "SOX10+Ki67+","SOX10+ZEB1+","CD8+",
                               "CD8+PD1+","CD4+",
                               "CD4+PD1+","CD4+Ki67+")

#Reponses et zeb1
reponse = readWorkbook("C:/Users/aliar/Desktop/INSA/Stage/Scatterplot 3D/reponses_traitements_melpredict.xlsx", "MELPREDICT + PAIR tot")
reponse<-as.data.frame(cbind(reponse$Reponse.ICP.1.an, reponse$`n_bloc`, reponse$ZEB1.statut))
colnames(reponse)<-c("Rep", "ID", "ZEB1")
reponse$Rep[reponse$Rep %in% c("LR", "LR*")] <- "R"
reponse$Rep[reponse$Rep %in% c("SR", "NR ")] <- "NR"
reponse = reponse[reponse$ID %in% id_tumeurs,]
reponse = reponse[order(reponse$ID),]
id_tumeurs = id_tumeurs[id_tumeurs %in% reponse$ID]

#Separation des repondeurs et des non repondeurs
lister = reponse[reponse$Rep=="R",]
listenr = reponse[reponse$Rep=="NR",]


#Boucle pour repondeurs
tableaudistancesr = as.data.frame(0)
for(i in 1:length(lister$ID)){
  distances = readWorkbook(paste0("D:melpredict_cd8/", lister$ID[i], "/Results.xlsx"), "Nearest Neighbors")
  distances<-as.data.frame(cbind(distances[2], distances$X3, distances$X4, distances$X6))
  colnames(distances) = c("type", "p1", "p2", "mean")
  distances = distances[distances$type=="Tumor",]
  tableaudistancesr = bind_rows(tableaudistancesr, distances)
#  dr = c(dr, distances$mean[(distances$p1==phenotype1 & distances$p2==phenotype2) | (distances$p2==phenotype1 & distances$p1==phenotype2)])
}

dr = c()
datar = as.data.frame(0)
for(i in 1:length(phenotypes)){
  for(j in 1:length(phenotypes)){
    dr = mean(na.omit(as.numeric(tableaudistancesr$mean[distances$p1==phenotypes[i] & distances$p2==phenotypes[j]])))
    ligner = cbind.data.frame(phenotypes[i], phenotypes[j], dr)
    datar = bind_rows(datar,ligner)
  }
}

#Boucle pour repondeurs
tableaudistancesnr = as.data.frame(0)
for(i in 1:length(listenr$ID)){
  distances = readWorkbook(paste0("D:melpredict_cd8/", listenr$ID[i], "/Results.xlsx"), "Nearest Neighbors")
  distances<-as.data.frame(cbind(distances[2], distances$X3, distances$X4, distances$X6))
  colnames(distances) = c("type", "p1", "p2", "mean")
  distances = distances[distances$type=="Tumor",]
  tableaudistancesnr = bind_rows(tableaudistancesnr, distances)
  #  dr = c(dr, distances$mean[(distances$p1==phenotype1 & distances$p2==phenotype2) | (distances$p2==phenotype1 & distances$p1==phenotype2)])
}

dnr = c()
datanr = as.data.frame(0)
for(i in 1:length(phenotypes)){
  for(j in 1:length(phenotypes)){
    dnr = mean(na.omit(as.numeric(tableaudistancesnr$mean[distances$p1==phenotypes[i] & distances$p2==phenotypes[j]])))
    lignenr = cbind.data.frame(phenotypes[i], phenotypes[j], dnr)
    datanr = bind_rows(datanr,lignenr)
  }
}

#Suppression des na parasites
datar = subset(datar, select=-`0`)
datar = datar[-1,]
datanr = subset(datanr, select=-`0`)
datanr = datanr[-1,]

#Jointure
tableaudedonnees = cbind(datar, subset(datanr, select=dnr))
colnames(tableaudedonnees) = c("Phenotype 1", "Phenotype 2", "Repondeurs", "Non repondeurs")


#Excel

wb = createWorkbook("C:/Users/aliar/Desktop/INSA/Stage/Distances_R_NR.xlsx")
addWorksheet(wb, "Moyenne des distances")
write.xlsx(x = reponse, file = "C:/Users/aliar/Desktop/INSA/Stage/Distances_R_NR.xlsx", sheetName = "Statut de réponse", append=T)
