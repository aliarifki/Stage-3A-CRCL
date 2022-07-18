#librairies
library(phenoptr)
library(phenoptrReports)

#noms des dossiers
id_tumeurs<-c('CONT 1_Scan1', 'CONT 2_Scan1', 'CONT 3_Scan1', 'CONT 4 BIS_Scan1', 'CONT 5_Scan1', 'CONT 6_Scan1', 'CONT 9_Scan1', 'CONT 11_Scan1', 'CONT 12_Scan1', 'CONT 13_Scan1', 'CONT 14_Scan1', 'CONT 15_Scan1', 'CONT 16_Scan1', 'CONT 17_Scan1', 'CONT 18_Scan1', 'CONT 19_Scan1', 'CONT 20_Scan1', 'CONT 21_Scan1', 'CONT 22_Scan1')
i = 1
par(mar = c(10, 5, 2, 2))
#boucle
for(i in 1:length(id_tumeurs)){
  a = read_cell_seg_data(paste0("D:Analyse_Inform/", id_tumeurs[i], "/Merge_cell_seg_data.txt"))
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
