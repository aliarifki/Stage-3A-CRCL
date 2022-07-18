install.packages('remotes')
library(remotes)

remotes::install_github("akoyabio/phenoptrReports")
install.packages('devtools')
library(devtools)
devtools::install_github('akoyabio/phenoptr')

library(phenoptrReports)
library(reshape)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(corrplot)
library(rlist)
library(dplyr)
require(akima)
library(latticeExtra)
library(phenoptr)
library(tidyverse)
library(rtree)
library(XLConnect)
library(stringr)
rm(list = ls())

id_tumeurs<-c('CONT 1_Scan1', 'CONT 2_Scan1', 'CONT 3_Scan1', 'CONT 4 BIS_Scan1', 'CONT 5_Scan1', 'CONT 6_Scan1', 'CONT 9_Scan1', 'CONT 11_Scan1', 'CONT 12_Scan1', 'CONT 13_Scan1', 'CONT 14_Scan1', 'CONT 15_Scan1', 'CONT 16_Scan1', 'CONT 17_Scan1', 'CONT 18_Scan1', 'CONT 19_Scan1', 'CONT 20_Scan1', 'CONT 21_Scan1', 'CONT 22_Scan1','DESMO 1_Scan1', 'DESMO 4_Scan1', 'DESMO 5_Scan1',  'DESMO 8_Scan1', 'DESMO 9_Scan1', 'DESMO 11_Scan1', 'DESMO 12_Scan1', 'DESMO 13 BIS_Scan1','DESMO 14_Scan1', 'DESMO 16_Scan1',  'DESMO 18_Scan1',  'DESMO 20_Scan1',  'DESMO 22_Scan1')

i = 4
for (i in 1:length(id_tumeurs)){
  path_merge<-paste0("D:/Analyse_Inform/", id_tumeurs[[i]],"/Merge_cell_seg_data.txt")
  csd<-read_cell_seg_data(path_merge)
  csd$`Cell X Position`
  x_tumeur<-filter(csd,str_to_lower(csd$`Tissue Category`) == "tumor")$`Cell X Position`
  y_tumeur<-filter(csd,str_to_lower(csd$`Tissue Category`) == "tumor")$`Cell Y Position`
  
  x_stroma<-filter(csd,str_to_lower(csd$`Tissue Category`) == "stroma")$`Cell X Position`
  y_stroma<-filter(csd,str_to_lower(csd$`Tissue Category`) == "stroma")$`Cell Y Position`
  
  x_other<-filter(csd,str_to_lower(csd$`Tissue Category`) == "other")$`Cell X Position`
  y_other<-filter(csd,str_to_lower(csd$`Tissue Category`) == "other")$`Cell Y Position`
  
  path2<-paste0("D:/Analyse_Inform/",id_tumeurs[i],"/Merge_cell_seg_data_summary.txt")
  csd2<-read_cell_seg_data(path2)
  surface<-sum(csd2$`Tissue Category Area (square microns)`)
  nb_cellules<-length(csd$Path)
  
  resolution = (1525424+0.03306167*surface)
  w=sqrt(resolution)
  h=w
  
  png(paste0("D:/Analyse_Inform/",id_tumeurs[[i]],"/check_segmentation",id_tumeurs[[i]],"_Zone_Tumorale", ".png"),width=w,height=h)
  par(bg = 'black', fg = 'white',col.main="white")
  c<-plot(x_other,y_other,main='Zone tumorale',cex=0.01,col='dodgerblue1', asp = 1,
          xlim=c(min(x_tumeur, x_stroma, x_other),max(x_tumeur, x_stroma, x_other)),
          ylim=c(min(y_tumeur, y_stroma, y_other),max(y_tumeur, y_stroma, y_other)))
  points(x_stroma,y_stroma,col='cornsilk',cex=0.01) 
  points(x_tumeur,y_tumeur,col='red1',cex=0.01) 
  dev.off()
}

