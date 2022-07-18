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

#####
#chargement des donnees (Melpredict)
datas<-read.table("Melpred-cd8.txt", sep = ",", header = TRUE)
id_tumeurs <- na.omit(datas$id_tumeurs)
path <- na.omit(datas$path)
marqueurs <- na.omit(datas$marqueurs)
opal <- na.omit(datas$opal)
zone_marqueurs <- na.omit(datas$zone_marqueurs)
dossier_data <- na.omit(datas$dossier_data)
phenotypes <- parse_phenotypes("SOX10+", "SOX10+Ki67+","SOX10+ZEB1+","CD8+",
                               "CD8+PD1+","CD8+Ki67+","CD4+",
                               "CD4+PD1+","CD4+Ki67+", "Total Cells")


for(i in 42:length(id_tumeurs)){
  # Read the consolidated data file
  csd_path = paste0("D:melpredict_cd8/", id_tumeurs[i],"/Consolidated_data.txt")
  csd = read_cell_seg_data(csd_path, col_select="phenoptrReports")
  
  # Make a table summarizing the number of fields per slide
  summary_table =   summarize(group_by(csd, `Slide ID`),`Number of fields`=n_distinct(`Annotation ID`))
  
  tissue_categories = c("Tumor", "tumor", "TUMOR")
  
  # Column to aggregate by
  .by = "Slide ID"
  
  # Count phenotypes per tissue category
  counts = count_phenotypes(csd, phenotypes, tissue_categories, .by=.by)
  percents = counts_to_percents(counts)
  
  expression_params = NULL
  
  # Summarize nearest neighbor distances
  nearest_detail_path = file.path(
    paste0("D:melpredict_cd8/", id_tumeurs[i],"/nearest_neighbors.txt"))
  nearest_neighbors = nearest_neighbor_summary(
    csd, phenotypes, tissue_categories, nearest_detail_path, .by=.by,
    extra_cols=expression_params)
  
  # Summary of cells within a specific distance
  radii = 15
  count_detail_path = file.path(
    paste0("D:melpredict_cd8/", id_tumeurs[i],"/count_within.txt"))
  count_within = count_within_summary(
    csd, radii, phenotypes, tissue_categories,
    count_detail_path, .by=.by, extra_cols=expression_params)
  
  # Write it all out to an Excel workbook
  wb = createWorkbook()
  write_summary_sheet(wb, summary_table)
  write_counts_sheet(wb, counts)
  write_percents_sheet(wb, percents)
  write_nearest_neighbor_summary_sheet(wb, nearest_neighbors)
  write_count_within_sheet(wb, count_within)
  
  workbook_path = file.path(
    paste0("D:melpredict_cd8/", id_tumeurs[i],"/Results.xlsx"))
  if (file.exists(workbook_path)) file.remove(workbook_path)
  saveWorkbook(wb, workbook_path)
  
  # Write summary charts
  #charts_path = file.path(
  #  paste0(path, dossier_data, id_tumeurs[i]),"/Charts.docx")
  #if (file.exists(charts_path)) file.remove(charts_path)
  #write_summary_charts(workbook_path, charts_path, .by=.by)
  
  # Save session info
  info_path = file.path(
    paste0("D:melpredict_cd8/",id_tumeurs[i],"/session_info.txt"))
  write_session_info(info_path)
}

#scatterplot 3D pour la distance entre deux phenotypes

d = 1
pp1 = 1
pp2 = 1

for(i in 1:3){
  distances = readWorkbook(paste0("D:melpredict_cd8/", id_tumeurs[i], "/Results.xlsx"), "Nearest Neighbors")
  proportion = readWorkbook(paste0("D:melpredict_cd8/", id_tumeurs[i], "/Results.xlsx"), "Cell Percents")
  distances<-as.data.frame(cbind(distances[2], distances$X3, distances$X4, distances$X6))
  colnames(distances) = c("type", "p1", "p2", "mean")
  distances = distances[distances$type=="Tumor"]
  colnames(proportions) = c(proportion[1,1], proportion[1,2], proportion[1,3], proportion[1,4], proportion[1,5], proportion[1,6], proportion[1,7], proportion[1,8], proportion[1,9], proportion[1,10], proportion[1,11], proportion[1,12], proportion[1,13])
  

  print(phenotypes)
  phenotype1 = readLines("Phenotype 1 : ")
  phenotype2 = readLines("Phenotype 2 : ")
  
  d = c(d, distances$mean[(distances$p1==phenotype1 & distances$p2==phenotype2) | (distances$p2==phenotype1 & distances$p1==phenotype2)])
  pp1 = c(pp1, c(proportion[2, phenotype1], proportion[3, phenotype1]))
  pp2 = c(pp2, c(proportion[2, phenotype2], proportion[3, phenotype2]))
  #là on a d, pp1 et pp2
}

#i = 42