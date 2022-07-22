library(openxlsx)
library(cutpointr)
library(dplyr)

densites = readWorkbook("C:/Users/aliar/Desktop/INSA/Stage/Scatterplot 3D/decompte_phenotypes_cd8.xlsx", "Cell Density (tumor)")

densitecd8 = subset(densites, select = CD8.)
densitecd8$pd1 = "no"
colnames(densitecd8) = c("densites", "pd1")
densitecd8pd1 = subset(densites, select = CD8.PD1.)
densitecd8pd1$pd1 = "yes"
colnames(densitecd8pd1) = c("densites", "pd1")

tableau = bind_rows(densitecd8, densitecd8pd1)
tableau$densites = as.numeric(tableau$densites)
tableau$pd1 = as.factor(tableau$pd1)

cp = cutpointr(tableau, densites, pd1, direction = ">=", method=maximize_metric(), metric = sum_sens_spec())
