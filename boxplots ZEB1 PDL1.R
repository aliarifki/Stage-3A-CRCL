library(phenoptr)
library(phenoptrReports)

#Importation
a = read_cell_seg_data("D:Analyse_Inform/CONT 1_Scan1/Merge_cell_seg_data.txt")
head(a)
head(a$`Tissue Category`)
print(a$`Tissue Category`)
a = a[a$`Tissue Category`=="Tumor", ]
#Noms des colonnes
names(a)
#Celles  qu'on veut c'est intensité et coordonnées

#Noms des dossiers
id_tumeurs_CONT<-c('CONT 1_Scan1', 'CONT 2_Scan1', 'CONT 3_Scan1', 'CONT 4 BIS_Scan1', 'CONT 5_Scan1', 'CONT 6_Scan1', 'CONT 9_Scan1', 'CONT 11_Scan1', 'CONT 12_Scan1', 'CONT 13_Scan1', 'CONT 14_Scan1', 'CONT 15_Scan1', 'CONT 16_Scan1', 'CONT 17_Scan1', 'CONT 18_Scan1', 'CONT 19_Scan1', 'CONT 20_Scan1', 'CONT 21_Scan1', 'CONT 22_Scan1')
id_tumeurs_DESMO<-c('DESMO 1_Scan1', 'DESMO 4_Scan1', 'DESMO 5_Scan1',  'DESMO 8_Scan1', 'DESMO 9_Scan1', 'DESMO 11_Scan1', 'DESMO 12_Scan1', 'DESMO 13 BIS_Scan1','DESMO 14_Scan1', 'DESMO 16_Scan1',  'DESMO 18_Scan1',  'DESMO 20_Scan1',  'DESMO 22_Scan1')

par(mfrow=c(1,2))

#Noms des dossiers
id_tumeurs_CONT<-c('CONT 1_Scan1', 'CONT 2_Scan1', 'CONT 3_Scan1', 'CONT 4 BIS_Scan1', 'CONT 5_Scan1', 'CONT 6_Scan1', 'CONT 9_Scan1', 'CONT 11_Scan1', 'CONT 12_Scan1', 'CONT 13_Scan1', 'CONT 14_Scan1', 'CONT 15_Scan1', 'CONT 16_Scan1', 'CONT 17_Scan1', 'CONT 18_Scan1', 'CONT 19_Scan1', 'CONT 20_Scan1', 'CONT 21_Scan1', 'CONT 22_Scan1')

#Création des listes
PDL1 = c()
ZEB1 = c()

#Boucle Controle
for(i in 1:length(id_tumeurs_CONT)){
  a = read_cell_seg_data(paste0("E:Analyse_Inform/", id_tumeurs_CONT[i], "/Merge_cell_seg_data.txt"))
  a = a[a$`Tissue Category`=="Tumor", ]
  PDL1 = c(PDL1, mean(a$`Entire Cell Opal 620 Mean`))
  ZEB1 = c(ZEB1, mean(a$`Entire Cell Opal 570 Mean`))
}


boxplot(ZEB1,PDL1, names=c("ZEB1", "PDL1"), main = "Contrôle")


PDL1b = 1
ZEB1b = 1
#Boucle Desmoplastique
for(i in 1:length(id_tumeurs_DESMO)){
  a = read_cell_seg_data(paste0("D:Analyse_Inform/", id_tumeurs_DESMO[i], "/Merge_cell_seg_data.txt"))
  a = a[a$`Tissue Category`=="Tumor", ]
  PDL1b = c(PDL1b, mean(a$`Entire Cell Opal 620 Mean`))
  ZEB1b = c(ZEB1b, mean(a$`Entire Cell Opal 570 Mean`))
}
PDL1b = PDL1b[-1]
ZEB1b = ZEB1b[-1]
boxplot(ZEB1b,PDL1b, names=c("ZEB1", "PDL1"), main = "Test")

boxplot(ZEB1, ZEB1b, names = c("Contrôle", "Test"), main = "ZEB1")

#test de corrélation
shapiro.test(ZEB1b)
shapiro.test(ZEB1)
t.test(ZEB1, ZEB1b)

boxplot(PDL1, PDL1b, names = c("Contrôle", "Test"), main = "PDL1")
#test de corrélation
shapiro.test(PDL1)
shapiro.test(PDL1b)
t.test(PDL1, PDL1b)

png("PDL1 en fonction de ZEB1 contrôle.png")
plot(ZEB1, PDL1) 
points(ZEB1b, PDL1b, col=2)
legend( "topleft", legend= c("Contrôle", "Test"), col=c(1, 2), pch=c(1,1))
dev.off()

#Controle
lmc = lm(PDL1~ZEB1)
summary(lmc)
anova(lmc)

#Hypothèses du modèle
#H1
plot(ZEB1,PDL1)
abline(lmc)
legend("topleft", legend = "Modèle", lty=1)

#H2
plot(lmc, 5)
res=residuals(lmc)
grubbs.test(res)
cooks.distance(lmc)

#H3
library(car)
durbinWatsonTest(lmc)

#H4
mean(lmc$residuals)

#H5
library(AER)
library(stats)
dispersiontest(glm(PDL1~ZEB1))
glm(PDL1~ZEB1)

#H6
shapiro.test(lmc$residuals) 
hist(PDL1)

#Desmoplastique
lmd = lm(PDL1b~ZEB1b)
summary(lmd)
anova(lmd)

#Hypothèses du modèle
#H1
plot(ZEB1b,PDL1b)
abline(lmd)

#H2
plot(lmd, 5)
res=residuals(lmc)
grubbs.test(res)
cooks.distance(lmd)

#H3
library(car)
durbinWatsonTest(lmd)

#H4
mean(lmd$residuals)

#H5
library(AER)
library(stats)
dispersiontest(glm(PDL1~ZEB1))
glm(PDL1~ZEB1)

#H6
shapiro.test(lmd$residuals) 
hist(PDL1b)

#Graphique avec toutes les cellules controles
plot(ZEB1, PDL1)
