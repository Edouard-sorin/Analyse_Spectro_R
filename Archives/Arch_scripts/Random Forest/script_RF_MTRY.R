# ----PACKAGES--------------------------------------------------------------------------------
	
	# Chargement des packages utiles pour notre script

# Usage de l'algorithme random Forest
library(randomForest)# permet de lancer l'algorithme de RF
library(rpart) #Cloisonnement recursif pour la classification = RF
library(rpart.plot) # met a echelle et ajuste automatiquement l'arbre = RF
library(caret) # fonction de conversion pour la confusion des classesMatrix
library(ModelMetrics) # fonction de conversion pour la confusion des classesMatrix

# Mise en forme des donnees 
library(tidyr) # transformation du format des donnees : pivot_longer
library(tidyverse)#faciliter l'installation et le chargement de plusieurs paquets "tidyverse"
library(dplyr) # travail sur les donnees pour avoir des analyses performantes
library(sf) #prise en charge de fonctions simples,standardisation de codage DATA vectorielles
library(foreign) # pour utiliser la fct read.dbf
library(openxlsx)
#----------------------------------------------------------------------------------------------
#  					            IMPORTATIONS / MANIPULATIONS BASE APPRENTISSAGE
# 						            pretraitements et traitements des donnees 
#----------------------------------------------------------------------------------------------

#importation de la couche de polygones _ nomee Data
# Post calibration et validation du model RF cree sur la base d'apprentissage, 
# nous predirons tous les autres polygones (Data) avec le modele construit sur les trainData (modele = poly_rf)

#importation de la couche vecteur Data au format .dbf : permet l'interpolation entre Qgis et R
Data <- read.dbf("donnees/Data_TOT.dbf")
summary(Data, maxsum = 7 )


#extraction de la base d'apprentissage, trainData
trainData <- Data[! is.na (Data$CLASSname),]
summary(trainData, maxsum = 7)

#------------------------------------------------------------------------------------------------
#                                  MODELE RANDOM FOREST 
# -----------------------------------------------------------------------------------------------

#Lancer le modele de prediction  
poly_rf <- randomForest(CLASSname ~ . - Id, data=trainData, proximity=T)# pas optimisé
poly_rf <- randomForest(CLASSname ~ . - Id, data=trainData, ntree= 400, proximity=T, mtry = 10)# optimisé 
print(poly_rf)

# comment selectionner le nombre de variables aleatoires utilisees dans chaque arbre?  (mtry)
# mtry <- tuneRFmedian(trainData[,-1], trainData$CLASSname, ntreeTry = 400, mtryStart = 1,
# 							 stepFactor = 2, improve = 0.2, trace = T, plot = T, doBest = T, proximity = T)
# opt.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]

test.RF <- sapply(1:15, function(x) median(randomForest(CLASSname ~ . - Id, data=trainData, ntree= 400, proximity = T, mtry = x)$err.rate[,"OOB"]))
names(test.RF) <- 1:15

pdf("figures/estimation mtry optimal.pdf")
plot(test.RF, pch = 16, las = 1, ylab = "erreur OOB médiane", xlab = "Valeur de Mtry")
lines(test.RF)
abline(v =11, col = "red", lty = 2)
dev.off()

temps.calcul <- sapply(1:15, function(x) system.time(randomForest(CLASSname ~ . - Id, data=trainData, ntree= 400, proximity = T, mtry = x)$err.rate[,"OOB"]))

pdf("figures/relation temps calcul erreur OOB pour mtry.pdf")
plot(test.RF, temps.calcul["elapsed",], type = "b")
text(test.RF-0.0005, temps.calcul["elapsed",], 1:15)
dev.off()

sortie <- rbind(test.RF, temps.calcul)
write.xlsx(sortie, file = "donnees/sorties RF.xlsx", row.names = T) 
