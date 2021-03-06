
rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(pls) # package pls
library(randomForest)# package RF
library(e1071) # package SVM

library(tidyr) # pivot longer & pivot wider

library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Prediction/Ct32/Fct_ConfusionMatrix_ct32.R") 

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Parametres MachinLearning Ct<32 ####

nb.simu <- 10  # Minimu 1000 simu
rep.max <- 6  # nombre de repetition SPIR sur les feuilles , maximum 6

#Time difference of 18.0745 mins pour 10 simu a 6 rep


Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T) # drop = T pour enlever les tiroirs vides !!


# Calcul parralelle

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)  # la library des packages utilisés
sfLibrary(pls)
sfLibrary(randomForest)
sfLibrary(e1071)
sfLibrary(caret)
sfExport("fct_ConfusionMatrix","Tirage","rep.max","nb.simu","trueP") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.32 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.ML.32 <- sfClusterApplyLB(rep(rep.max, each = nb.simu), fct_ConfusionMatrix, seuil.ct = 32 , list.feuilles= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'à mis l'operation 

intermed.32 <- as.data.frame(do.call(rbind,res.ML.32))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

ML_global.32 <- pivot_longer((intermed.32), cols = 1:9, names_to = "critere", values_to = "valeurs")

# On moyenne tout les parametres pour chaque type de machin learning

ML.32 <-  aggregate(valeurs ~ critere, ML_global.32, mean)

names(ML.32)[2] <- "Moyenne"

# Meme chose avec ecart type et mediane

ML.32$et <- aggregate(valeurs ~ critere, ML_global.32, sd)$valeurs
ML.32
