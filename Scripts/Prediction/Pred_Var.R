rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # package SVM

library(tidyr) # pivot longer & pivot wider

library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Prediction/Fct_Var.R") 

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/data_SPIR_Ed.Rdata")

# Prediction Var ####

nb.simu <- 5  # Minimu 1000 simu
rep.max <- 6  # nombre de repetition SPIR sur les feuilles , maximum 6


Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071)
sfLibrary(caret)
sfExport("fct_var","Tirage","rep.max","nb.simu") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.32 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.svm.var <- sfClusterApplyLB(rep(rep.max, each = nb.simu), fct_var , list.feuilles= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'à mis l'operation 


## I.b) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct32_6rep_100simu_3cpu_sfClusterApplyLB.Rdata")

intermed.32 <- as.data.frame(do.call(rbind, res.svm.var ))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

data_global.32 <- pivot_longer((intermed.32), cols = 1:3, names_to = "critere", values_to = "valeurs")

var.32 <-  aggregate(valeurs ~ critere, data_global.32, mean) 

names(var.32)[2] <- "moyenne"

var.32$et <- aggregate(valeurs ~ critere, data_global.32, sd)$valeurs

var.32
