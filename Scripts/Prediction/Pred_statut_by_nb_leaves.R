rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(pls) # PLS

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider


library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Prediction/Fct_Feuilles.R") # fonction feuille avec SVM


#source(file = "Scripts/Prediction/Fct_Feuilles_PLS.R") #  fonction feuille avec PLS

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# I) Parametres SVM Ct<36 ####

nb.simu <- 100  # Minimu 1000 simu  Time difference of 15.39098 hours pour 1000 simu
nb.feuille <- 10  # nombre de feuilles echantillonees par arbre , maximum 10
seuil.ct <- 36

intermed.arbre <- aggregate(reflectance ~ code_ech_feuille + code_ech_arbre + qPCR_32 + qPCR_36 + lambda, data =  data_long_Ed, mean)

mean.arbre <- pivot_wider(intermed.arbre, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

Tirage <- split(mean.arbre, mean.arbre$code_ech_arbre, drop = T) # drop = T pour enlever les tiroirs vide


## I.a) Calcul parallele ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071)
sfLibrary(pls)
sfLibrary(caret)
sfExport("fct_feuille","Tirage","nb.feuille","nb.simu") # les elements exterieur a la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.36 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 36 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.svm.36 <- sfClusterApplyLB(rep(1:nb.feuille, each = nb.simu), fct_feuille, seuil.ct = 36 ,list.arbres= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'a mis l'operation 


## I.b) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu_3cpu_sfClusterApplyLB.Rdata")

intermed.36 <- as.data.frame(do.call(rbind, res.svm.36))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

data_global.36 <- pivot_longer((intermed.36), cols = 1:3, names_to = "critere", values_to = "valeurs")

data_global.36$nb.rep <- rep(1:nb.feuille, each = (nb.simu*3))

#save(res.svm.36,data_global.36,  file = "Sauvegardes_objet_R.data/Pred_nb_feuille_ct36_100_simu_19_04.Rdata")

# save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_10simu_para") 

# write.table(x =data_global.36 , file = "Donnees/SVM_ct36_6rep_100simu.csv" , sep = ';')

## I.c) Graph calculé avec stat_summary ####

ggplot(data = data_global.36) +
  aes(x = nb.rep, y = valeurs, color = critere, group = critere)+
  stat_summary(geom = "pointrange", fun.data = function(x) mean_se(x, mult = qt(0.975, length(x) - 1))) +
  stat_summary(geom = "line", fun = mean) +
  scale_x_continuous(breaks=seq(1:10)) +  
  labs(x = "Nombre de feuilles échantillonés sur chaque arbre ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de feuilles échantillonés sur chaque arbre ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 1000 SVM", color = "Paramètres de robustesse en SVM ") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


#ggsave("Graphiques/Graph_param_SVM/Pred_nb_feuille_ct36_100_simu_19_04.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

#save(data_global.36, file = "Sauvegardes_objet_R.data/Pred_SVM/Param_SVM_nb_feuille_ct36_1000_simu")

#write.table(x = data_global.36  , file = "Donnees/Simu_SVM//Param_SVM_nb_feuille_ct36_1000_simu.csv" , sep = ';')
