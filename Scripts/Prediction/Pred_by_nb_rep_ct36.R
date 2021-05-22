rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(pls) # PLS

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider
library(caret) # fonction confusionMatrix

library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Prediction/Fct_SVM.R")   # pour calculer avec fonction svm

#source(file = "Scripts/Prediction/Fct_PLS.R")  # pour calculer avec fonction pls

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# I) Parametres SVM Ct<36 ####

nb.simu <- 100  # Minimu 1000 simu
rep.max <- 6  # nombre de repetition SPIR sur les feuilles , maximum 6


Tirage <- split(data_SPIR_Ed[,-c(1:2)], data_SPIR_Ed$code_ech_feuille, drop = T) # drop = T pour enlever les tiroirs vides !!


## I.a) Calcul parallele ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071) # pour fct_svm
sfLibrary(pls) # pour fct_pls
sfLibrary(caret)
sfExport("fct_svm","Tirage","rep.max","nb.simu") # les elements exterieur a la fonction 
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.36 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 36 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.svm.36 <- sfClusterApplyLB(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 36 , list.feuilles= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'a mis l'operation 


## I.b) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu_3cpu_sfClusterApplyLB.Rdata")

intermed.36 <- as.data.frame(do.call(rbind, res.svm.36))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

data_global.36 <- pivot_longer((intermed.36), cols = 1:3, names_to = "critere", values_to = "valeurs")

data_global.36$nb.rep <- rep(1:rep.max, each = (rep.max*nb.simu/2))


#save(res.svm.36,data_global.36,  file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu_18_02.Rdata")

#save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_10simu_para") 

#write.table(x =data_global.36 , file = "Donnees/SVM_ct36_6rep_100simu_18_02.csv" , sep = ';')

## I.c) Graph calculé avec stat_summary ####

ggplot(data = data_global.36) +
  aes(x = nb.rep, y = valeurs, color = critere, group = critere)+
  stat_summary(geom = "pointrange", fun.data = function(x) mean_se(x, mult = qt(0.975, length(x) - 1))) +
  stat_summary(geom = "line", fun = mean) +
  scale_x_continuous(breaks=seq(1:6)) +  
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


#ggsave(filename = "Graphiques/SVM_ct36_6rep_100simu_18_02.png", plot = last_plot() ,width = 20 ,height = 15 )

# ## I.d) Graph calculé avec IC ####
# 
# tab.36 <-  aggregate(valeurs ~ nb.rep + critere, data_global.36, mean) 
# names(tab.36)[3] <- "moyenne"
# tab.36$et <- aggregate(valeurs ~ nb.rep + critere, data_global.36, sd)$valeurs
# tab.36$mediane <- aggregate(valeurs ~ nb.rep + critere, data_global.36, median )$valeurs
# tab.36$nb <-  aggregate(valeurs ~ nb.rep + critere, data_global.36, length )$valeurs
# 
# tab.36$IC <- qt(0.975, tab.36$nb -1 )*tab.36$et/sqrt(tab.36$nb)
# 
# ggplot(data = tab.36) +
#   aes(x = nb.rep, y = moyenne, color = critere, group = critere)+
#   geom_errorbar(aes(ymin = moyenne - IC , ymax = moyenne + IC),width=0.05, lwd = 1.1 )+
#   scale_x_continuous(breaks=seq(1:6)) +  
#   geom_line(lwd = 1.1) +
#   labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
#   dark_theme_gray() +
#   theme(legend.position = "right")+
#   scale_colour_viridis_d() +
#   theme(panel.grid.major.y = element_line(colour = "grey20"))
# 
# 
# 
# 
# rm (fct_svm,Tirage,nb.simu,rep.max,IC)
# 
