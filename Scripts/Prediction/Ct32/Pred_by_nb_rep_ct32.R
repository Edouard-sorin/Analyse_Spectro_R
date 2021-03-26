rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider
library(caret) # fonction confusionMatrix

library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Prediction/Fct_SVM.R") 


# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# I) Parametres SVM Ct<32 ####

nb.simu <- 5  # Minimu 1000 simu
rep.max <- 3  # nombre de repetition SPIR sur les feuilles , maximum 6


Tirage <- split(data_SPIR_Ed[,-1], data_SPIR_Ed$code_ech_feuille, drop = T) # drop = T pour enlever les tiroirs vides !!


## I.a) Calcul parallele ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071)
sfLibrary(caret)
sfExport("fct_svm","Tirage","rep.max","nb.simu") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.32 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.svm.32 <- sfClusterApplyLB(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'à mis l'operation 


## I.b) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct32_6rep_100simu_3cpu_sfClusterApplyLB.Rdata")

intermed.32 <- as.data.frame(do.call(rbind, res.svm.32))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

data_global.32 <- pivot_longer((intermed.32), cols = 1:3, names_to = "critere", values_to = "valeurs")

data_global.32$nb.rep <- rep(1:rep.max, each = (rep.max*nb.simu/2))


#save(res.svm.32,data_global.32,  file = "Sauvegardes_objet_R.data/SVM_ct32_6rep_100simu_3cpu_Sapply.Rdata")

# save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct32_6rep_10simu_para") 

# write.table(x =data_global.32 , file = "Donnees/SVM_ct32_6rep_100simu.csv" , sep = ';')

## I.c) Graph calculé avec stat_summary ####

ggplot(data = data_global.32) +
  aes(x = nb.rep, y = valeurs, color = critere, group = critere)+
  stat_summary(geom = "pointrange", fun.data = function(x) mean_se(x, mult = qt(0.975, length(x) - 1))) +
  stat_summary(geom = "line", fun = mean) +
  scale_x_continuous(breaks=seq(1:6)) +  
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<32), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


# ggsave(filename = "Graphiques/SVM_ct32_6rep_100simu.png", plot = last_plot() ,width = 900 ,height = 680 )

## I.d) graph calculé avec IC ####

tab.32 <-  aggregate(valeurs ~ nb.rep + critere, data_global.32, mean) 
names(tab.32)[3] <- "moyenne"
tab.32$et <- aggregate(valeurs ~ nb.rep + critere, data_global.32, sd)$valeurs
tab.32$mediane <- aggregate(valeurs ~ nb.rep + critere, data_global.32, median )$valeurs
tab.32$nb <-  aggregate(valeurs ~ nb.rep + critere, data_global.32, length )$valeurs

tab.32$IC <- qt(0.975, tab.32$nb -1 )*tab.32$et/sqrt(tab.32$nb)

ggplot(data = tab.32) +
  aes(x = nb.rep, y = moyenne, color = critere, group = critere)+
  geom_errorbar(aes(ymin = moyenne - IC , ymax = moyenne + IC),width=0.05, lwd = 1.1 )+
  scale_x_continuous(breaks=seq(1:6)) +  
  geom_line(lwd = 1.1) +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<32), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))




rm (fct_svm,Tirage,nb.simu,rep.max,IC)

