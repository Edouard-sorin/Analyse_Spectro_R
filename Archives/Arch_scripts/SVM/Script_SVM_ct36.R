rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(tidyr) # Utilisation de la fonction pivot longer
library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider


library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/Fonction_SVM.R") 

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/data_SPIR_Ed.Rdata")

# Importation du jeu de donnee de Hoarau

load("Sauvegardes_objet_R.data/SPIR_Ho.Rdata")

# Importation du jeu de donnee de Pothin

load("Sauvegardes_objet_R.data/SPIR_Po.Rdata")

# ------------------------------ Support Vector Machine en kernel linear à Ct<36 ------------------------------
# I) Parametres SVM Ct<36 #### 

nb.simu <- 10  # Minimu 1000 simu
rep.max <- 6  # nombre de repetition SPIR sur les feuilles , maximum 6


Tirage <- split(data_SPIR_Ed[,-1], data_SPIR_Ed$code_ech_feuille, drop = T) # drop = T pour enlever les tiroirs vides !!

## I.a) Calcul parallele ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071)
sfLibrary(caret)
sfExport("fct_svm","Tirage","rep.max","nb.simu") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

#res.svm.36 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 36 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

res.svm.36 <- sfClusterApplyLB(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 36 , list.feuilles= Tirage ) # Le plus rapide # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'à mis l'operation 

## I.b) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu_3cpu_sfClusterApplyLB.Rdata")

intermed.36 <- as.data.frame(do.call(rbind, res.svm.36))  # permet de basculer de la liste à la data.frame pour le resultat issu de sfClusterApplyLB

data_global.36 <- pivot_longer((intermed.36), cols = 1:3, names_to = "critere", values_to = "valeurs")

data_global.36$nb.rep <- rep(1:rep.max, each = (rep.max*nb.simu/2))


#save(res.svm.36,data_global.36,  file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu_3cpu_Sapply.Rdata")

# save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_10simu_para") 

# write.table(x =data_global.36 , file = "Donnees/SVM_ct36_6rep_100simu.csv" , sep = ';')


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

## I.d) Graph calculé avec IC ####

tab.36 <-  aggregate(valeurs ~ nb.rep + critere, data_global.36, mean) 
names(tab.36)[3] <- "moyenne"
tab.36$et <- aggregate(valeurs ~ nb.rep + critere, data_global.36, sd)$valeurs
tab.36$mediane <- aggregate(valeurs ~ nb.rep + critere, data_global.36, median )$valeurs
tab.36$nb <-  aggregate(valeurs ~ nb.rep + critere, data_global.36, length )$valeurs

tab.36$IC <- qt(0.975, tab.36$nb -1 )*tab.36$et/sqrt(tab.36$nb)

ggplot(data = tab.36) +
  aes(x = nb.rep, y = moyenne, color = critere, group = critere)+
  geom_errorbar(aes(ymin = moyenne - IC , ymax = moyenne + IC),width=0.05, lwd = 1.1 )+
  scale_x_continuous(breaks=seq(1:6)) +  
  geom_line(lwd = 1.1) +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))




rm (fct_svm,Tirage,nb.simu,rep.max,IC)


# II) Prediction ct<36 ####

# qd pred avec data_hoarau

#data_long_Ed <- data_long_Ho
#trueP <- truep_Ho

seuil.ct <- 36
data_long_Ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(data_long_Ed[[paste0("qPCR_", seuil.ct)]] ))
feuille_spir_36 <- aggregate(data_long_Ed[c(paste0("qPCR_", seuil.ct), "reflectance")] , data_long_Ed[, c("code_ech_arbre", "code_ech_feuille", "lambda")], mean  )

data_pred_36 <- pivot_wider(feuille_spir_36, names_from = "lambda", values_from = "reflectance", names_prefix = "X")
mod.complet.36 <- svm(y = data_pred_36[[paste0("qPCR_", seuil.ct)]]               # ici on prend les seuil donc sois 36 sois 36 , avec "paste0" colle sans separateurs
                      , x = data_pred_36[, grep("^X", names(data_pred_36))]      #  ici on prend ttes les longueurs d'ondes
                      , type = 'C-classification'
                      , kernel = 'linear'
) 

## II.a) Prediction a partir de la base d'apprentissage ####

# Prédiction sur les feuilles de la base d'apprentissage

data_pred_36$svm_pred_36 <- predict(mod.complet.36, decision.values = T) 

# Conversion des resultats de la prediction en numerique

data_pred_36$svm_pred_36 = as.numeric(as.character(data_pred_36$svm_pred_36))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(data_pred_36))
data_pred_36 <- data_pred_36[,c(names(data_pred_36)[-select.lambda], names(data_pred_36)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

svm_pred_arbres = aggregate(svm_pred_36 ~ code_ech_arbre, data = data_pred_36, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.3

# Parametrage pour presentation graphique et la matrice de confusion

svm_pred_arbres$statut_pred <- "Indéterminé"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_36 >= crit.pos  ] <- "Positif"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_36 <= crit.neg  ] <- "Négatif"

svm_pred_arbres$pred_36 <- 0.5
svm_pred_arbres$pred_36[svm_pred_arbres$svm_pred_36 >= crit.pos  ] <- 1
svm_pred_arbres$pred_36[svm_pred_arbres$svm_pred_36 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.32 <- NULL

svm_pred_arbres[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(svm_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion

svm_pred_arbres$pred_36 <- factor(svm_pred_arbres$pred_36 ,levels= c(0,1))

confusion_matrix_36 <- ftable( qPCR_36 ~ pred_36 , data = svm_pred_arbres )


confusion_matrix_36 <- as.matrix(confusion_matrix_36)

TP <- confusion_matrix_36[1,1]

TN <- confusion_matrix_36[2,2]

FN <- confusion_matrix_36[2,1]

FP <- confusion_matrix_36[1,2]

# Calcul des parametres SVM de la matrice de confusion

Accuracy <- TP+TN / TP+TN+FN+FP

Precision <- TP / TP+FP

Sensitivity <- TP / TP+FN

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

# Graphique de la prediction à ct 36 

ggplot(svm_pred_arbres)+
  aes(y = svm_pred_36 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_36)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Résultat de la prédiction hoareau SVM ct36 " ,subtitle = "T = positif au HLB à Ct36 en qPCR", color = "Statut") +
  annotate(geom = "text", x = trueP$seuil.36, y = 0.8, label = "T")

rm (feuille_spir_new,feuille_spir_36,mod.complet.36,select.lambda,data_SPIR_Ed,data_long_Ed,truep_Ho,truep_Po,SPIR_Ho,SPIR_Po,data_long_Ho,data_long_Po, data_pred_36)


## II.b) Prediction a partir de nouvelles donnees ####

data_long_new <- data_long_Po
truep_new <- truep_Po

data_long_new[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(data_long_new[[paste0("qPCR_", seuil.ct)]] ))
feuille_spir_new <- aggregate(data_long_new[c(paste0("qPCR_", seuil.ct), "reflectance")] , data_long_new[, c("code_ech_arbre", "code_ech_feuille", "lambda")], mean  )

data_pred_new <- pivot_wider(feuille_spir_new, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

# Prediction sur des donnees differente de la base d'apprentissage

data_pred_new$svm_pred_new <- predict(mod.complet.36, newdata = data_pred_new[grep("^X", names(data_pred_new))],  decision.values = T)

# Conversion les feuilles de la prediction en numerique

data_pred_new$svm_pred_new = as.numeric(as.character(data_pred_new$svm_pred_new))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(data_pred_new))
data_pred_new <- data_pred_new[,c(names(data_pred_new)[-select.lambda], names(data_pred_new)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

svm_pred_arbres = aggregate(svm_pred_new ~ code_ech_arbre, data = data_pred_new, mean, na.rm = T)


# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.3

# Parametrage pour presentation graphique et la matrice de confusion

svm_pred_arbres$statut_pred <- "Indéterminé"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_new >= crit.pos  ] <- "Positif"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_new <= crit.neg  ] <- "Négatif"

svm_pred_arbres$pred_new <- 0.5
svm_pred_arbres$pred_new[svm_pred_arbres$svm_pred_new >= crit.pos  ] <- 1
svm_pred_arbres$pred_new[svm_pred_arbres$svm_pred_new <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.32 <- NULL

svm_pred_arbres[paste("qPCR", seuil.ct , sep = "_")] <- lapply(truep_new, function(x)
  as.numeric(svm_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion

svm_pred_arbres$pred_new <- factor(svm_pred_arbres$pred_new ,levels= c(0,1))

confusion_matrix_36 <- ftable( qPCR_36 ~ pred_new , data = svm_pred_arbres )


confusion_matrix_36 <- as.matrix(confusion_matrix_36)

VP <- confusion_matrix_36[1,1]

VN <- confusion_matrix_36[2,2]

FN <- confusion_matrix_36[2,1]

FP <- confusion_matrix_36[1,2]

# Calcul des parametres SVM de la matrice de confusion

Accuracy <- VP+VN / VP+VN+FN+FP

Precision <- VP / VP+FP

Sensitivity <- VP / VP+FN

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("Vrai Négatif (VN)","Faux Négatif (FN)(Erreur type 2)","Faux Positif (FP)(Erreur type 1)","Vrai Positif (VP)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Négatif confirmé","Positif confirmé")
rownames(Confusion_matrix) <- c("Négatif prédit","Positif prédit")

print(Confusion_matrix)

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif


ggplot(svm_pred_arbres)+
  aes(y = svm_pred_new , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_new)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Résultat de la prédiction hoareau SVM ct32 " ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +
  annotate(geom = "text", x = truep_new$seuil.36, y = 0.8, label = "T")

rm (feuille_spir_new,feuille_spir_36,mod.complet.36,select.lambda,data_SPIR_Ed,data_long_Ed,truep_Ho,truep_Po,SPIR_Ho,SPIR_Po,data_long_Ho,data_long_Po, data_pred_36)
