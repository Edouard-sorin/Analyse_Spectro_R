
rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider


library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation des fonctions utiles

source(file = "Scripts/SVM/Fonction_SVM.R") 

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/data_SPIR_Ed.Rdata")

# Importation du jeu de donnee de Hoarau

load("Sauvegardes_objet_R.data/SPIR_Ho.Rdata")

# Importation du jeu de donnee de Pothin

load("Sauvegardes_objet_R.data/SPIR_Po.Rdata")

# Importation du jeu de donnee Global Pred_svm_32 produit en II)


load("Sauvegardes_objet_R.data/Pred_svm_32.Rdata") # Global


load("Sauvegardes_objet_R.data/Pred_svm_32_Hoarau.Rdata") # Hoarau 


load("Sauvegardes_objet_R.data/Pred_svm_32_Pothin.Rdata") # Pothin

# ------------------------------Support Vector Machine en kernel linear à Ct<32------------------------------
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


# II) Prediction Ct<32 ####

# Base d'apprentissage sur le jeu de donnee global 

# qd pred avec data_hoarau

#data_long_Ed <- data_long_Po
#trueP <- truep_Po

seuil.ct <- 32

data_long_Ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(data_long_Ed[[paste0("qPCR_", seuil.ct)]] ))

feuille_spir_32 <- aggregate(data_long_Ed[c(paste0("qPCR_", seuil.ct), "reflectance")] , data_long_Ed[, c("code_ech_arbre", "code_ech_feuille","code_variete", "lambda")], mean  )

data_pred_32 <- pivot_wider(feuille_spir_32, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

#save(feuille_spir_32, data_pred_32,  file = "Sauvegardes_objet_R.data/Pred_svm_32_Pothin.Rdata")


mod.complet.32 <- svm(y = data_pred_32[[paste0("qPCR_", seuil.ct)]]       # ici on prend les seuil donc sois 32 sois 36 , avec "paste0" colle sans separateurs
                      , x = data_pred_32[, grep("^X", names(data_pred_32))]      #  ici on prend ttes les longueurs d'ondes
                      , type = 'C-classification'
                      , kernel = 'linear'
) 

## II.a) Base d'apprentissage en fonction des varietes ####

#var <- "Citron"
#var <- "Zanzibar"
var <- "Tangor"

# Construction d'un masque logique pour filtrer en fonction des variete 

Masque_var <- data_pred_32[, grep("code_variete", names(data_pred_32))] == var

# Construction de la base d'apprentissage avec la variete choisis

data_pred_32 <- data_pred_32[Masque_var,]

# Modification de trueP ?

#truep_var <- trueP

#truep_var$seuil.32 <- lapply(Citron_32$code_ech_arbre, function(x) factor(truep_var$seuil.32 %in% x ) )


## II.b) Prediction a partir de la base d'apprentissage ####


# Prediction sur les feuilles de la base d'apprentissage   ( Obligatoirement le meme nbr de lignes, Comment avoir une prediction pour n'importe quelle taille de jeu de donne ?)

data_pred_32$svm_pred_32 <- predict(mod.complet.32, decision.values = T)

# Conversion des resultats de la prediction en numerique

data_pred_32$svm_pred_32 = as.numeric(as.character(data_pred_32$svm_pred_32))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(data_pred_32))
data_pred_32 <- data_pred_32[,c(names(data_pred_32)[-select.lambda], names(data_pred_32)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

svm_pred_arbres = aggregate(svm_pred_32 ~ code_ech_arbre, data = data_pred_32, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.3

# Parametrage pour presentation graphique et la matrice de confusion

svm_pred_arbres$statut_pred <- "Indéterminé"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_32 >= crit.pos  ] <- "Positif"
svm_pred_arbres$statut_pred[svm_pred_arbres$svm_pred_32 <= crit.neg  ] <- "Négatif"

svm_pred_arbres$pred_32 <- 0.5
svm_pred_arbres$pred_32[svm_pred_arbres$svm_pred_32 >= crit.pos  ] <- 1
svm_pred_arbres$pred_32[svm_pred_arbres$svm_pred_32 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.36 <- NULL

svm_pred_arbres[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(svm_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

svm_pred_arbres$pred_32 <- factor(svm_pred_arbres$pred_32 ,levels= c(0,1))

confusion_matrix_32 <- ftable( qPCR_32 ~ pred_32 , data = svm_pred_arbres )

confusion_matrix_32 <- as.matrix(confusion_matrix_32)

TP <- confusion_matrix_32[1,1]

TN <- confusion_matrix_32[2,2]

FN <- confusion_matrix_32[2,1]

FP <- confusion_matrix_32[1,2]

# Calcul des parametres SVM de la matrice de confusion

Accuracy <- ((TP+TN) / (TP+TN+FN+FP)*100)

Precision <- ((TP / (TP+FP)*100))

Sensitivity <- ((TP / (TP+FN)*100))

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif


ggplot(svm_pred_arbres)+
  aes(y = svm_pred_32 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_32)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction SVM Global a ct32  " ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +
  annotate(geom = "text", x = trueP$seuil.32, y = 0.8, label = "T")

#ggsave("Graphiques/Prédiction SVM Global ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

rm (feuille_spir_32,mod.complet.32,select.lambda,crit.pos,crit.neg,seuil.ct,data_SPIR_Ed,data_long_Ed,trueP)


## II.c) Prediction a partir de de nouvelles donnees ####

data_long_new <- data_long_Ho
truep_new <- truep_Ho

data_long_new[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(data_long_new[[paste0("qPCR_", seuil.ct)]] ))
feuille_spir_new <- aggregate(data_long_new[c(paste0("qPCR_", seuil.ct), "reflectance")] , data_long_new[, c("code_ech_arbre", "code_ech_feuille", "lambda")], mean  )

data_pred_new <- pivot_wider(feuille_spir_new, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

# Prediction sur des donnees differente de la base d'apprentissage

data_pred_new$svm_pred_new <- predict(mod.complet.32, newdata = data_pred_new[grep("^X", names(data_pred_new))],  decision.values = T)

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

truep_new$seuil.36 <- NULL

svm_pred_arbres[paste("qPCR", seuil.ct, sep = "_")] <- lapply(truep_new, function(x)
  as.numeric(svm_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion

svm_pred_arbres$pred_new <- factor(svm_pred_arbres$pred_new ,levels= c(0,1))

confusion_matrix_32 <- ftable( qPCR_32 ~ pred_new , data = svm_pred_arbres )

confusion_matrix_32 <- as.matrix(confusion_matrix_32)

VP <- confusion_matrix_32[1,1]

VN <- confusion_matrix_32[2,2]

FN <- confusion_matrix_32[2,1]

FP <- confusion_matrix_32[1,2]

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
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction SVM à Ct32" ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +
  annotate(geom = "text", x = truep_new$seuil.32, y = 0.8, label = "T")

#ggsave("Graphiques/Prédiction SVM Hoarau avec base d'apprentissage Pothin à Ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

rm (feuille_spir_new,feuille_spir_32,mod.complet.32,select.lambda,data_SPIR_Ed,data_long_Ed,truep_Ho,truep_Po,SPIR_Ho,SPIR_Po,data_long_Ho,data_long_Po, data_pred_32)


# III) Inflence de la variete sur la prediction ? ####

data_long_Ed_var <- data_long_Ho

data_long_Ed_var$var <- 1
data_long_Ed_var$var[data_long_Ed_var$code_variete == "Tangor"  ] <- 2
data_long_Ed_var$var[data_long_Ed_var$code_variete == "Zanzibar"] <- 3


intermed.var <- aggregate(data_long_Ed_var[c("var", "reflectance")] , data_long_Ed_var[, c("code_variete","code_ech_arbre", "code_ech_feuille", "lambda")], mean  )

data_pred_var <- pivot_wider(intermed.var, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

mod.complet.var <- svm(y = data_pred_var$var               
                      , x = data_pred_var[, grep("^X", names(data_pred_var))]      
                      , type = 'C-classification'
                      , kernel = 'linear'
) 


## III.a) Prediction des variete a partir de la base d'apprentissage ####



#data_pred_var$svm_pred_var <- predict(mod.complet.var, decision.values = T)
data_pred_var$svm_pred_var <- predict(mod.complet.var,data_pred_var[, grep("^X", names(data_pred_var))] , decision.values = T)

data_pred_var$svm_pred_var = as.numeric(as.character(data_pred_var$svm_pred_var))

select.lambda <- grep("^X", names(data_pred_var))
data_pred_var <- data_pred_var[,c(names(data_pred_var)[-select.lambda], names(data_pred_var)[select.lambda] )]


svm_var = aggregate(svm_pred_var ~ code_ech_arbre, data = data_pred_var, mean, na.rm = T)

var.Tang <- 1.5
var.Zanz <- 2.5

svm_var$pred_var <- "Citron"
svm_var$pred_var[svm_var$svm_pred_var >= var.Tang  ] <- "Tangor"
svm_var$pred_var[svm_var$svm_pred_var >= var.Zanz  ] <- "Zanzibar"

# confirmation de la prediction 

conf_var = aggregate(var ~ code_variete + code_ech_arbre, data = data_pred_var, mean, na.rm = T)

conf_var$code_variete <- factor(substr(conf_var$code_variete, 1, 1))

rm (data_long_Ed_var,select.lambda,mod.complet.var,data_long_Ed,trueP,data_SPIR_Ed,intermed.var)

# Graphique de la prediction des varietes a partir de la base d'apprentissage


ggplot(svm_var)+
  aes(y = svm_pred_var , x = code_ech_arbre, color = pred_var) +
  #geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_var)) +
  scale_color_manual(values = c(Citron = "yellow", Tangor = "orange", Zanzibar = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen de la prediction des variete entre 1 et 3", title = "Résultat de la prédiction " ,subtitle = "Citron = 1, Tangor = 2, Zanzibar = 3", color = "Variete predite") +
  
  #facet_wrap(~ statut, scales = "free_x")
  annotate(geom = "text", x = conf_var$code_ech_arbre, y = 3.2, label = conf_var$code_variete)

## III.b) Prediction sur les feuilles à partir de nouvelles donnees ####

data_long_new <- data_long_Po

data_long_new$var <- 1
data_long_new$var[data_long_new$code_variete == "Tangor"  ] <- 2
data_long_new$var[data_long_new$code_variete == "Zanzibar"] <- 3


intermed.new <- aggregate(data_long_new[c("var", "reflectance")] , data_long_new[, c("code_variete","code_ech_arbre", "code_ech_feuille", "lambda")], mean  )

pred.new.data <- pivot_wider(intermed.new, names_from = "lambda", values_from = "reflectance", names_prefix = "X")

# essayer de predir new data !

pred.new.data$svm_pred_var <- predict(mod.complet.var, newdata = pred.new.data[grep("^X", names(pred.new.data))],  decision.values = T)

pred.new.data$svm_pred_var = as.numeric(as.character(pred.new.data$svm_pred_var))

select.lambda <- grep("^X", names(pred.new.data))
pred.new.data <- pred.new.data[,c(names(pred.new.data)[-select.lambda], names(pred.new.data)[select.lambda] )]


svm_var = aggregate(svm_pred_var ~ code_ech_arbre, data = pred.new.data, mean, na.rm = T)

var.Tang <- 1.5
var.Zanz <- 2.5

svm_var$pred_var <- "Citron"
svm_var$pred_var[svm_var$svm_pred_var >= var.Tang  ] <- "Tangor"
svm_var$pred_var[svm_var$svm_pred_var >= var.Zanz  ] <- "Zanzibar"

# confirmation de la prediction 

conf_var = aggregate(var ~ code_variete + code_ech_arbre, data = pred.new.data, mean, na.rm = T)

conf_var$code_variete <- factor(substr(conf_var$code_variete, 1, 1))

rm (data_long_Ed_var,data_long_new,data_long_Ho,data_long_Po,select.lambda,truep_Po,truep_Ho,SPIR_Ho,SPIR_Po,intermed.var,intermed.new)

# Graphique de la prediction des varietes a partir de la base d'apprentissage


ggplot(svm_var)+
  aes(y = svm_pred_var , x = code_ech_arbre, color = pred_var) +
  #geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_var)) +
  scale_color_manual(values = c(Citron = "yellow", Tangor = "orange", Zanzibar = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen de la prediction des variete entre 1 et 3", title = "Résultat de la prédiction " ,subtitle = "Citron = 1, Tangor = 2, Zanzibar = 3", color = "Variete predite") +
  annotate(geom = "text", x = conf_var$code_ech_arbre, y = 3.2, label = conf_var$code_variete)
