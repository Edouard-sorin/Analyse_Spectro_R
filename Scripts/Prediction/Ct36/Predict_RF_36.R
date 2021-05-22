rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####


# Usage de l'algorithme random Forest
library(randomForest)# permet de lancer l'algorithme de RF
library(rpart) #Cloisonnement recursif pour la classification = RF
library(rpart.plot) # met a echelle et ajuste automatiquement l'arbre = RF
library(caret) # fonction de conversion pour la confusion des classesMatrix
library(ModelMetrics) # fonction de conversion pour la confusion des classesMatrix
library(party)  # Plot l'arbre de decision en Random Forest

# Mise en forme des donnees 
library(tidyr) # transformation du format des donnees : pivot_longer
library(tidyverse)#faciliter l'installation et le chargement de plusieurs paquets "tidyverse"
library(dplyr) # travail sur les donnees pour avoir des analyses performantes
library(sf) #prise en charge de fonctions simples,standardisation de codage DATA vectorielles
library(foreign) # pour utiliser la fct read.dbf
library(caTools) # sample.split

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir

# Importation du jeu de donnees

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Preparation du training_rf et test_rf  ####

nb.rep <- 6
seuil.ct <- 36

Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)


maliste <- lapply(Tirage,function(feuille){ 
  list_rf <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
  # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
  code <- grep("^code_ech", names(list_rf))
  qPCR <- grep("^qPCR_", names(list_rf))
  sortie <- cbind.data.frame(unique(list_rf[c(code,qPCR)])
                             , matrix(apply(list_rf[-(which(colnames(list_rf) == "code_variete"):which(colnames(list_rf) == "qPCR_36"))], 2, mean)
                                      , nr = 1, dimnames = list(NULL, names(list_rf)[-(which(colnames(list_rf) == "code_variete"):which(colnames(list_rf) == "qPCR_36"))])))
  sortie
})

test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 

decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
train_rf <- test_ed[decoup,]
test_rf <- test_ed[!decoup,] 


# Prediction a partir du training set et du test set a ct36 ####

#set.seed(1)  # For reproducibility

# Definition of out of bag error (OOB)

# l’échantillon bootstrap, est les données choisies pour être «dans le sac» par échantillonnage avec remplacement. L’ensemble hors sac est toutes les données qui ne sont pas choisies dans le processus d’échantillonnage.L’erreur hors sac (OOB) est l’erreur moyenne pour chaque calcul à l’aide de prédictions provenant des arbres qui ne contiennent pas dans leur échantillon bootstrap respectif

# pour connaitre les paramètre d'erreus en fonction du nombre d'arbres effectués


model_rf_36 <- randomForest(
  x = train_rf[, grep("^X", names(train_rf))], 
  y = train_rf[[paste0("qPCR_", seuil.ct)]],
  #xtest = test_rf[,Pred],
  ntree = 1000
)


print(model_rf_36)
plot(model_rf_36)

rf_pred <- test_rf  # Prediction sur le jeu de donne test qui a le meme level pour code arbre et predi donc sur tout les abres

# Prediction sur les feuilles de la base d'apprentissage

rf_pred$rf_pred_36 <- predict(model_rf_36,newdata = test_rf, decision.values = T)

# Conversion des resultats de la prediction en numerique

rf_pred$rf_pred_36 = as.numeric(as.character(rf_pred$rf_pred_36))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(rf_pred))
rf_pred <- rf_pred[,c(names(rf_pred)[-select.lambda], names(rf_pred)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

rf_pred_arbres_36 = aggregate(rf_pred_36 ~ code_ech_arbre, data = rf_pred, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.35

# Parametrage pour presentation graphique et la matrice de confusion

rf_pred_arbres_36$statut_pred <- "Indéterminé"
rf_pred_arbres_36$statut_pred[rf_pred_arbres_36$rf_pred_36 >= crit.pos  ] <- "Positif"
rf_pred_arbres_36$statut_pred[rf_pred_arbres_36$rf_pred_36 <= crit.neg  ] <- "Négatif"

rf_pred_arbres_36$crit <- 0.5
rf_pred_arbres_36$crit[rf_pred_arbres_36$rf_pred_36 >= crit.pos  ] <- 1
rf_pred_arbres_36$crit[rf_pred_arbres_36$rf_pred_36 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.32 <- NULL

rf_pred_arbres_36[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(rf_pred_arbres_36$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

rf_pred_arbres_36$crit <- factor(rf_pred_arbres_36$crit ,levels= c(0,1))

rf_confusion_matrix_36 <- ftable( qPCR_36 ~ crit , data = rf_pred_arbres_36 )

rf_confusion_matrix_36 <- as.matrix(rf_confusion_matrix_36)

TP <- rf_confusion_matrix_36[1,1]

TN <- rf_confusion_matrix_36[2,2]

FN <- rf_confusion_matrix_36[2,1]

FP <- rf_confusion_matrix_36[1,2]

#rf_confusion_matrix_36

# Calcul des parametres rf de la matrice de confusion

Accuracy_rf36 <- ((TP+TN) / (TP+TN+FN+FP)*100)

Precision_rf36 <- ((TP / (TP+FP)*100))

Sensitivity_rf36 <- ((TP / (TP+FN)*100))

Parametre_rf_36 <- rbind(Accuracy_rf36,Precision_rf36,Sensitivity_rf36)

names(Parametre_rf_36) <- c("Accuracy","Precision","Sensitivity")

rf_pred_arbres_36$Vrai = "FAUX"
rf_pred_arbres_36$Vrai[rf_pred_arbres_36$code_ech_arbre %in% trueP$seuil.36] = "VRAI"

#save(rf_confusion_matrix_36,Accuracy_rf36,Precision_rf36,Sensitivity_rf36,rf_pred_arbres_36,Parametre_rf_36,  file = "Sauvegardes_objet_R.data/Pred_RF/parameter_RF_ct36_1000trees.Rdata")

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif

rf_pred <- ggplot(rf_pred_arbres_36)+
  aes(y = rf_pred_36 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T7", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point(aes(shape = Vrai, size = Vrai), fill ="blue") + 
  scale_size_manual(values =c(2.1,2.9))+
  scale_shape_manual(values =c(4,21))+
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = rf_pred_36)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction rf Global a ct36" ,subtitle = "T = positif au HLB à Ct36 en qPCR", color = "Statut") +dark_theme_gray()

rf_pred

rf_confusion_matrix_36

Parametre_rf_36

ggsave("Graphiques/Graph_RF/Prediction_rf_Global_ct36.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)
