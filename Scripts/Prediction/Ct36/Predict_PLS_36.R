rm(list=ls())  # nettoyage des listes de l'environnement de travail

T1 <- Sys.time()

# Library ####


# Usage de l'algorithme partial least square

library(pls)
library(plsdepot) # Partial Least Square

# Mise en forme des donnees 
library(tidyr) # transformation du format des donnees : pivot_longer
library(tidyverse)#faciliter l'installation et le chargement de plusieurs paquets "tidyverse"

library(caTools) # sample.split

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir

# Importation du jeu de donnees

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Preparation du training_pls et test_pls  ####

nb.rep <- 6
seuil.ct <- 36

Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)


maliste <- lapply(Tirage,function(feuille){ 
  list_pls <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
  # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
  code <- grep("^code_ech", names(list_pls))
  qPCR <- grep("^qPCR_", names(list_pls))
  sortie <- cbind.data.frame(unique(list_pls[c(code,qPCR)])
                             , matrix(apply(list_pls[-(which(colnames(list_pls) == "code_variete"):which(colnames(list_pls) == "qPCR_36"))], 2, mean)
                                      , nr = 1, dimnames = list(NULL, names(list_pls)[-(which(colnames(list_pls) == "code_variete"):which(colnames(list_pls) == "qPCR_36"))])))
  sortie
})

test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 

test_ed <- test_ed[-(which(colnames(test_ed) == "qPCR_32"))]

test_ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(test_ed[[paste0("qPCR_", seuil.ct)]] ))

decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
train_pls <- test_ed[decoup,]
test_pls <- test_ed[!decoup,] 


# Prediction a partir du training set et du test set a ct36 ####


model_pls_36 <-  plsr(
  
  train_pls[[paste0("qPCR_", seuil.ct)]] ~ . ,
  data = train_pls[, grep("^X", names(train_pls))], 
  scale = TRUE, 
  validation = "CV"
  
)

# Prediction

plot( model_pls_36,main="Test Dataset", xlab="observed", ylab="PLS Predicted")

print(model_pls_36)

pls_pred <- test_pls

# Prediction sur les feuilles de la base d'apprentissage

pls_pred$pls_pred_36 <- predict(model_pls_36, newdata = test_pls[, grep("^X", names(test_pls))], decision.values = T, ncomp=100)

# Conversion des resultats de la prediction en numerique

pls_pred$pls_pred_36 = as.numeric(as.character(pls_pred$pls_pred_36))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(pls_pred))
pls_pred <- pls_pred[,c(names(pls_pred)[-select.lambda], names(pls_pred)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

pls_pred_arbres_36 = aggregate(pls_pred_36 ~ code_ech_arbre, data = pls_pred, mean, na.rm = T)

pls_pred_arbres_36$pls_pred_36[pls_pred_arbres_36$pls_pred_36 < 0] = 0   # Pour enlever les valeurs négatives

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.35

# Parametrage pour presentation graphique et la matrice de confusion

pls_pred_arbres_36$statut_pred <- "Indéterminé"
pls_pred_arbres_36$statut_pred[pls_pred_arbres_36$pls_pred_36 >= crit.pos  ] <- "Positif"
pls_pred_arbres_36$statut_pred[pls_pred_arbres_36$pls_pred_36 <= crit.neg  ] <- "Négatif"

pls_pred_arbres_36$crit <- 0.5
pls_pred_arbres_36$crit[pls_pred_arbres_36$pls_pred_36 >= crit.pos  ] <- 1
pls_pred_arbres_36$crit[pls_pred_arbres_36$pls_pred_36 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.32 <- NULL

pls_pred_arbres_36[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(pls_pred_arbres_36$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

pls_pred_arbres_36$crit <- factor(pls_pred_arbres_36$crit ,levels= c(0,1))

pls_confusion_matrix_36 <- ftable( qPCR_36 ~ crit , data = pls_pred_arbres_36 )

pls_confusion_matrix_36 <- as.matrix(pls_confusion_matrix_36)

TP <- pls_confusion_matrix_36[1,1]

TN <- pls_confusion_matrix_36[2,2]

FN <- pls_confusion_matrix_36[2,1]

FP <- pls_confusion_matrix_36[1,2]

#pls_confusion_matrix_36

# Calcul des parametres pls de la matrice de confusion

Accuracy_pls36 <- ((TP+TN) / (TP+TN+FN+FP)*100)

Precision_pls36 <- ((TP / (TP+FP)*100))

Sensitivity_pls36 <- ((TP / (TP+FN)*100))

Parametre_pls_36 <- rbind(Accuracy_pls36,Precision_pls36,Sensitivity_pls36)

pls_pred_arbres_36$Pred.correct = "FAUX"
pls_pred_arbres_36$Pred.correct[pls_pred_arbres_36$code_ech_arbre %in% trueP$seuil.36] = "VRAI"

#save(pls_confusion_matrix_36,Accuracy_pls36,Precision_pls36,Sensitivity_pls36,pls_pred_arbres_36,Parametre_pls_36,  file = "Sauvegardes_objet_R.data/Pred_pls/parameter_pls_ct36.Rdata")

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif

pls_pred <- ggplot(pls_pred_arbres_36)+
  aes(y = pls_pred_36 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T7", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point(aes(shape = Pred.correct, size = Pred.correct), fill ="blue") + 
  scale_size_manual(values =c(2.1,2.9))+
  scale_shape_manual(values =c(4,21))+
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = pls_pred_36)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction pls Global a ct36  " ,subtitle = "T = positif au HLB à Ct36 en qPCR", color = "Statut") +dark_theme_gray()

pls_pred

pls_confusion_matrix_36

Parametre_pls_36

T2 <- Sys.time()

difftime(T2,T1)

#ggsave("Graphiques/Graph_pls/Prediction_pls_Global_ct36.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)


