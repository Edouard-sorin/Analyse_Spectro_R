rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(e1071) # SVM

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(caTools) # sample.split
library(tidyr) # pivot longer & pivot wider


# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Preparation du training_svm et test_svm ####


nb.rep <- 6
seuil.ct <- 32

Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)


maliste <- lapply(Tirage,function(feuille){ 
  list_svm <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
  # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
  code <- grep("^code_ech", names(list_svm))
  qPCR <- grep("^qPCR_", names(list_svm))
  sortie <- cbind.data.frame(unique(list_svm[c(code,qPCR)])
                             , matrix(apply(list_svm[-(which(colnames(list_svm) == "code_variete"):which(colnames(list_svm) == "qPCR_36"))], 2, mean)
                                      , nr = 1, dimnames = list(NULL, names(list_svm)[-(which(colnames(list_svm) == "code_variete"):which(colnames(list_svm) == "qPCR_36"))])))
  sortie
})

test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 

test_ed <- test_ed[-(which(colnames(test_ed) == "qPCR_36"))]

test_ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(test_ed[[paste0("qPCR_", seuil.ct)]] ))


decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.5) # on decoupe le jeu de donné en training set et test set
train_svm <- test_ed[decoup,]
test_svm <- test_ed[!decoup,] 


# Prediction a partir du training set et du test set a ct32 ####

#set.seed(1)  # For reproducibility


model_SVM_32 <- svm(  y = train_svm[[paste0("qPCR_", seuil.ct)]]       #   on fix la valeur de ct à predir    
                      , x = train_svm[, grep("^X", names(train_svm))] #   on prend ttes les longueurs d'ondes
                      , type = 'C-classification'
                      , kernel = 'linear'
) 


print(model_SVM_32)

svm_pred <- test_svm


# Prediction sur les feuilles de la base d'apprentissage

svm_pred$svm_pred_32 <- predict(model_SVM_32 ,newdata = test_svm[, grep("^X", names(test_svm))] , decision.values = T)

# Conversion des resultats de la prediction en numerique

svm_pred$svm_pred_32 = as.numeric(as.character(svm_pred$svm_pred_32))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(svm_pred))
svm_pred <- svm_pred[,c(names(svm_pred)[-select.lambda], names(svm_pred)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

svm_pred_arbres_32 = aggregate(svm_pred_32 ~ code_ech_arbre, data = svm_pred, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.35

# Parametrage pour presentation graphique et la matrice de confusion

svm_pred_arbres_32$statut_pred <- "Indéterminé"
svm_pred_arbres_32$statut_pred[svm_pred_arbres_32$svm_pred_32 >= crit.pos  ] <- "Positif"
svm_pred_arbres_32$statut_pred[svm_pred_arbres_32$svm_pred_32 <= crit.neg  ] <- "Négatif"

svm_pred_arbres_32$crit <- 0.5
svm_pred_arbres_32$crit[svm_pred_arbres_32$svm_pred_32 >= crit.pos  ] <- 1
svm_pred_arbres_32$crit[svm_pred_arbres_32$svm_pred_32 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.36 <- NULL

svm_pred_arbres_32[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(svm_pred_arbres_32$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

svm_pred_arbres_32$crit <- factor(svm_pred_arbres_32$crit ,levels= c(0,1))

svm_confusion_matrix_32 <- ftable( qPCR_32 ~ crit , data = svm_pred_arbres_32 )

svm_confusion_matrix_32 <- as.matrix(svm_confusion_matrix_32)

TP <- svm_confusion_matrix_32[1,1]

TN <- svm_confusion_matrix_32[2,2]

FN <- svm_confusion_matrix_32[2,1]

FP <- svm_confusion_matrix_32[1,2]

svm_confusion_matrix_32

# Calcul des parametres svm de la matrice de confusion

Accuracy_svm32 <- ((TP+TN) / (TP+TN+FN+FP)*100)

Precision_svm32 <- ((TP / (TP+FP)*100))

Sensitivity_svm32 <- ((TP / (TP+FN)*100))

Parametre_svm_32 <- rbind(Accuracy_svm32,Precision_svm32,Sensitivity_svm32)

svm_pred_arbres_32$Vrai = "FAUX"
svm_pred_arbres_32$Vrai[svm_pred_arbres_32$code_ech_arbre %in% trueP$seuil.32] = "VRAI"

#save(svm_confusion_matrix_32,Accuracy_svm32,Precision_svm32,Sensitivity_svm32,svm_pred_arbres_32,Parametre_svm_32,  file = "Sauvegardes_objet_R.data/Pred_svm/parameter_svm_ct32.Rdata")

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif


svm_pred <- ggplot(svm_pred_arbres_32)+
  aes(y = svm_pred_32 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point(aes(shape = Vrai, size = Vrai), fill ="blue") + 
  scale_size_manual(values =c(2.1,2.9))+
  scale_shape_manual(values =c(4,21))+
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = svm_pred_32)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction svm svm Global a ct32  " ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +dark_theme_gray()

svm_pred

#ggsave("Graphiques/Graph_pred_SVM-ct32/Prediction_SVM_Global_ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)
