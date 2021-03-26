
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

library(snowfall) # Utilisation du calcul paralell pour optimiser la vitesse de calcul

# Importation du jeu de donnees

load("Sauvegardes_objet_R.data/data_SPIR_Ed.Rdata")

# Preparation du training_rf et test_rf  ####

nb.rep <- 6
seuil.ct <- 32

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

decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.5) # on decoupe le jeu de donné en training set et test set
train_rf <- test_ed[decoup,]
test_rf <- test_ed[!decoup,] 

# Package Random Forest ####

# Identify the predictor columns

Pred <- grep("^X", names(test_ed))


# Identity the response column

Resp = "qPCR_32"




## Essais_1 sur train et test rf ####

#set.seed(1)  # For reproducibility

# Definition of out of bag error (OOB)

# l’échantillon bootstrap, est les données choisies pour être «dans le sac» par échantillonnage avec remplacement. L’ensemble hors sac est toutes les données qui ne sont pas choisies dans le processus d’échantillonnage.L’erreur hors sac (OOB) est l’erreur moyenne pour chaque calcul à l’aide de prédictions provenant des arbres qui ne contiennent pas dans leur échantillon bootstrap respectif

# pour connaitre les paramètre d'erreus en fonction du nombre d'arbres effectués


  model_rf1 <- randomForest(
    x = train_rf[,Pred], 
    y = train_rf[,Resp],
    #xtest = test_rf[,Pred],
    ntree = 1000
  )


print(model_rf1)
plot(model_rf1)

rf_pred <- test_rf

# Prediction sur les feuilles de la base d'apprentissage

rf_pred$rf_pred_32 <- predict(model_rf1,newdata = test_rf, decision.values = T)

# Conversion des resultats de la prediction en numerique

rf_pred$rf_pred_32 = as.numeric(as.character(rf_pred$rf_pred_32))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(rf_pred))
rf_pred <- rf_pred[,c(names(rf_pred)[-select.lambda], names(rf_pred)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

rf_pred_arbres = aggregate(rf_pred_32 ~ code_ech_arbre, data = rf_pred, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.35

# Parametrage pour presentation graphique et la matrice de confusion

rf_pred_arbres$statut_pred <- "Indéterminé"
rf_pred_arbres$statut_pred[rf_pred_arbres$rf_pred_32 >= crit.pos  ] <- "Positif"
rf_pred_arbres$statut_pred[rf_pred_arbres$rf_pred_32 <= crit.neg  ] <- "Négatif"

rf_pred_arbres$pred_32 <- 0.5
rf_pred_arbres$pred_32[rf_pred_arbres$rf_pred_32 >= crit.pos  ] <- 1
rf_pred_arbres$pred_32[rf_pred_arbres$rf_pred_32 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.36 <- NULL

rf_pred_arbres[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(rf_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

rf_pred_arbres$pred_32 <- factor(rf_pred_arbres$pred_32 ,levels= c(0,1))

confusion_matrix_32 <- ftable( qPCR_32 ~ pred_32 , data = rf_pred_arbres )

confusion_matrix_32 <- as.matrix(confusion_matrix_32)

TP <- confusion_matrix_32[1,1]

TN <- confusion_matrix_32[2,2]

FN <- confusion_matrix_32[2,1]

FP <- confusion_matrix_32[1,2]

confusion_matrix_32

# Calcul des parametres rf de la matrice de confusion

Accuracy <- ((TP+TN) / (TP+TN+FN+FP)*100)

Precision <- ((TP / (TP+FP)*100))

Sensitivity <- ((TP / (TP+FN)*100))

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif


ggplot(rf_pred_arbres)+
  aes(y = rf_pred_32 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = rf_pred_32)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction rf Global a ct32  " ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +
  annotate(geom = "text", x = trueP$seuil.32, y = 0.8, label = "T")

#ggsave("Graphiques/Prédiction rf Global ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

## Essais_2 sur tte les donnes de test_ed ####

mod.complet.32 <- randomForest(
  x = test_ed[,Pred], 
  y = test_ed[,Resp],
  ntree = 1000
)

# Creation de data_pred_32

data_long_Ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(data_long_Ed[[paste0("qPCR_", seuil.ct)]] ))

feuille_spir_32 <- aggregate(data_long_Ed[c(paste0("qPCR_", seuil.ct), "reflectance")] , data_long_Ed[, c("code_ech_arbre", "code_ech_feuille","code_variete", "lambda")], mean  )

data_pred_32 <- pivot_wider(feuille_spir_32, names_from = "lambda", values_from = "reflectance", names_prefix = "X")


# Prediction sur les feuilles de la base d'apprentissage

data_pred_32$rf_pred_32 <- predict(mod.complet.32,newdata = data_pred_32, decision.values = T)

# Conversion des resultats de la prediction en numerique

data_pred_32$rf_pred_32 = as.numeric(as.character(data_pred_32$rf_pred_32))

# Re arrangemant des donnees

select.lambda <- grep("^X", names(data_pred_32))
data_pred_32 <- data_pred_32[,c(names(data_pred_32)[-select.lambda], names(data_pred_32)[select.lambda] )]

# Moyennage des resultats de la prediction pour chaque arbres

rf_pred_arbres = aggregate(rf_pred_32 ~ code_ech_arbre, data = data_pred_32, mean, na.rm = T)

# Critere choisi sur les observations graphiques

crit.pos <- 0.4
crit.neg <- 0.3

# Parametrage pour presentation graphique et la matrice de confusion

rf_pred_arbres$statut_pred <- "Indéterminé"
rf_pred_arbres$statut_pred[rf_pred_arbres$rf_pred_32 >= crit.pos  ] <- "Positif"
rf_pred_arbres$statut_pred[rf_pred_arbres$rf_pred_32 <= crit.neg  ] <- "Négatif"

rf_pred_arbres$pred_32 <- 0.5
rf_pred_arbres$pred_32[rf_pred_arbres$rf_pred_32 >= crit.pos  ] <- 1
rf_pred_arbres$pred_32[rf_pred_arbres$rf_pred_32 <= crit.neg  ] <- 0

# Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

trueP$seuil.36 <- NULL

rf_pred_arbres[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(rf_pred_arbres$code_ech_arbre %in% x ))

# Résultats de la prédiction sous forme de matrice de confusion 

rf_pred_arbres$pred_32 <- factor(rf_pred_arbres$pred_32 ,levels= c(0,1))

confusion_matrix_32 <- ftable( qPCR_32 ~ pred_32 , data = rf_pred_arbres )

confusion_matrix_32 <- as.matrix(confusion_matrix_32)

TP <- confusion_matrix_32[1,1]

TN <- confusion_matrix_32[2,2]

FN <- confusion_matrix_32[2,1]

FP <- confusion_matrix_32[1,2]

# Calcul des parametres rf de la matrice de confusion

Accuracy <- TP+TN / TP+TN+FN+FP

Precision <- TP / TP+FP

Sensitivity <- TP / TP+FN

# présentation de la Matrice de confusion 

Confusion_matrix <-  matrix(c("True Positive (TP)","False Positive (FP)(Type error 1)","False Negative (FN)(Type error 2)","True Negative (TN)"), ncol = 2 , byrow = T)

colnames(Confusion_matrix) <- c("Actual Positive","Actual Negative")
rownames(Confusion_matrix) <- c("Predicted Positive","Predicted Negative")

Confusion_matrix

# Organiser resultat dans un tableau à 4 cases avec Vrai postif, Vrai négatif , Faux positif, Faut négatif


ggplot(rf_pred_arbres)+
  aes(y = rf_pred_32 , x = code_ech_arbre, color = statut_pred) +
  geom_rect(aes(xmin = "A1", xmax = "T+", ymin = crit.neg, ymax = crit.pos), lwd = 1 , fill = "gray50", color = "gray50") +
  geom_point() + 
  geom_segment(aes(xend = code_ech_arbre, y = 0, yend = rf_pred_32)) +
  scale_color_manual(values = c(Négatif = "red", Indéterminé = "gray20", Positif = "green")) +
  labs(x = "Code de l'arbre", y = "Résultat moyen des feuilles entre 0 et 1", title = "Prédiction rf Global a ct32  " ,subtitle = "T = positif au HLB à Ct32 en qPCR", color = "Statut") +
  annotate(geom = "text", x = trueP$seuil.32, y = 0.8, label = "T")

#ggsave("Graphiques/Prédiction rf Global ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

rm (feuille_spir_32,mod.complet.32,select.lambda,crit.pos,crit.neg,seuil.ct,data_SPIR_Ed,data_long_Ed,trueP)

## Essais_3 Calcul paralell pour obtenir la meilleur ####

Mtry <- 1:1000  # Nombre de variables échantillonnées au hasard en tant que candidats à chaque fractionnement (noeud)

#  5.586323 hours de temps de calcul pour Mtry <- 1:1000  à 1000 arbres

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
                              # la library des packages utilisés
sfLibrary(randomForest)
sfLibrary(caret)
sfExport("Pred","Resp","Mtry","train_rf") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser


#res.rf.32 <- sfClusterApplyLB(Mtry , function(x) median(randomForest(train_rf[,Resp]~ . , data=train_rf[,Pred], ntree= 400, proximity = T, mtry = x)$err.rate[,"OOB"])) 

Erreur_OOB_mediane <- sfSapply(Mtry , function(x) median(randomForest(train_rf[,Resp]~ . , data=train_rf[,Pred], ntree= 1000, proximity = T, mtry = x)$err.rate[,"OOB"])) 

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes

difftime(T2,T1) # temps de calcul

#test.RF <- sapply(1:15, function(x) median(randomForest(train_rf[,Resp]~ . , data=train_rf[,Pred], ntree= 400, proximity = T, mtry = x)$err.rate[,"OOB"]))

#names(res.rf.32) <- rep

res.rf.32 <- as.data.frame(Erreur_OOB_mediane) 

res.rf.32$Mtry <- Mtry

plot(Erreur_OOB_mediane, pch = 16, las = 1, ylab = "erreur OOB médiane", xlab = "Valeur de Mtry")

#esquisse::esquisser() 

lines(Erreur_OOB_mediane)
abline(v =11, col = "red", lty = 2)
dev.off()

#save(Erreur_OOB_mediane,res.rf.32,  file = "Sauvegardes_objet_R.data/Erreur_OOB_mediane_RF_1000simu_1000arbres.Rdata")

#write.table(x = res.rf.32  , file = "Donnees/Erreur_OOB_mediane_RF_1000simu_1000arbres.csv" , sep = ';')

#load("Sauvegardes_objet_R.data/Erreur_OOB_mediane_RF_1000simu_1000arbres.Rdata")


g1 <- ggplot(data = res.rf.32) +
  aes(x = Mtry, y = Erreur_OOB_mediane)+
  #stat_summary(geom = "pointrange", fun.data = function(x) mean_se(x, mult = qt(0.975, length(x) - 1))) +
  #stat_summary(geom = "line", fun = mean) +
  geom_point(lwd = 1.1) +
  labs(x = "Nombre de Mtry effectué ", y = "Erreur_OOB_mediane", title = "Prédiction des paramètres erreur OOB medianede en RF ", subtitle = "Détection du statut HLB à (Ct<32), obtenu après avoir fait tournée 1000 Mtry à 1000 arbres") +
  dark_theme_gray() +
  theme(legend.position = "right")+
  #scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))

#ggsave("Graphiques/Graph_RF/Erreur_OOB_mediane_RF_1000simu_1000arbres_points.png",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)


## Essais_4 model_carret ####
 
model_carret <- train(
  x = train_rf[,Pred], 
  y = train_rf[,Resp],
  method = "parRF",
  preProcess = NULL,
  weights = NULL,
  metric = "Accuracy",
  maximize = TRUE,
  trControl = trainControl(),
  tuneGrid = NULL,
  tuneLength = 3
)

rf_pred_feuille <-  data.frame(predict(model_carret, newdata = test_rf[,grep("^X", names(test_rf))], decision.values = T))

# Package party ####
# Essais_5 Arbre de decision sur l'ensemble des donnees

gtree_glob <- ctree( test_ed[,Pred] ~ . ,  data=test_ed[[paste0("qPCR_", seuil.ct)]] )


RF_ct32 <- data_SPIR_Ed
RF_ct32$qPCR_32 <- as.numeric(RF_ct32$qPCR_32)
masque_numeric <- sapply(RF_ct32, is.numeric)
RF_ct32 <- RF_ct32[,masque_numeric]

gtree_glob <- ctree( qPCR_32 ~ . , data=RF_ct32)

print(gtree_glob)

# plot 

plot(gtree_glob ,)

plot(gtree_glob, inner_panel = node_barplot,
     edge_panel = function(ctreeobj, ...) { function(...) invisible() },
     tnex = 1)

nodes(gtree_glob, 1)


