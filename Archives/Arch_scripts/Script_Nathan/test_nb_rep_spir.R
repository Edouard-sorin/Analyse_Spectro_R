getwd()

library(ggplot2)
library(readxl)
library("imager")

library("e1071") 
library(plyr)
library(tidyr)
library(tidyverse)
library(caret)
library(pls)
library(caTools) 
library(reshape2)


#### Mise en forme du jeu de données:
spir <- read.csv(file ="Donnees/SPIR/spir_rep_trees_1_2.csv", sep=";")
obs_N1 <- as.data.frame(read_excel("Donnees/SPIR/spir arbres 1-2.xlsx", sheet = "N1")[,2:8])
obs_N2 <- as.data.frame(read_excel("Donnees/SPIR/spir arbres 1-2.xlsx", sheet = "N2")[,2:8])
obs_P1 <- as.data.frame(read_excel("Donnees/SPIR/spir arbres 1-2.xlsx", sheet = "P1")[,2:8])
obs_P2 <- as.data.frame(read_excel("Donnees/SPIR/spir arbres 1-2.xlsx", sheet = "P2")[,2:8])

obs <- list(N1=obs_N1[1:30,], N2=obs_N2[1:30,], P1=obs_P1[1:30,], P2=obs_P2[1:30,])
obs <- lapply(X = obs, FUN = function(x){replace(x, is.na(x), 0)}) # on convertit les NA en 0

obs_global <- data.frame()
for (i in 1:length(obs)) {
  obs_global <- rbind(obs_global, obs[[i]])
}
for (i in 1:ncol(obs_global)) {
  obs_global[,i] <- factor(obs_global[,i])
}

for (i in 1:nrow(spir)) {
  code_labo <- spir$X[i]
  spir$HLB[i] <- substr(code_labo, 1, 1)
  spir$rep_arbre[i] <- substr(code_labo, 2, 2)
  spir$rep_lot[i] <- substr(code_labo, 3, 3)
  spir$rep_feuille[i] <- floor( (as.numeric(substr(code_labo, 6, 8)) -1) /10) + 1
}
spir <- spir[,c(2154:2157,1:2153)]
for (i in 1:4) {
  spir[,i] <- factor(spir[,i])
}

test <- merge(obs_global,  spir, by= c("HLB", "rep_arbre", "rep_lot", "rep_feuille"), all.y = T ) # on ajoute les données d'obs sur les feuilles
for (i in 1:8) {
  test[,i] <- factor(test[,i]) # on convertit les 8 premières col d'obs en facteur
}

test <- na.omit(test)


spir <- test[,-c(8,9)]
spir$HLB <- mapvalues(spir$HLB, from=c("N", "P"), to=c(0,1))
class(spir)

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) \n et à la validation du modèle (test_set):

tirage <- do.call(rbind, lapply(1:nrow(spir), function(x) {
  data.frame(tirage= paste(unlist(spir[x,1:4]), collapse = ""))
}))
spir <- cbind(spir, tirage)

# on convertit les var pheno en numérique:
spir$chlorose <- as.numeric(as.character(spir$chlorose))
spir$chancre <- as.numeric(as.character(spir$chancre))
spir$flush <- as.numeric(as.character(spir$flush))

Nmeas <- c(1:10)
data_global <- data.frame(accuracy= rep(NA,10), sensi=rep(NA,10), pospred=rep(NA,10))
data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())
for (r in 1:10) {
  print(r)
  for (i in 1:length(Nmeas)) {
    if (i < 10) {
      set.seed(123) 
      split_five_meas = sample.split(spir$tirage , SplitRatio = Nmeas[i]/10) 
      five_meas <- subset(spir, split_five_meas == TRUE) 
      # on supprime colonne tirage: 
      five_meas <- five_meas[,-2159]
    } 
    ## on réduit le nb de mesures à x mesures:
    
    
    #↓ on crée training et test set:
    set.seed(123) 
    split = sample.split(five_meas$HLB , SplitRatio = 0.75) 
    training_set <- subset(five_meas[,-c(2:4)], split == TRUE) 
    test_set <- subset(five_meas[,-c(2:4)], split == FALSE) 
    #centrage et réduction des données
    training_set[,-1] <- scale(training_set[, -1] ) 
    test_set[, -1] <- scale(test_set[, -1] ) 
    
    
    folds = createFolds(training_set$HLB, k = 10)
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    for (j in 1:length(folds)) {
      x= folds[[j]]
      training_fold = training_set[-x, ] 
      test_fold = training_set[x, ] 
      # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
      classifier = svm(y = training_fold$HLB, x = training_fold[, -1],
                       data = training_fold,
                       type = 'C-classification',
                       kernel = 'linear')
      
      # Désormais on prédit grâce au modèle classifier le statut HLB en utilisant le test_FOLD et on calcule la précision du modèle
      y_pred = predict(classifier, newdata = test_fold[, -1])
      matrice <- confusionMatrix(test_fold$HLB, y_pred)
      accuracy <- accuracy + matrice$overall[1] /10
      sensi <- sensi + matrice$byClass[1] /10
      pospred <- pospred + matrice$byClass[3] /10
    }
    data_global[i,1] <- accuracy
    data_global[i,2] <- sensi
    data_global[i,3] <- pospred
  }
  data_global_rep <- rbind(data_global_rep, data_global)
}

data_global$nb <- seq(1:9)  

dataglobal2 <- melt(data_global,id.vars = "nb", measure.var=colnames(data_global[1:3]) )

ggplot(data=dataglobal2, aes(nb, value, fill= variable, colour= variable))+
  scale_x_continuous(breaks=seq(1:9))+
  geom_point()+
  geom_line()

library(ggplot2)
write.csv(data_global_rep, "datarep.csv")
data_global_rep$rep <- rep(seq(1:10), each=10)


