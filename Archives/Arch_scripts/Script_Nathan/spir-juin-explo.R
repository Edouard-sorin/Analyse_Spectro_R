



# Chargement des packages: ------------------------------------------------

library(ggplot2)
library(readxl)
library("imager")
library(data.table)
library(plyr)
library(tidyr)
library(tidyverse)
library(caret)
library(pls)
library(caTools) 
library(reshape2)
library(e1071)
library(gt)
library(ROCR)

# Mise en forme du jeu de données: ----------------------------------------

spir <- read_excel(path = "Donnees/Donnees_Nathan/SPIRjuin.xlsx")
spir <- as.data.frame(spir)
mypath_var_pheno="Donnees/SPIR/spir2-juin/var_pheno"#select the path where your resu data files are
filenames_var_pheno<-list.files(path = mypath_var_pheno, full.names = TRUE)
var_pheno_all<-rbindlist(lapply(filenames_var_pheno,read.csv2))
var_pheno_all[is.na(var_pheno_all)] <- 0
var_pheno_all <- var_pheno_all[,-8]
var_pheno_all <- lapply(var_pheno_all, factor)
var_pheno_all <- as.data.frame(var_pheno_all)
colnames(var_pheno_all)[1] <- "hyp_HLB"
var_pheno_all$chlorose[which(var_pheno_all$chlorose==11)] <- 1


# Identification des lots positifs ----------------------------------------
# plusieurs conditions: - CT <32 ou 36
#                       - au moins 1 lot + sur les 3 => les 3 lots sont positifs ou non

qPCR <- read.csv(file=  "Donnees/Donnees_Nathan/results_qPCR.csv", sep=";")
trueP <- list(levels(factor(qPCR$Sample.Name[which(qPCR$C..Mean <32)])), levels(factor(qPCR$Sample.Name[which(qPCR$C..Mean <36)])) ) # on stocke les noms d'échantillons positifs selon nos 2 seuils de Ct

arbres <- list(c(), c() )
for (a in 1:length(arbres)) {
  for (i in 1:length(trueP[[a]])) {  arbres[[a]] <- append(arbres[[a]], substr(trueP[[a]][i], 1,2))}
  arbres[[a]] <- unique(arbres[[a]])
}
trueP[[3]] <- arbres[[1]]
trueP[[4]] <- arbres[[2]]
remove(arbres)
names(trueP) <- c("ech_ct<32_lot","ech_ct<36_lot", "ech_ct<32_arbre","ech_ct<36_arbre" )


testHLB <- list(data.frame(hyp_HLB= rep(c(rep("N",3), rep("P", 3)), 7), rep_arbre= rep(1:7, each=6), rep_lot= rep(1:3, times=14), qPCR_HLB= NA))
testHLB[[1]]$code_ech <- paste(testHLB[[1]]$hyp_HLB, testHLB[[1]]$rep_arbre, testHLB[[1]]$rep_lot, sep="")
testHLB[[2]] <- testHLB[[1]]
testHLB[[3]] <- testHLB[[1]]
testHLB[[4]] <- testHLB[[1]]
names(testHLB) <- names(trueP)
for (l in 1:2) { #on applique la boucle à chaque élément de testHLB, soit on détermine les posi selon nos deux seuils de Ct
  for (i in 1:nrow(testHLB[[l]])) {
    if ( testHLB[[l]]$code_ech[i] %in%  trueP[[l]]){ 
      testHLB[[l]]$qPCR_HLB[i] <- 1}
    else { testHLB[[l]]$qPCR_HLB[i] <- 0}
  }
  testHLB[[l]] <- lapply(testHLB[[l]], factor)
  testHLB[[l]] <- as.data.frame(testHLB[[l]])
}
for (l in 3:4) { #on applique la boucle à chaque élément de testHLB, soit on détermine les posi selon nos deux seuils de Ct
  for (i in 1:nrow(testHLB[[l]])) {
    if ( substr(testHLB[[l]]$code_ech[i], 1, 2) %in%  trueP[[l]]){ 
      testHLB[[l]]$qPCR_HLB[i] <- 1}
    else { testHLB[[l]]$qPCR_HLB[i] <- 0}
  }
  testHLB[[l]] <- lapply(testHLB[[l]], factor)
  testHLB[[l]] <- as.data.frame(testHLB[[l]])
}





for (i in 1:nrow(spir)) {
  code_labo <- spir$code_labo[i]
  code_labo <- str_replace(code_labo, "NEGA", "N")
  code_labo <- str_replace(code_labo, "HLB", "P")
  spir$hyp_HLB[i] <- substr(code_labo, 1, 1)
  spir$rep_arbre[i] <- substr(code_labo, 2, 2)
  spir$rep_lot[i] <- substr(code_labo, 3, 3)
  spir$rep_feuille[i] <- floor((as.numeric(substr(code_labo, 6, 8))) /10) + 1
}
spir <- spir[,c(2154:2157,1:2153)]
spir$code_ech <- factor(paste(spir$hyp_HLB, spir$rep_arbre, spir$rep_lot, sep=""))

for (i in 1:length(testHLB)) { #on joint le résultat du test HLB au jdd de spir (on a 4 colonnes de test HLB: voir leurs noms)
  spir <- merge(spir, testHLB[[i]][4:5], by= "code_ech", all.y= F, suffixes = c("", names(testHLB)[i]) )
}
colnames(spir)[2159] <- "qPCR_HLBech_ct<32_lot" # pourquoi le suffixe ne fonctionne pas pour le premier merge avec i=1 ?



spir$hyp_HLB <- factor(spir$hyp_HLB)
levels(spir$hyp_HLB) <- c("0", "1")
spir <- merge(var_pheno_all,  spir, by= c("hyp_HLB", "rep_arbre", "rep_lot", "rep_feuille"), all.y = T ) # on ajoute les données d'obs sur les feuilles
spir <- spir[,-c(1:4,8:10)]
colnames(spir)[3] <- "flush"


# prédiction du statut HLB ------------------------------------------------


#on convertit les var pheno en numérique:
spir$chlorose <- as.numeric(as.character(spir$chlorose))
spir$chancre <- as.numeric(as.character(spir$chancre))
spir$flush <- as.numeric(as.character(spir$flush))
  
spirsave <- spir #on sauvegarde une version de spir
spirsave_2 <- spirsave
#on essaie en ayant juste 2 modalités pour les varphéno (partie à ne pas exécuter si on veut prédire avec les var phéno en catégorielles)
for (i in 1:3) {
  for (j in 1:nrow(spirsave)) {
    if (spirsave[j,i] >= 1) {
      spirsave[j,i] <- 1
    }
    else{
      spirsave[j,i] <- 0
    }
  }
}

dataglobal_binaire <- data.frame(condition= rep(NA,4), accuracy= rep(NA,4), kappa= rep(NA,4), speci= rep(NA,4), sensi= rep(NA,4), pospred= rep(NA,4), negpred= rep(NA,4)) # on crée le dataframe qui contiendra les performances des classif pour chaque conditions

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) et à la validation du modèle (test_set):
for (c in 1:4) { # nos 4 conditions (2 sur le Ct seuil * 2 sur lot/arbre considéré en positif)
  print(c)
  spir <- spirsave[,c(1:2154, 2154+c)] # on construit le jdd avec les variables explicatives et la variable à expliquer (2154+i)ème colonne de 'spirsave'
  colnames(spir)[2155] <- "HLB" #on renomme la dernière colonne
  spir <- spir[,c(2155, 1:2154)] # on réorganise l'ordre des colonnes
  # split = sample.split(spir$HLB, SplitRatio = 0.75) 
  # training_set <- subset(spir, split == TRUE) 
  # test_set <- subset(spir, split == FALSE) 
  training_set <- spir # en fait on n'a pas besoin de diviser le jddd vu qu'on fait des cross validation derrière!!!
  training_set[c(2:2155)] <- scale(training_set[c(2:2155)] ) #centrage et réduction des données
  test_set[c(2:2155)] <- scale(test_set[c(2:2155)] ) 
  
  folds = createFolds(training_set$HLB, k = 10)
  accuracy <- 0
  kappa <- 0
  speci <- 0
  pospred <- 0
  negpred <- 0
  sensi <- 0
  
  for (j in 1:length(folds)) {
    print(j)
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
    matrice <- confusionMatrix(test_fold$HLB, y_pred, positive = '1')
    accuracy <- accuracy + matrice$overall[1] /10
    kappa <- kappa + matrice$overall[2] /10
    sensi <- sensi + matrice$byClass[1] /10
    speci <- speci + matrice$byClass[2] /10
    pospred <- pospred + matrice$byClass[3] /10
    negpred <- negpred + matrice$byClass[4] /10
  }
  dataglobal_binaire$condition[c] <-  colnames(spirsave)[2154+c]
  dataglobal_binaire$accuracy[c] <- accuracy
  dataglobal_binaire$sensi[c] <- sensi
  dataglobal_binaire$pospred[c] <- pospred
  dataglobal_binaire$kappa[c] <- kappa
  dataglobal_binaire$speci[c] <- speci
  dataglobal_binaire$negpred[c] <- negpred
}

library(dplyr)
library(gt)


dataglobal_binaire$'CT seuil' <- rep(c(32,36), 2)
dataglobal_binaire$'Lots du même arbre indépendants' <- rep(c("oui", "non"), each=2)
dataglobal_binaire <- dataglobal_binaire[,-1]
dataglobal_binaire <- dataglobal_binaire[, c(7,8, 1:6)]

dataglobal_categorie$'CT seuil' <- rep(c(32,36), 2)
dataglobal_categorie$'Lots du même arbre indépendants' <- rep(c("oui", "non"), each=2)
dataglobal_categorie <- dataglobal_categorie[,-1]
dataglobal_categorie <- dataglobal_categorie[, c(7,8, 1:6)]

dataglobal_binaire %>%
  gt() %>%
  tab_header(
    title = "Evaluation de la prédiction du statut HLB",
    subtitle = "Variables phénotypiques binaires")
dataglobal_categorie %>%
  gt() %>%
  tab_header(
    title = "Evaluation de la prédiction du statut HLB",
    subtitle = "Variables phénotypiques catégorielles")




gt::  
tab_header(data= dataglobal_binaire,
    title = "Evaluation de la prédiction du statut HLB",
    subtitle = "Variables phénotypiques binaires"
  )
gt(data= dataglobal_binaire)






# prédiction sans les var phénotypiques -----------------------------------


dataglobal_spectre <- data.frame(condition= rep(NA,4), accuracy= rep(NA,4), kappa= rep(NA,4), speci= rep(NA,4), sensi= rep(NA,4), pospred= rep(NA,4), negpred= rep(NA,4)) # on crée le dataframe qui contiendra les performances des classif pour chaque conditions

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) et à la validation du modèle (test_set):
for (c in 1:4) { # nos 4 conditions (2 sur le Ct seuil * 2 sur lot/arbre considéré en positif)
  print(c)
  spir <- spirsave[,c(1:2154, 2154+c)] # on construit le jdd avec les variables explicatives et la variable à expliquer (2154+i)ème colonne de 'spirsave'
  colnames(spir)[2155] <- "HLB" #on renomme la dernière colonne
  spir <- spir[,c(2155, 1:2154)] # on réorganise l'ordre des colonnes
  # split = sample.split(spir$HLB, SplitRatio = 0.75) 
  # training_set <- subset(spir, split == TRUE) 
  # test_set <- subset(spir, split == FALSE) 
  training_set <- spir # en fait on n'a pas besoin de diviser le jddd vu qu'on fait des cross validation derrière!!!
  training_set[c(2:2155)] <- scale(training_set[c(2:2155)] ) #centrage et réduction des données
  test_set[c(2:2155)] <- scale(test_set[c(2:2155)] )
  
  folds = createFolds(training_set$HLB, k = 10)
  accuracy <- 0
  kappa <- 0
  speci <- 0
  pospred <- 0
  negpred <- 0
  sensi <- 0
  
  for (j in 1:length(folds)) {
    print(j)
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
    matrice <- confusionMatrix(test_fold$HLB, y_pred, positive = '1')
    accuracy <- accuracy + matrice$overall[1] /10
    kappa <- kappa + matrice$overall[2] /10
    sensi <- sensi + matrice$byClass[1] /10
    speci <- speci + matrice$byClass[2] /10
    pospred <- pospred + matrice$byClass[3] /10
    negpred <- negpred + matrice$byClass[4] /10
  }
  dataglobal_spectre$condition[c] <-  colnames(spirsave)[2154+c]
  dataglobal_spectre$accuracy[c] <- accuracy
  dataglobal_spectre$sensi[c] <- sensi
  dataglobal_spectre$pospred[c] <- pospred
  dataglobal_spectre$kappa[c] <- kappa
  dataglobal_spectre$speci[c] <- speci
  dataglobal_spectre$negpred[c] <- negpred
}
dataglobal_spectre$'CT seuil' <- rep(c(32,36), 2)
dataglobal_spectre$'Lots du même arbre indépendants' <- rep(c("oui", "non"), each=2)
dataglobal_spectre <- dataglobal_spectre[,-1]
dataglobal_spectre <- dataglobal_spectre[, c(7,8, 1:6)]

dataglobal_spectre %>%
  gt() %>%
  tab_header(
    title = "Evaluation de la prédiction du statut HLB",
    subtitle = "Données spectrales uniquement")



# importance des variables ------------------------------------------------

library(rminer)
M <- fit(HLB~., data= training_set, model="svm", kpar=list(sigma=0.10), C=2)
svm.imp <- Importance(M, data=training_set)

L=list(runs=1,sen=t(svm.imp$imp),
      sresponses=svm.imp$sresponses)

mgraph(L,graph="IMP",leg=names(training_set),col="gray",Grid=0, sort=TRUE)

importance <-  data.frame(imp= svm.imp$imp, variable= colnames(spir))
importance_lambda <- importance[5:2155,]
importance_pheno <- importance[2:4,]
importance_pheno$variable <- factor(importance_pheno$variable)
importance_lambda$variable <- as.numeric(as.character(importance_lambda$variable))
importance_lambda <- importance_lambda[order(importance_lambda$variable, decreasing = F),]
most_important <- importance_lambda[1:100,]

ggplot(importance_lambda, aes(x=variable, y=imp)) +
  labs(x="Longueurs d'onde (en nm)", y= "Importance")+
  geom_col( fill="orange", width = 0.5) +
  scale_x_continuous(breaks = seq(350, 2500, by=100)) +
  labs(title = "Importance des valeurs de réflectance pour chaque longueur d'onde")+
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size= 7.5), panel.background = element_blank())

ggplot(importance_pheno, aes(x=variable, y=imp, fill=variable)) +
  labs(x="Variable phénotypique", y= "Importance")+
  geom_col( width=0.5) +
  ylim (0, 0.0006)+
  scale_fill_manual(values = c("#CD853F", "#FFD700", "#00CD00"))+
  labs(title = "Importance des variables phénotypiques")+
  theme(axis.line = element_line(colour = "black"),  panel.background = element_blank())



# prédiction avec seulement les var phénotypique --------------------------

pheno_only <- data.frame(condition= rep(NA,4), Accuracy= rep(NA,4),  Precision= rep(NA,4), Specificity= rep(NA,4), Error_rate= rep(NA,4),  Sensitivity= rep(NA,4),auc= rep(NA,4)) # on crée le dataframe qui contiendra les performances des classif pour chaque conditions
resu_table_global <- list()
# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) et à la validation du modèle (test_set):
for (c in 1:4) { # nos 4 conditions (2 sur le Ct seuil * 2 sur lot/arbre considéré en positif)
  print(c)
  spir <- spirsave[,c(1:3,2154 + c)] # on construit le jdd avec les variables explicatives et la variable à expliquer (2154+i)ème colonne de 'spirsave'
  colnames(spir)[4] <- "HLB" #on renomme la dernière colonne
  spir <- spir[,c(4, 1:3)] # on réorganise l'ordre des colonnes et on garde que les var phénotypiques (colonnes 1 à 3)
  # split = sample.split(spir$HLB, SplitRatio = 0.75) 
  # training_set <- subset(spir, split == TRUE) 
  # test_set <- subset(spir, split == FALSE) 
  training_set <- spir
  training_set[c(2:ncol(spir))] <- scale(training_set[c(2:ncol(spir))] ) #centrage et réduction des données
  folds = createFolds(training_set$HLB, k = 10)
  # accuracy < <-  0
  # kappa <- 0
  # speci <- 0
  # pospred <- 0
  # negpred <- 0
  # sensi <- 0
  auc <- 0
  resu_table <- ftable (factor(c(),levels=c(0,1)), y_pred=factor(c(),levels=c(0,1)))
  
  for (j in 1:length(folds)) {
    print(j)
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
    pred_ROCR <- prediction(predictions = as.numeric(as.character(y_pred)), labels = as.numeric(as.character(test_fold[,1])) )
    auc_ROCR <- performance(pred_ROCR, measure = "auc")
    auc_ROCR <- auc_ROCR@y.values[[1]]
    
    resu <- ftable(test_fold$HLB, y_pred)
    # matrice <- confusionMatrix(test_fold$HLB, y_pred, positive = '1')
    # accuracy <- accuracy + matrice$overall[1] /10
    # kappa <- kappa + matrice$overall[2] /10
    # sensi <- sensi + matrice$byClass[1] /10
    # speci <- speci + matrice$byClass[2] /10
    # pospred <- pospred + matrice$byClass[3] /10
    # negpred <- negpred + matrice$byClass[4] /10
    
    auc <- auc + auc_ROCR/10
    resu_table <- resu_table + resu/10 
  }
  pheno_only$condition[c] <-  colnames(spirsave)[2154+c]
  # pheno_only$accuracy[c] <- accuracy
  # pheno_only$sensi[c] <- sensi
  # pheno_only$pospred[c] <- pospred
  # pheno_only$kappa[c] <- kappa
  # pheno_only$speci[c] <- speci
  # pheno_only$negpred[c] <- negpred
  pheno_only$auc[c] <- round(auc, digits=3)
  resu_table_global[[c]] <- resu_table 
}
for (i in 1:4) {
  toto <-  resu_table_global[[i]]
  pheno_only$Accuracy[i] <-round((toto[2,2]+toto[1,1])/ sum(toto), digits=3)
  pheno_only$Precision[i] <- round(toto[2,2]/sum(toto[,2]), digits=3)
  pheno_only$Specificity[i] <- round(toto[1,1]/sum(toto[1,]), digits=3)
  pheno_only$Error_rate[i] <- round((toto[1,2]+toto[2,1])/sum(toto),digits=3)
  pheno_only$Sensitivity[i] <- round(toto[2,2] / sum(toto[2,]), digits=3)
}

pheno_only$CTseuil <- c(32,36,32,36)
pheno_only$lots_indépendants <- rep(c("oui", "non (arbre)"), each=2)
pheno_only <- pheno_only[,-1]

pheno_only %>%
  gt() %>%
  tab_header(
    title = "Evaluation de la prédiction du statut HLB",
    subtitle = "variables phénotypiques uniquement")


# Visualisation des graphes: ----------------------------------------------

##spectre moyen 350-2500 nm:
spir <- spirsave
nega <- c()
posi <- c()
lambda <- c()
for (i in (4:(ncol(spir)-4)) ) {
  calcul <- aggregate( formula= spir[,i] ~ spir[,(ncol(spir)-2)], data= spir , FUN = mean) #on prédit en fonction du statut HLB ct32 et arbre
  nega <- append( nega, calcul[1,2] )
  posi <- append( posi, calcul[2,2] )
  lambda <- append( lambda, as.numeric(colnames(spir)[i] ) )}


spectra_N <- data.frame(lambda, HLB=factor("N"), Tr= nega)
spectra_P <- data.frame(lambda, HLB=factor("P"), Tr= posi)

mean_spectra_N_P <- rbind.data.frame(spectra_N, spectra_P)


graph1 <- ggplot(mean_spectra_N_P, aes(x= lambda, y=Tr, color=HLB)) + 
  geom_point(size=1) +
  scale_color_manual( values = c("#00CD00", "#CD3333"), labels=c("Négatif", "Positif"))+
  labs(title = "Spectre moyen en fonction du statut HLB +/- des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est positif (Ct<32)",x="Longeur d'onde (en nm)", y="Absorbance", color= "Résultat du test HLB" ) +
  scale_x_continuous(breaks = seq(350, 2500, by=100))+
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size= 7.5), panel.background = element_blank())

#zoom
graph2 <- ggplot(mean_spectra_N_P, aes(x= lambda, y=Tr, color=HLB)) + 
  geom_point(size=1) +
  scale_color_manual( values = c("#00CD00", "#CD3333"), labels=c("Négatif", "Positif"))+
  labs(title = "Spectre moyen en fonction du statut HLB +/- des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est positif (Ct<32)",x="Longeur d'onde (en nm)", y="Absorbance", color= "Résultat du test HLB" ) +
  scale_x_continuous(limits = c(400,680),breaks = seq(400, 650, by=50))+
  scale_y_continuous(limits=c(NA, 0.15))+
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size= 7.5), panel.background = element_blank())

##Lzoom2
graph3 <- ggplot(mean_spectra_N_P, aes(x= lambda, y=Tr, color=HLB)) + 
  geom_point(size=1) +
  scale_color_manual( values = c("#00CD00", "#CD3333"), labels=c("Négatif", "Positif"))+
  labs(title = "Spectre moyen en fonction du statut HLB +/- des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est positif (Ct<32)",x="Longeur d'onde (en nm)", y="Absorbance", color= "Résultat du test HLB" ) +
  scale_x_continuous(limits = c(700,1000), breaks = seq(700, 1000, by=50)) +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size= 7.5), panel.background = element_blank())

graph4 <- ggplot(mean_spectra_N_P, aes(x= lambda, y=Tr, color=HLB)) + 
  geom_point(size=1) +
  scale_color_manual( values = c("#00CD00", "#CD3333"), labels=c("Négatif", "Positif"))+
  labs(title = "Spectre moyen en fonction du statut HLB +/- des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est positif (Ct<32)",x="Longeur d'onde (en nm)", y="Absorbance", color= "Résultat du test HLB" ) +
  scale_x_continuous(limits = c(900,1050), breaks = seq(900, 1050, by=50)) +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size= 7.5), panel.background = element_blank())



# ##chlorose:
# nega <- c()
# posi <- c()
# lambda <- c()
# chlorose <- c()
# 
# for (i in (4:(ncol(spir)-4)) ) {
#   # print(i)
#   calcul <- aggregate( formula= spir[,i] ~ spir[,(ncol(spir)-2)] + spir$chlorose, data= spir , FUN = mean) 
#   nega <- append(nega, calcul[which(calcul[,1]==0), 3])
#   posi<- append(posi, calcul[which(calcul[,1]==1), 3])
#   chlorose <- append(chlorose, c(1,2,3))
#   lambda <- append( lambda, as.numeric(colnames(spir)[i] ) )}
# 
# spectra_N <- data.frame(lambda, HLB=factor("N"), Tr= nega, chlorose= chlorose)
# spectra_P <- data.frame(lambda, HLB=factor("P"), Tr= posi)
# 
# mean_spectra_N_P <- rbind.data.frame(spectra_N, spectra_P)
# 



# Analyse tableau de résultats prédiction ---------------------------------

tablx <- read.csv("resultatsSPIRaucSVM.csv", sep=";")
names(tablx)[13] <- "poids_var_pheno"
tablx <- tablx[order(tablx$AUC, decreasing = T),]
for (j in 3:9){
  for (i in 1:nrow(tablx)) {
    tablx[i,j] <- round(tablx[i,j], digits=3)
  }
}

gt_tablx <-  tablx[,3:13 ] %>%
  gt() %>%
  tab_header(
    title = "Evaluation de la prediction du statut HLB des spectres proches infra-rouge",
    subtitle = "Modele SVM, kernel linéaire")

gt_tablx

