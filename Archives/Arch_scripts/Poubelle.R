



# selec <- selec[,-(which(colnames(selec) == "num"))]

A = substr(selec$rep_arbre, 1, 1)
B = substr(selec$rep_arbre, 2, 2)
selec$rep_arbre3 = "a"
for( i in 1:length(selec[,3])){
  
  if(A[i] == 2){
    selec$rep_arbre3[i] = paste("N", B[i], sep = "") 
  } 
  
  if(A[i] == 1){
    selec$rep_arbre3[i] = paste("P", B[i], sep = "")
  } 
  
  
}



b = list()
for (j in 1:length(selec[1,])){
  
  b = c(selec[,j], b)
  
}




# selec$hyp_arbre <- as.numeric(substr(selec$rep_arbre, 2, 2))
# 
# selec$nom_arbre <- as.numeric(substr(selec$code_ech, 1, 1))
# 
# selec$code_ech_arbre <- factor(paste(selec$nom_arbre, selec$hyp_arbre, sep=""))  # creation d'un code echantillon  numeric POUR LES ARBRESpour le passage en format long 
# 
# selec <- selec[, -(which(colnames(selec) =="hyp_arbre" ))]


rm(list=ls())

library(ggplot2)
library(ggdark)
library(esquisse)
library(colourpicker)
library(readxl)
library(imager)  
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
library(dplyr)
library(rminer)

# Importation des donnees data_SPIR_Ed ####

data_SPIR_Ed <- read.table(file = "Donnees/Donnees_Nathan/SPIRjuin_ED.csv"
                           
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )

# Attention des qu'il y a du texte(ou des virgules) tt le jeu de donnee est compris comme etant non numerique , il faut donc le mettre en numerique !


# data_SPIR_Ed
data_SPIR_Ed$X <- rownames(data_SPIR_Ed)
for (i in 1:nrow(data_SPIR_Ed )) {
  code_labo <- data_SPIR_Ed$X[i]                   # X devient temporairement code_labo
  code_ech <- str_replace(code_labo, "NG", "N")
  code_ech <- str_replace(code_labo, "HB", "P")
  data_SPIR_Ed$hyp_HLB[i] <- substr(code_ech, 1, 1)
  data_SPIR_Ed$HLB[i] <- substr(code_labo, 3, 3)
  data_SPIR_Ed$rep_arbre[i] <- (as.numeric(substr(code_labo, 3, 4)))
  data_SPIR_Ed$hyp_arbre[i] <- substr(code_labo, 4, 4)
  data_SPIR_Ed$rep_lot[i] <- (as.numeric(substr(code_labo, 5, 5)))
  data_SPIR_Ed$hyp_feuille[i] <- substr(code_labo, 10, 10)
  data_SPIR_Ed$rep_feuille[i] <- floor( (as.numeric(substr(code_labo, 8, 10))) /10) 
}
data_SPIR_Ed$X <-NULL
data_SPIR_Ed$HLB <- (as.numeric(mapvalues(data_SPIR_Ed$HLB, from=c(1, 2), to=c(1,0))))

# disponible pour le passage en format long jusqu'ici  !!! 

data_SPIR_Ed$code_ech_arbre <- factor(paste(data_SPIR_Ed$hyp_HLB, data_SPIR_Ed$hyp_arbre, sep=""))


data_SPIR_Ed$code_ech_lot <- factor(paste(data_SPIR_Ed$hyp_HLB, data_SPIR_Ed$hyp_arbre, data_SPIR_Ed$rep_lot, sep=""))

data_SPIR_Ed$code_ech_feuille <- factor(paste(data_SPIR_Ed$hyp_HLB, data_SPIR_Ed$hyp_arbre, data_SPIR_Ed$rep_lot, data_SPIR_Ed$hyp_feuille ,sep=""))
# length(levels(data_SPIR_Ed$code_ech)) correspondant à 2x HLB 3x lot et 7x arbres soit 42

#data_SPIR_Ed$hyp_HLB<-NULL

data_SPIR_Ed <- data_SPIR_Ed[,c(-(which(colnames(data_SPIR_Ed) == "hyp_HLB")),-(which(colnames(data_SPIR_Ed) == "hyp_arbre")),-(which(colnames(data_SPIR_Ed) == "hyp_feuille")))]

# disponible pour combiner avec les resultats Qpcr


# Importation et preparation des resultat de la Qpcr ####

data_Qpcr_Ed <- read.table(file = "Donnees/Donnees_Nathan/results_qPCR.csv"
                           
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA"))
data_Qpcr_Ed

trueP <- list(levels(factor(data_Qpcr_Ed$Sample.Name[which(data_Qpcr_Ed$C..Mean <32)])), levels(factor(data_Qpcr_Ed$Sample.Name[which(data_Qpcr_Ed$C..Mean <36)])) ) # on stocke les noms d'échantillons positifs selon nos 2 seuils de Ct

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
testHLB[[1]]$code_ech_lot <- paste(testHLB[[1]]$hyp_HLB, testHLB[[1]]$rep_arbre, testHLB[[1]]$rep_lot, sep="")
testHLB[[2]] <- testHLB[[1]]
testHLB[[3]] <- testHLB[[1]]
testHLB[[4]] <- testHLB[[1]]
names(testHLB) <- names(trueP)

for (l in 1:2) { #on applique la boucle à chaque élément "lot" de testHLB, soit on détermine les posi selon nos deux seuils de Ct
  for (i in 1:nrow(testHLB[[l]])) {
    if ( testHLB[[l]]$code_ech_lot[i] %in%  trueP[[l]]){ 
      testHLB[[l]]$qPCR_HLB[i] <- 1}
    else { testHLB[[l]]$qPCR_HLB[i] <- 0}
  }
  testHLB[[l]] <- lapply(testHLB[[l]], factor)
  testHLB[[l]] <- as.data.frame(testHLB[[l]])
}



for (l in 3:4) { #on applique la boucle à chaque élément de "arbre" testHLB , soit on détermine les posi selon nos deux seuils de Ct
  for (i in 1:nrow(testHLB[[l]])) {
    if ( substr(testHLB[[l]]$code_ech_lot[i], 1, 2) %in%  trueP[[l]]){ 
      testHLB[[l]]$qPCR_HLB[i] <- 1}
    else { testHLB[[l]]$qPCR_HLB[i] <- 0}
  }
  testHLB[[l]] <- lapply(testHLB[[l]], factor)
  testHLB[[l]] <- as.data.frame(testHLB[[l]])
}

# on obtient testHLB avec les positions des éléments positifs ou négatif au HLB


for (i in 1:length(testHLB)) { #on joint le résultat du test HLB au jdd de data_SPIR_Ed (on a 4 colonnes de test HLB: voir leurs noms)
  
  # !!!! Code_ech différent...Ecrire data_SPIR_ed avec code_ech de la Qpcr
  
  data_SPIR_Ed <- merge(data_SPIR_Ed, testHLB[[i]][4:5], by= "code_ech_lot", all.y= F, suffixes = c("", names(testHLB)[i]) )
}

data_SPIR_Ed <- data_SPIR_Ed[,c(1,2154:2163,3:2153)]
data_SPIR_Ed <- data_SPIR_Ed[,c(7,1,6,3:5,8:2162)]

# Il faut changer les noms car ftable impossible sinon...

names(data_SPIR_Ed )[match("qPCR_HLB",names(data_SPIR_Ed ))] <- "qPCR_HLBech_ct32_lot"

names(data_SPIR_Ed )[match("qPCR_HLBech_ct<36_lot",names(data_SPIR_Ed ))] <- "qPCR_HLBech_ct36_lot"

names(data_SPIR_Ed )[match("qPCR_HLBech_ct<32_arbre",names(data_SPIR_Ed ))] <- "qPCR_HLBech_ct32_arbre"

names(data_SPIR_Ed )[match("qPCR_HLBech_ct<36_arbre",names(data_SPIR_Ed ))] <- "qPCR_HLBech_ct36_arbre"

# CT<32_ARBRE ####

# Attribution  de la valeur 0 (negatif) ou 1 (positif) de la variable HLB en fonction des qPCR<32 ou qPCR<36

Masque_32_P <- selec[ ,(which(colnames(selec) == "qPCR_HLBech_ct32_arbre"))] == 1 # Masque logique 

Masque_32_N <- !Masque_32_P

selec_32_P <- selec[Masque_32_P,]

selec_32_N <- selec[Masque_32_N,]

ct32_arbre <- rbind.data.frame(selec_32_P,selec_32_N)   # fusion des jeu de donnee positif et negatif


# supression des variable inutiles

ct32_arbre <- ((ct32_arbre [, -(which(colnames(ct32_arbre) == "rep_arbre" | colnames(ct32_arbre) =="rep_lot" | colnames(ct32_arbre) =="rep_feuille" | colnames(ct32_arbre) =="HLB" | colnames(ct32_arbre) =="qPCR_HLBech_ct36_arbre"| colnames(ct32_arbre) =="qPCR_HLBech_ct32_lot"| colnames(ct32_arbre) =="qPCR_HLBech_ct36_lot"))])) 

ct32_arbre <- ct32_arbre[,c(3:2152,1:2)]

# Passage en format long de ct32_arbre

test_ct_32 <- pivot_longer( data =ct32_arbre, cols = 1:(ncol(ct32_arbre)-2) , values_to = "reflectance", names_to = "lambda"  )

test_ct_32 $lambda <- as.numeric(substring(test_ct_32 $lambda,2))

summary(test_ct_32)

length(test_ct_32$code_ech[test_ct_32$code_ech== "N51"])


# Calcul de la moyenne par arbre de ct32_arbre


Mean_mask_P<- ct32_arbre [ ,(which(colnames(ct32_arbre ) == "qPCR_HLBech_ct32_arbre"))] == 1 # nom des ligne positifs au HLB 

Mean_P_32 <- colMeans(ct32_arbre[Mean_mask_P]) # formule fausse ! car s'arret à 900 ligne la ou il y a le dernier TRUE

Mean_N_32 <- colMeans(ct32_arbre[!Mean_mask_P])  # negatif au HLB

Mean_HLB_32 <- rbind.data.frame(Mean_P_32,Mean_HLB_0)


# Passage en format long de la moyenne des lots

# Passage en format long de la moyenne des feuille


# arbre <- do.call(rbind, lapply(1:nrow(test_ed), function(x) {
#   data.frame(arbre= paste(unlist(test_ed[x,1:2]), collapse = ""))  # arbre correspond à l'ensemble (HLB + rep_arbre) combiné en 1 seul colonne
# }))

#test_ed <- cbind(test_ed, arbre)


summary(test_ed)



# entrainement et test sur data_spir_ed ------------------------------------------------------------------


## On construit le modèle sur le training set:

classifier = svm(y = data_SPIR_Ed$qPCR_HLBech_ct36_arbre , x=data_SPIR_Ed[10:2160], 
                 
                 type = 'C-classification', 
                 kernel = 'linear') 

## Et on le fait prédire sur le data_SPIR_Ed:
y_pred = predict(classifier, newdata = data_SPIR_Ed[c(10:2160)]) 


## On construit la matrice de confusion associée à cette prédiction:
cm = table(data_SPIR_Ed$qPCR_HLBech_ct36_arbre, y_pred) 
confusionMatrix(data_SPIR_Ed$qPCR_HLBech_ct36_arbre, y_pred)

## On vérifie qu'il n'y a pas de surapprentissage du modèle en réalisant une prédiction sur le data_SPIR_Ed:
y_pred_OF = predict(classifier, newdata = data_SPIR_Ed[c(10:2160)]) 
cmOF = table(data_SPIR_Ed$qPCR_HLBech_ct36_arbre, y_pred_OF)
cmOF #la prédiction n'est pas meilleure donc il n'y a pas de problème de surapprentissage

## K-Folds Cross-Validation: 
# On sépare aléatoirement le data_SPIR_Ed en 10 sous-jeu de données (ou folds), et on construit le modèle sur chacun d'eux
folds = createFolds(data_SPIR_Ed$qPCR_HLBech_ct36_arbre, k = 10)

## La fonction cv applique une fonction à chaque sous-jeud de données ou fold
cv = lapply(folds, function(x) { # début de la function
  # le data_SPIR_Ed est séparé en 10 grâce au lapply
  training_fold = data_SPIR_Ed[-x, ] 
  test_fold = data_SPIR_Ed[x, ] 
  # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
  classifier = svm(y = training_fold$qPCR_HLBech_ct36_arbre, x = training_fold[10:2160],
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'linear')
  
  # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
  y_pred = predict(classifier, newdata = test_fold[10:2160])
  cm = table(test_fold$qPCR_HLBech_ct36_arbre, y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

#on calcule la précision moyenne du modèle sur nos 10 folds:
accuracy_mean = mean(as.numeric(cv))
accuracy_mean

# accuracy  0.752381


# Test avec 1 seul rep de 10 feuilles ####

#ftable (rep_feuille + rep_lot ~ code_ech_arbre , data = data_SPIR_Ed) # Comment avoir 1 rep_feuille et 1 seul lot

selec_rep_lot_1 <- data_SPIR_Ed[,(which(colnames(data_SPIR_Ed) == "rep_lot"))]== 1 # Masque logique qui ne prend que les rep_lot =1

data_1rep_lot <- data_SPIR_Ed[selec_rep_lot_1,] #application du masque logique sur les LIGNES !!!

ftable (rep_feuille + rep_lot ~ code_ech_arbre , data = data_1rep_lot) # résultat

# Test fold aves 1 rep et 1 lot à ct<36

tirage <- do.call(rbind, lapply(1:nrow(data_1rep_lot), function(x) {
  data.frame(tirage= paste(unlist(data_1rep_lot[x, c((which(colnames(data_1rep_lot) == "qPCR_HLBech_ct36_arbre")),(which(colnames(data_1rep_lot) == "rep_arbre")):(which(colnames(data_1rep_lot) == "rep_feuille")))]), collapse = ""))  # tirage correspond à l'ensemble (qPCR_HLBech_ct36_arbre + rep_arbre + rep_lot + rep_feuille) combiné en 1 seul colonne
  
  
}))

data_1rep_lot <- cbind(data_1rep_lot, tirage)

Nmeas <- c(1:10)
data_global <- data.frame(accuracy= rep(NA,10), sensi=rep(NA,10), pospred=rep(NA,10))
#data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())

for (r in 1:10) {
  print(r)
  for (i in 1:length(Nmeas)) {
    if (i < 10) {
      set.seed(123) 
      split_five_meas = sample.split(data_1rep_lot$tirage , SplitRatio = Nmeas[i]/10) 
      five_meas <- subset(data_1rep_lot, split_five_meas == TRUE) 
      # on supprime colonne tirage: 
      five_meas <- five_meas[,-(which(colnames(five_meas) == "tirage"))]        #supression de la colonne tirage à la position ; which(colnames(five_meas) == "tirage")
    } 
    
    
    ## on réduit le nb de mesures à x mesures:
    
    
    #↓ on crée training et test set:
    set.seed(123) 
    split = sample.split(five_meas$qPCR_HLBech_ct36_arbre , SplitRatio = 0.75) 
    
    training_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct32_arbre"   ))]), split == TRUE) 
    #subset(five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" ))])
    test_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct32_arbre"   ))]), split == FALSE) 
    
    # marche pas car plusieurs tirage s'accumule...
    
    #centrage et réduction des données 
    
    # training_set  <- training_set[sapply(training_set, is.numeric)]  # tout doit être numerique
    #
    # View(sapply(training_set, is.numeric) )
    #
    #as.numeric(as.character(training_set))
    
    training_set[,-(which(colnames(training_set) == "qPCR_HLBech_ct36_arbre"))] <- scale(training_set[, -(which(colnames(training_set) == "qPCR_HLBech_ct36_arbre"))] )  # Il faut enlever tt ce qui n'est pas numerique dans five means ou enlever five_means...
    test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct36_arbre"))] <- scale(test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct36_arbre"))] ) 
    
    
    folds = createFolds(training_set$qPCR_HLBech_ct36_arbre, k = 10)
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    for (j in 1:length(folds)) {
      x= folds[[j]]
      training_fold = training_set[-x, ] 
      test_fold = training_set[x, ] 
      # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
      
      
      classifier = svm(y = training_fold$qPCR_HLBech_ct36_arbre, x = training_fold[,-(which(colnames(training_fold) == "qPCR_HLBech_ct36_arbre"))],
                       data = training_fold,
                       type = 'C-classification',
                       kernel = 'linear')
      
      #Marche jusque là !
      
      # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
      
      y_pred = predict(classifier, newdata = test_fold[,-(which(colnames(test_fold) == "qPCR_HLBech_ct36_arbre"))])
      test_fold$qPCR_HLBech_ct36_arbre <- factor(test_fold$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
      matrice <- confusionMatrix(test_fold$qPCR_HLBech_ct36_arbre, y_pred)
      accuracy <- accuracy + matrice$overall[1] /10
      sensi <- sensi + matrice$byClass[1] /10
      pospred <- pospred + matrice$byClass[3] /10
    }
    data_global[i,1] <- accuracy
    data_global[i,2] <- sensi
    data_global[i,3] <- pospred
  }
  #data_global_rep <- rbind(data_global_rep, data_global)
}

#write.csv2(data_global, "SVM_Juin_Nat.csv")

data_global <- data_global[-10,]  # suppression de la dernier ligne

data_global$nb <- seq(1:8)  # mettre les nombres d'observations

dataglobal2 <- melt(data_global,id.vars = "nb", measure.var=colnames(data_global[1:3]) ) # compile les variables

#save(list = ls(), file = "Sauvegardes_objet_R.data/Objets_SVM_Ed_CT36_1rep_1lot") 

# Graphique avec data_global transformé en dataglobal2

graph_1 <- ggplot(data=dataglobal2, aes(nb, value, fill= variable, colour= variable))+
  scale_x_continuous(breaks=seq(1:8))+
  geom_point()+
  # stat_summary(color = "black"
  #              , fun.data = function(x) mean_se(x, mult = qt(0.975,
  #                                                            length(x) - 1))
  #              , aes(shape = "IC à 95%")
  #              , position = position_dodge(width = 0.8)) +
  geom_line()+ labs(x = "Nombre de mesures de SPIR par feuille", y = "Performance de l'indicateur", title = "Qualité de la prédiction du statut HLB à CT<36_1rep & 1lot ", color = "Indicateur du test SVM") +
  theme_gray() 

graph_1


#esquisse::esquisser() 

# ggsave(filename = "Graphiques/Parametres_SVM_SPIR_Ed_ct36_1rep_1lot.png", plot = last_plot())

# 


# 

# R markdown poubelle ####


tmb=split(data_SPIR_Ed,(data_SPIR_Ed$code_ech_feuille)) 



l_c = 100    # nombre de creation training_set et test_set 
l_a = 6      # nombdre de rep par feuille (sur un total de 6 repetition)
l_Z = 420    # nombre  de feuilles (sur un total de 420 feuilles)



data_global <- data.frame(accuracy= rep(NA,l_c), sensi=rep(NA,l_c), pospred=rep(NA,l_c))
var_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
sd_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
mean_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c())


progress = function(a, length_a){
  print(a/length_a)
}

r = 0
stock = 0



for (a in 1:l_a) {
  
  tmp <- sample(x = tmb,size = l_Z ,replace = FALSE,prob = NULL)
  
  maliste =lapply(1:length(tmp),function(b){ 
    tmp[[b]][sample(1:length(tmp[[b]]$code_ech_feuille),a),]   
  })
  
  test_ed=do.call(rbind,maliste) 
  
  for (c in 1:l_c ){
    if (c == 1){ T_1 = Sys.time()}
    if (c == 2){ T_2 = Sys.time();
    diff = T_2 - T_1
    
    for (j in 1:l_a){r = j+r}
    temp_total = as.numeric(diff) * l_c * r
    ttm = temp_total/3600
    hours = floor(ttm)
    mt = (ttm - hours)*60
    minutes = floor(mt)
    st = (mt-minutes)*60 ; secondes = floor(st)
    print(paste("Il reste", hours, "heures", minutes, "minutes", secondes, "secondes", sep = " "))
    }
    
    
    print_a = paste("Le programme en est à ", substr(stock + 100*(c-1)/l_c/(r-a+1), start = 1, stop = 4) , "% de sa progression", sep = "")
    if (c == l_c){stock = 100*(c-1)/l_c/(r-a+1)+stock}
    print(print_a)
    
    #set.seed(123)
    
    split = sample.split(test_ed$qPCR_HLBech_ct36_arbre, SplitRatio = 0.75) 
    
    training_set <- subset(test_ed, split == TRUE) 
    test_set <- subset(test_ed, split == FALSE) 
    
    
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    
    classifier = svm(y = training_set$qPCR_HLBech_ct36_arbre, x = training_set[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))],
                     data = training_set,
                     type = 'C-classification',
                     kernel = 'linear')
    
    
    
    y_pred = predict(classifier, newdata = test_set[,c ((which(colnames(test_set) == "X350" )): (which(colnames(test_set) == "X2500" )))])
    test_set$qPCR_HLBech_ct36_arbre <- factor(test_set$qPCR_HLBech_ct36_arbre)           
    matrice <- confusionMatrix(test_set$qPCR_HLBech_ct36_arbre, y_pred)
    
    
    accuracy <- accuracy + matrice$overall[1] 
    sensi <- sensi + matrice$byClass[1] 
    pospred <- pospred + matrice$byClass[3] 
    
    data_global[c,1] <- accuracy
    data_global[c,2] <- sensi
    data_global[c,3] <- pospred
  }
  
  sd_rep  <-   rbind(sd_rep,c(sd(data_global$accuracy), sd(data_global$pospred), sd(data_global$sensi)))   
  mean_rep <- rbind(mean_rep,colMeans(data_global))
  print(a) # ce qui est renvoye a l'ecran et correspond au nombre de rep_feuille
  
} # fermeture a





names(sd_rep)[1:3] <- c("a_sd","s_sd","p_sd")

names(mean_rep)[1:3] <- c("Accuracy","Sensitivity","Precision")




mean_rep$nb_rep <- seq(nrow(mean_rep)) 

mean_rep_long <- melt(mean_rep,id.vars = "nb_rep", measure.var=colnames(mean_rep[1:3]))

sd_rep$ns <- seq(nrow(sd_rep))

sd_rep_long <- melt(sd_rep,id.vars = "ns", measure.var=colnames(sd_rep[1:3]))


data_global_rep <- cbind(mean_rep_long ,sd_rep_long )

data_global_rep<- data_global_rep[,-(which(colnames(data_global_rep) == "ns" ))]

names(data_global_rep)[2:5] <- c("Mean_values","Mean","Sd_values", "Sd")




IC <- (1.96*data_global_rep$Sd)/sqrt(l_c)

# Importation des R.data car le temps de calcul est long

load("C:/Users/esori/Desktop/Stage La reunion/R/Analyse_Spectro_R/Sauvegardes_objet_R.data/SVM_ct36_1lot_6rep_100byRep_Ed_ss_test_fold_420f")

# Graphique des paramètre de SVM

graph_rep_420f <- ggplot(data=data_global_rep, aes(nb_rep, Mean, colour= Mean_values))+
  scale_x_continuous(breaks=seq(1:6))+
  geom_errorbar(aes(ymin = Mean-IC , ymax = Mean+IC),width=0.05, lwd = 1.1 )+
  geom_point(size = 3)+
  geom_line(lwd = 1.1) +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction du statut HLB  des abres par SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM pour chaque répétition de feuille (420f)", color = "Paramètres de robustesse en SVM ") +
  
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


graph_rep_420f



# SVM : TEST FOLD CT<32####

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) \n et à la validation du modèle (test_set):

tirage <- do.call(rbind, lapply(1:nrow(data_SPIR_Ed), function(x) {
  data.frame(tirage= paste(unlist(data_SPIR_Ed[x, c((which(colnames(data_SPIR_Ed) == "qPCR_HLBech_ct32_arbre")),(which(colnames(data_SPIR_Ed) == "rep_arbre")):(which(colnames(data_SPIR_Ed) == "rep_feuille")))]), collapse = ""))  # tirage correspond à l'ensemble (qPCR_HLBech_ct32_arbre + rep_arbre + rep_lot + rep_feuille) combiné en 1 seul colonne
  
  
}))

data_SPIR_Ed <- cbind(data_SPIR_Ed, tirage)

Nmeas <- c(1:10)
data_global <- data.frame(accuracy= rep(NA,10), sensi=rep(NA,10), pospred=rep(NA,10))
#data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())

for (r in 1:10) {
  print(r)
  for (i in 1:length(Nmeas)) {
    if (i < 10) {
      set.seed(123) 
      split_five_meas = sample.split(data_SPIR_Ed$tirage , SplitRatio = Nmeas[i]/10) 
      five_meas <- subset(data_SPIR_Ed, split_five_meas == TRUE) 
      # on supprime colonne tirage: 
      five_meas <- five_meas[,-(which(colnames(five_meas) == "tirage"))]        #supression de la colonne tirage à la position ; which(colnames(five_meas) == "tirage")
    } 
    
    
    ## on réduit le nb de mesures à x mesures:
    
    
    #↓ on crée training et test set:
    set.seed(123) 
    split = sample.split(five_meas$qPCR_HLBech_ct32_arbre , SplitRatio = 0.75) 
    
    training_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_arbre"   ))]), split == TRUE) 
    #subset(five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" ))])
    test_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_arbre"   ))]), split == FALSE) 
    
    # marche pas car plusieurs tirage s'accumule...
    
    #centrage et réduction des données 
    
    # training_set  <- training_set[sapply(training_set, is.numeric)]  # tout doit être numerique
    #
    # View(sapply(training_set, is.numeric) )
    #
    #as.numeric(as.character(training_set))
    
    training_set[,-(which(colnames(training_set) == "qPCR_HLBech_ct32_arbre"))] <- scale(training_set[, -(which(colnames(training_set) == "qPCR_HLBech_ct32_arbre"))] )  # Il faut enlever tt ce qui n'est pas numerique dans five means ou enlever five_means...
    test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct32_arbre"))] <- scale(test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct32_arbre"))] ) 
    
    
    folds = createFolds(training_set$qPCR_HLBech_ct32_arbre, k = 10)
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    for (j in 1:length(folds)) {
      x= folds[[j]]
      training_fold = training_set[-x, ] 
      test_fold = training_set[x, ] 
      # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
      
      
      classifier = svm(y = training_fold$qPCR_HLBech_ct32_arbre, x = training_fold[,-(which(colnames(training_fold) == "qPCR_HLBech_ct32_arbre"))],
                       data = training_fold,
                       type = 'C-classification',
                       kernel = 'linear')
      
      #Marche jusque là !
      
      # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct32_arbre en utilisant le test_FOLD et on calcule la précision du modèle
      
      y_pred = predict(classifier, newdata = test_fold[,-(which(colnames(test_fold) == "qPCR_HLBech_ct32_arbre"))])
      test_fold$qPCR_HLBech_ct32_arbre <- factor(test_fold$qPCR_HLBech_ct32_arbre)            # pour mettre test_fold$qPCR_HLBech_ct32_arbre en facteur pour la fonction cufusionMatrix
      matrice <- confusionMatrix(test_fold$qPCR_HLBech_ct32_arbre, y_pred)
      accuracy <- accuracy + matrice$overall[1] /10
      sensi <- sensi + matrice$byClass[1] /10
      pospred <- pospred + matrice$byClass[3] /10
    }
    data_global[i,1] <- accuracy
    data_global[i,2] <- sensi
    data_global[i,3] <- pospred
  }
  #data_global_rep <- rbind(data_global_rep, data_global)
}

#write.csv2(data_global, "SVM_Juin_Nat.csv")

data_global <- data_global[-10,]  # suppression de la dernier ligne

data_global$nb <- seq(1:9)  # mettre les nombres d'observations

dataglobal2 <- melt(data_global,id.vars = "nb", measure.var=colnames(data_global[1:3]) ) # compile les variables

#save(list = ls(), file = "Sauvegardes_objet_R.data/Objets_SVM_Ed_CT32") 

# Graphique avec data_global transformé en dataglobal2

graph_1 <- ggplot(data=dataglobal2, aes(nb, value, fill= variable, colour= variable))+
  scale_x_continuous(breaks=seq(1:9))+
  geom_point()+
  # stat_summary(color = "black"
  #              , fun.data = function(x) mean_se(x, mult = qt(0.975,
  #                                                            length(x) - 1))
  #              , aes(shape = "IC à 95%")
  #              , position = position_dodge(width = 0.8)) +
  geom_line()+ labs(x = "Nombre de mesures de SPIR par feuille", y = "Performance de l'indicateur", title = "Qualité de la prédiction du statut HLB à CT<32 ", color = "Indicateur du test SVM") +
  theme_gray() 

graph_1


#esquisse::esquisser() 

# ggsave(filename = "Graphiques/Parametres_SVM_SPIR_Ed_ct32.png", plot = last_plot())

# 

# SVM : TEST FOLD CT<36####

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) \n et à la validation du modèle (test_set):



tirage <- do.call(rbind, lapply(1:nrow(data_SPIR_Ed), function(x) {
  data.frame(tirage= paste(unlist(data_SPIR_Ed[x, c((which(colnames(data_SPIR_Ed) == "qPCR_HLBech_ct36_arbre")),(which(colnames(data_SPIR_Ed) == "rep_arbre")):(which(colnames(data_SPIR_Ed) == "rep_feuille")))]), collapse = ""))  # tirage correspond à l'ensemble (qPCR_HLBech_ct36_arbre + rep_arbre + rep_lot + rep_feuille) combiné en 1 seul colonne
  
  
}))

data_SPIR_Ed <- cbind(data_SPIR_Ed, tirage)

Nmeas <- c(1:10)
data_global <- data.frame(accuracy= rep(NA,10), sensi=rep(NA,10), pospred=rep(NA,10))
#data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())

#for (r in 1:10)  r de 1 a 3 sinon ça met trop longtemp

for (r in 1:3) {
  print(r) # ce qui est renvoye a l'ecran
  for (i in 1:length(Nmeas)) {
    if (i < 10) {
      set.seed(123) 
      split_five_meas = sample.split(data_SPIR_Ed$tirage , SplitRatio = Nmeas[i]/10) # facteur aleatoire ?
      five_meas <- subset(data_SPIR_Ed, split_five_meas == TRUE) 
      # on supprime colonne tirage: 
      five_meas <- five_meas[,-(which(colnames(five_meas) == "tirage"))]        #supression de la colonne tirage à la position ; which(colnames(five_meas) == "tirage")
    } 
    
    
    ## on réduit le nb de mesures à x mesures:
    
    
    #↓ on crée training et test set:
    set.seed(123) 
    split = sample.split(five_meas$qPCR_HLBech_ct36_arbre , SplitRatio = 0.75) # facteur aleatoire ?
    
    training_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct32_arbre"   ))]), split == TRUE) 
    #subset(five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" ))])
    test_set <- subset((five_meas [, -(which(colnames(five_meas) == "rep_arbre" | colnames(five_meas) =="rep_lot" | colnames(five_meas) =="rep_feuille" | colnames(five_meas) =="code_ech_lot" | colnames(five_meas) =="code_ech_arbre" | colnames(five_meas) =="qPCR_HLBech_ct32_lot" | colnames(five_meas) =="qPCR_HLBech_ct36_lot" | colnames(five_meas) =="qPCR_HLBech_ct32_arbre"   ))]), split == FALSE) 
    
    # marche pas car plusieurs tirage s'accumule...
    
    #centrage et réduction des données 
    
    # training_set  <- training_set[sapply(training_set, is.numeric)]  # tout doit être numerique
    #
    # View(sapply(training_set, is.numeric) )
    #
    #as.numeric(as.character(training_set))
    
    training_set[,-(which(colnames(training_set) == "qPCR_HLBech_ct36_arbre"))] <- scale(training_set[, -(which(colnames(training_set) == "qPCR_HLBech_ct36_arbre"))] )  # Il faut enlever tt ce qui n'est pas numerique dans five means ou enlever five_means...
    test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct36_arbre"))] <- scale(test_set[, -(which(colnames(test_set) == "qPCR_HLBech_ct36_arbre"))] ) 
    
    
    folds = createFolds(training_set$qPCR_HLBech_ct36_arbre, k = 10)# On sépare aléatoirement le training_set en 10 sous-jeu de données (ou folds), et on construit le modèle sur chacun d'eux
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    for (j in 1:length(folds)) {
      x= folds[[j]]
      training_fold = training_set[-x, ] 
      test_fold = training_set[x, ] 
      # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
      
      
      classifier = svm(y = training_fold$qPCR_HLBech_ct36_arbre, x = training_fold[,-(which(colnames(training_fold) == "qPCR_HLBech_ct36_arbre"))],
                       data = training_fold,
                       type = 'C-classification',
                       kernel = 'linear')
      
      #Marche jusque là !
      
      # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
      
      y_pred = predict(classifier, newdata = test_fold[,-(which(colnames(test_fold) == "qPCR_HLBech_ct36_arbre"))])
      test_fold$qPCR_HLBech_ct36_arbre <- factor(test_fold$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
      matrice <- confusionMatrix(test_fold$qPCR_HLBech_ct36_arbre, y_pred)
      return(matrice)
      accuracy <- accuracy + matrice$overall[1] /10
      sensi <- sensi + matrice$byClass[1] /10
      pospred <- pospred + matrice$byClass[3] /10
    }
    data_global[i,1] <- accuracy
    data_global[i,2] <- sensi
    data_global[i,3] <- pospred
  }
  #data_global_rep <- rbind(data_global_rep, data_global)
}

#write.csv2(data_global, "SVM_Juin_Nat.csv")

data_global <- data_global[-10,]  # suppression de la dernier ligne

data_global$nb <- seq(1:9)  # mettre les nombres d'observations

dataglobal2 <- melt(data_global,id.vars = "nb", measure.var=colnames(data_global[1:3]) ) # compile les variables

#save(list = ls(), file = "Sauvegardes_objet_R.data/Objets_SVM_Ed") 

# Graphique avec data_global transformé en dataglobal2

graph_1 <- ggplot(data=dataglobal2, aes(nb, value, fill= variable, colour= variable))+
  scale_x_continuous(breaks=seq(1:9))+
  geom_point()+
  # stat_summary(color = "black"
  #              , fun.data = function(x) mean_se(x, mult = qt(0.975,
  #                                                            length(x) - 1))
  #              , aes(shape = "IC à 95%")
  #              , position = position_dodge(width = 0.8)) +
  geom_line()+ labs(x = "Nombre de mesures de SPIR par feuille", y = "Performance de l'indicateur", title = "Qualité de la prédiction du statut HLB à CT<36 ", color = "Indicateur du test SVM") +
  theme_gray() 

graph_1


#esquisse::esquisser() 

# ggsave(filename = "Graphiques/Parametres_SVM_SPIR_Ed_ct36.png", plot = last_plot())

# 

# SVM Kernel linéaire ct36 ####

# On divise le jeu de données en 2 sous-jeux servant à l'apprentissage (training_set) \n et à la validation du modèle (test_set):
set.seed(123) 
split = sample.split(data_SPIR_Ed$qPCR_HLBech_ct36_lot, SplitRatio = 0.75) 

training_set <- subset(data_SPIR_Ed, split == TRUE) 
test_set <- subset(data_SPIR_Ed, split == FALSE) 

## On construit le modèle sur le training set:
classifier = svm(y = training_set$qPCR_HLBech_ct36_arbre , x=training_set[10:1890], 
                 
                 type = 'C-classification', 
                 kernel = 'linear') 

## Et on le fait prédire sur le test_set:

rownames(test_set) <- sapply(1:length(test_set$code_ech_arbre), function(i){paste(test_set$code_ech_arbre[i],i,sep=".")})

y_pred = predict(classifier, newdata = test_set[c(10:1890)]) 


## On construit la matrice de confusion associée à cette prédiction:
cm = table(test_set$qPCR_HLBech_ct36_arbre, y_pred) 
confusionMatrix(test_set$qPCR_HLBech_ct36_arbre, y_pred)

## On vérifie qu'il n'y a pas de surapprentissage du modèle en réalisant une prédiction sur le training_set:
y_pred_OF = predict(classifier, newdata = training_set[c(10:1890)]) 
cmOF = table(training_set$qPCR_HLBech_ct36_arbre, y_pred_OF)
cmOF #la prédiction n'est pas meilleure donc il n'y a pas de problème de surapprentissage

## K-Folds Cross-Validation: 
# On sépare aléatoirement le training_set en 10 sous-jeu de données (ou folds), et on construit le modèle sur chacun d'eux
folds = createFolds(training_set$qPCR_HLBech_ct36_arbre, k = 10)

## La fonction cv applique une fonction à chaque sous-jeud de données ou fold
cv = lapply(folds, function(x) { # début de la function
  # le training_set est séparé en 10 grâce au lapply
  training_fold = training_set[-x, ] 
  test_fold = training_set[x, ] 
  # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
  classifier = svm(y = training_fold$qPCR_HLBech_ct36_arbre, x = training_fold[10:1890],
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'linear')
  
  # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
  y_pred = predict(classifier, newdata = test_fold[10:1890])
  cm = table(test_fold$qPCR_HLBech_ct36_arbre, y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

#on calcule la précision moyenne du modèle sur nos 10 folds:
accuracy = mean(as.numeric(cv))
accuracy


