
rm(list=ls())

# LIBRARY ####

library(ggplot2)
library(ggparty)
library(partykit)
library(cowplot)
library(gridExtra)
library(ggdark)
library(esquisse)
library(colourpicker)
library(readxl)
library(imager)  
library(data.table)
library(plyr)
library(tidyr)
library(tidyverse)  
library(caret) # package Machine learning
library(pls)     
library(caTools) 
library(reshape2)
library(e1071) # SVM
library(gt) 
library(ROCR)
library(MASS) # package LDA
library(mda)
library(klaR)
library(dplyr)
library(rminer)
library(party)  # Random Forest
library(plsdepot) # Partial Least Square

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

# Exploration du jeu de donnees     ####


ftable (rep_lot ~ rep_arbre , data = data_SPIR_Ed)  # repartition en fonction des lots

ftable (rep_feuille ~ rep_arbre , data = data_SPIR_Ed) # repartition en fonction des feuilles

ftable (rep_feuille + rep_lot ~ code_ech_arbre , data = data_SPIR_Ed) # repartition en fonction des feuilles,lot et arbres

View(ftable (code_ech_feuille ~ rep_feuille , data = data_SPIR_Ed) )


length(levels(data_SPIR_Ed$code_ech_feuille))# le code est bon car on a bien 420 feuilles echantilloné en tt

View (ftable (code_ech_feuille  ~ qPCR_HLBech_ct32_lot , data = data_SPIR_Ed) )

# Exploration avec les resultats qPCR ####



# repartition des lots en fonction des resultat de la qPCR à ct<32_lot

ftable (rep_lot  ~ qPCR_HLBech_ct32_lot + rep_arbre   , data = data_SPIR_Ed) 
ftable (rep_arbre  ~ qPCR_HLBech_ct32_lot , data = data_SPIR_Ed) 


# repartition des lots en fonction des resultat de la qPCR à ct<32_arbre

ftable (code_ech_arbre  ~ qPCR_HLBech_ct32_arbre , data = data_SPIR_Ed) 


# repartition des lots en fonction des resultat de la qPCR à ct<36_lot

ftable (rep_lot  ~ qPCR_HLBech_ct36_lot + rep_arbre   , data = data_SPIR_Ed) 
ftable (rep_arbre  ~ qPCR_HLBech_ct36_lot , data = data_SPIR_Ed) 

# repartition des lots en fonction des resultat de la qPCR à ct<36_arbre

ftable (code_ech_arbre  ~ qPCR_HLBech_ct36_arbre , data = data_SPIR_Ed) 

#save(list = ls(), file = "Sauvegardes_objet_R.data/SPIR_qPCR_ED") 

#-> Passage en format long ####

# importation des R.data = SPIR_qPCR_ED #save(list = ls(), file = "Sauvegardes_objet_R.data/SPIR_qPCR_ED") 

selec <- data_SPIR_Ed  


selec <- selec[,c(which(colnames(selec) == "X350"):which(colnames(selec) == "X2500"),which(colnames(selec) == "code_ech_feuille"):which(colnames(selec) == "code_ech_arbre"),which(colnames(selec) == "qPCR_HLBech_ct32_lot"):which(colnames(selec) == "qPCR_HLBech_ct36_arbre"))]  # on enleve tt les rep_...

test_ct <- pivot_longer( data =selec, cols = 1:(ncol(selec)-7) , values_to = "reflectance", names_to = "lambda"  )

test_ct $lambda <- as.numeric(substring(test_ct $lambda,2))

summary(test_ct)

# Nouveau tableau avec cette fois ci seulement les longueurs d'onde de 400 à 680 nm 

selec_400_680 <- selec[,c((which(colnames(selec) == "X400")):(which(colnames(selec) == "X680")),(which(colnames(selec) == "code_ech_lot")):(which(colnames(selec) == "qPCR_HLBech_ct36_arbre")))]

test_400_680 <- pivot_longer( data =selec_400_680, cols = 1:(ncol(selec_400_680)-6) , values_to = "reflectance", names_to = "lambda"  )

test_400_680 $lambda <- as.numeric(substring(test_400_680 $lambda,2))

# format long avec moyennes pour chaque arbre en ct_32 + graph ####

test_ct_32 <- aggregate( test_ct$reflectance ,
                         list(test_ct$qPCR_HLBech_ct32_arbre , test_ct$code_ech_arbre, test_ct$lambda)   , mean) 

names(test_ct_32) = c("ct<32","code_ech_arbre","lambda","reflectance")

ggplot(test_ct_32) +
  aes(x = lambda, y = reflectance, colour = code_ech_arbre, size = `ct<32`) +
  geom_point() +
  scale_color_hue() +
  theme_minimal()

#esquisse::esquisser() 

# just avec la moyenne des courbes positives et des coubres negatives au HLB à ct<32

test_HLB_32 <- aggregate( test_ct$reflectance ,
                          list(test_ct$qPCR_HLBech_ct32_arbre , test_ct$lambda)   , mean) 

names(test_HLB_32) = c("HLB_32","lambda","reflectance")

ggplot(test_HLB_32) +
  aes(x = lambda, y = reflectance, colour = HLB_32) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est  positif (Ct<32)", color = "Resultat du test HLB à Ct<32") +
  theme_gray() 

# graphique zoomé sur les longueurs d'ondes 400 à 680 nm a ct<32

test_zoom_32 <-  aggregate( test_400_680$reflectance ,
                            list(test_400_680$qPCR_HLBech_ct32_arbre , test_400_680$lambda)   , mean) 

names(test_zoom_32) = c("HLB_32","lambda","reflectance")

ggplot(test_zoom_32) +
  aes(x = lambda, y = reflectance, colour = HLB_32) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est  positif (Ct<32)", color = "Resultat du test HLB à Ct<32") +
  theme_gray() 

# format long avec moyennes pour chaque arbre en ct_36 + graph ####

test_ct_36 <- aggregate( test_ct$reflectance ,
                         list(test_ct$qPCR_HLBech_ct36_arbre , test_ct$code_ech_arbre, test_ct$lambda)   , mean) 

names(test_ct_36) = c("ct<36","code_ech_arbre","lambda","reflectance")

ggplot(test_ct_36) +
  aes(x = lambda, y = reflectance, colour = code_ech_arbre, size = `ct<36`) +
  geom_point() +
  scale_color_hue() +
  theme_minimal()

# just avec la moyenne des courbes positives et des coubres negatives au HLB à ct<36

test_HLB_36 <- aggregate( test_ct$reflectance ,
                          list(test_ct$qPCR_HLBech_ct36_arbre , test_ct$lambda)   , mean) 

names(test_HLB_36) = c("HLB_36","lambda","reflectance_moyenne")

ggplot(test_HLB_36) +
  aes(x = lambda, y = reflectance_moyenne, colour = HLB_36) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est  positif (Ct<36)", color = "Resultat du test HLB à Ct<36") +
  #dark_theme_gray() 
  theme_gray() 

# graphique zoomé sur les longueurs d'ondes 400 à 680 nm a ct<36

test_zoom_36 <-  aggregate( test_400_680$reflectance ,
                            list(test_400_680$qPCR_HLBech_ct36_arbre , test_400_680$lambda)   , mean) 

names(test_zoom_36) = c("HLB_36","lambda","reflectance")

ggplot(test_zoom_36) +
  aes(x = lambda, y = reflectance, colour = HLB_36) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "Un arbre est considéré comme positif dès qu'un lot sur les 3 est  positif (Ct<36)", color = "Resultat du test HLB à Ct<36") +
  theme_gray() 

#ggsave(filename = "Graphiques/Spectre 400_680 moyen en fonction du statut HLB des arbres à ct36.png", plot = last_plot())

# Calcul de la Discriminabilité  via D = (Mean0-Mean1)/((Sd1 + Sd0)/2) ; 0 : Negatif au HLB , 1 : Positif au HLB ####

# Passage en format long

selec <- data_SPIR_Ed  


selec <- selec[,c(which(colnames(selec) == "X350"):which(colnames(selec) == "X2500"),which(colnames(selec) == "code_ech_feuille"):which(colnames(selec) == "code_ech_arbre"),which(colnames(selec) == "qPCR_HLBech_ct32_lot"):which(colnames(selec) == "qPCR_HLBech_ct36_arbre"))]  # on enleve tt les rep_...

test_ct <- pivot_longer( data =selec, cols = 1:(ncol(selec)-7) , values_to = "reflectance", names_to = "lambda"  )

test_ct $lambda <- as.numeric(substring(test_ct $lambda,2))

# Moyenne des spectres en fonction des lambda et du statu HLB

test_HLB_36 <- aggregate( test_ct$reflectance ,
                          list(test_ct$qPCR_HLBech_ct36_arbre , test_ct$lambda)   , mean) 

names(test_HLB_36) = c("HLB_36","lambda","reflectance_moyenne")

# Ecart type des spectres en fonction des lambda et du statu HLB

test_sd_36 <- aggregate( test_ct$reflectance ,
                         list(test_ct$qPCR_HLBech_ct36_arbre , test_ct$lambda)   , sd) 

names(test_sd_36) = c("HLB_36","lambda","reflectance_Sd")

# Fusion des deux jeu moyenne et sd

test_ref <- cbind(test_HLB_36,test_sd_36[,3])

names(test_ref) = c("HLB_36","lambda","reflectance_moyenne","reflectance_Sd")

# Calcul de discriminabilité

test_discrim <-  test_ref[test_ref$HLB_36 == 0,]

test_discrim$refelectance_m1 <-  test_ref$reflectance_moyenne[test_ref$HLB_36 == 1]

test_discrim$refelectance_sd1 <-  test_ref$reflectance_Sd[test_ref$HLB_36 == 1]

names(test_discrim) = c("HLB_36","lambda","reflectance_m0","reflectance_sd0","reflectance_m1","reflectance_sd1")

test_discrim$discrim <-  (test_discrim$reflectance_m1 - test_discrim$reflectance_m0) / ((test_discrim$reflectance_sd1 + test_discrim$reflectance_sd0)/2)

test_discrim$Sratio <- test_discrim$reflectance_m1 / test_discrim$reflectance_m0

test_discrim$diff <- test_discrim$reflectance_m1 - test_discrim$reflectance_m0



# graph_discrim

graph_discrim <- ggplot(test_discrim) +
  
  geom_line(aes(x = lambda , y = discrim)) +
  
  geom_hline(yintercept = 0.45344 , linetype = "dotted" , size = 0.1 , color = "red" ) +
  
  geom_point(aes(x = 1476 , y = 0.45344 ), shape = 3 , fill="darkred" ,color = "darkred" , size = 6) +
  
  annotate("text", label = "X1476", x = 1476, y = 0.473, size = 4, colour = "red") +
  
  geom_point(aes(x = 712 , y = 0.2216931 ), shape = 3 , fill="darkred" ,color = "darkred" , size = 6)+
  
  annotate("text", label = "X712", x = 735, y = 0.24, size = 4, colour = "red") +
  
  geom_point(aes(x = 465 , y = 0.3235764 ), shape = 3 , fill="darkred" ,color = "darkred" , size = 6)+
  
  annotate("text", label = "X465", x = 465, y = 0.34, size = 4, colour = "red") +
  
  labs(x = "Longueur d'onde (en nm)", y = "Valeur de discriminabilité", title = "Valeurs des discriminabilités", subtitle = "Qpcr à Ct36", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20"))

graph_discrim


#dev.off()

# graph_Sratio

graph_Sratio <- ggplot(test_discrim) +
  
  geom_line(aes(x = lambda , y = Sratio)) +
  
  geom_hline(yintercept = 1.1464 , linetype = "dotted" , size = 0.1 , color = "red" ) +
  
  geom_point(aes(x = 633 , y = 1.1464  ), shape = 3 , fill="darkred" ,color = "darkred" , size = 6) +
  
  annotate("text", label = "X633", x = 633, y = 1.151, size = 4, colour = "red") +
  
  labs(x = "Longueur d'onde (en nm)", y = "Ratio Spectral", title = "Valeurs des ratios spectrals", subtitle = "Qpcr à Ct36", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20"))

graph_Sratio

# graph_diff

graph_diff <- ggplot(test_discrim) +
  
  geom_line(aes(x = lambda , y = diff)) +
  
  geom_hline(yintercept = 0.024064 , linetype = "dotted" , size = 0.1 , color = "red" ) +
  
  geom_point(aes(x = 718 , y = 0.024064  ), shape = 3 , fill="darkred" ,color = "darkred" , size = 6) +
  
  annotate("text", label = "X718", x = 718, y = 0.0248, size = 4, colour = "red") +
  
  labs(x = "Longueur d'onde (en nm)", y = "Difference m1-m0", title = "Valeurs des différences entre les réflectences moyennes", subtitle = "Qpcr à Ct36", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20"))

graph_diff

# Groupement des 3 graphs

# graph_grid <- grid.arrange(graph_discrim, graph_diff,graph_Sratio, nrow = 1)
# 
# graph_grid

grid_diff <- plot_grid(graph_discrim , graph_diff , graph_Sratio , ncol = 1)

grid_diff

#png("graph_grid.png")

#ggsave(filename = "Graphiques/Discriminabilite/Graphiques groupés des mesures de disciminabilité.png", plot = last_plot())

# Exploration

summary(test_discrim)

# Max_discrim

Masque_discrim <- test_discrim[,(which(colnames(test_discrim) == "discrim"))] > 0.453  # Masque logique

Max_discrim <- test_discrim[Masque_discrim,]

# Max_Sratio

Masque_Sratio <- test_discrim[,(which(colnames(test_discrim) == "Sratio"))] > 1.1464   # Masque logique

Max_Sratio <- test_discrim[Masque_Sratio,]

# Max_diff

Masque_diff <- test_discrim[,(which(colnames(test_discrim) == "diff"))] > 0.024064  # Masque logique

Max_diff <- test_discrim[Masque_diff,]



#ggsave(filename = "Graphiques/Discriminabilite/Valeurs des différences entre les réflectences moyennes.png", plot = last_plot())


#write.table(x =test_discrim  , file = "Donnees/test_discrim.xls" ,sep = ",")


# Extraction des paramètre SVM du statut HLB à CT36 ####

tmp=split(data_SPIR_Ed,(data_SPIR_Ed$code_ech_feuille)) 
# maliste=lapply(1:length(tmp),function(i){ 
#   tmp[[i]][sample(1:length(tmp[[i]]$code_ech_feuille),2),]   # on fait ça pour tte les feuilles mais tjrs en choisissant 2 rep tire aleatoirement
# })
# 
# test_ed=do.call(rbind,maliste) 

# set.seed(123) # fix la valeur aleatoire du tirage, mettre en commentaire pour avoir un tirage aleatoire
# 
# split = sample.split(test_ed$qPCR_HLBech_ct36_lot, SplitRatio = 0.75) 
# 
# training_set <- subset(test_ed, split == TRUE) 
# test_set <- subset(test_ed, split == FALSE) 

l_c = 10  # nombre de creation training_set et test_set 
l_a = 6   # nombdre de rep par feuille feuille

# Verifier que la taille de data_global corresponde à la taille de l_c avant chaque lancement !!!!!

data_global <- data.frame(accuracy= rep(NA,l_c), sensi=rep(NA,l_c), pospred=rep(NA,l_c))
var_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
sd_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
mean_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c())
#data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())





progress = function(a, length_a){
  print(a/length_a)
}

r = 0
stock = 0

for (a in 1:l_a) {
  
  
  maliste =lapply(1:length(tmp),function(b){ 
    tmp[[b]][sample(1:length(tmp[[b]]$code_ech_feuille),a),]   # on fait ça pour tte les feuilles mais tjrs en choisissant 2 rep tire aleatoirement
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
    
    #folds = createFolds(training_set$qPCR_HLBech_ct36_arbre, k = 10)# On sépare aléatoirement le training_set en 10 sous-jeu de données (ou folds), et on construit le modèle sur chacun d'eux
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    #for (d in 1:length(folds)) {
    #x= folds[[d]]
    #training_fold = training_set[-x, ] 
    #test_fold = training_set[x, ] 
    # On fait apprendre le modèle sur chaque training_FOLD, le modèle est stocké dans classifier 
    
    classifier = svm(y = training_set$qPCR_HLBech_ct36_arbre, x = training_set[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))],
                     data = training_set,
                     type = 'C-classification',
                     kernel = 'linear')
    
    # classifier = svm(y = training_fold$qPCR_HLBech_ct36_arbre, x = training_fold[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))],
    #                  data = training_fold,
    #                  type = 'C-classification',
    #                  kernel = 'linear')
    
    #Marche jusque là !
    
    # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
    
    y_pred = predict(classifier, newdata = test_set[,c ((which(colnames(test_set) == "X350" )): (which(colnames(test_set) == "X2500" )))])
    test_set$qPCR_HLBech_ct36_arbre <- factor(test_set$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
    matrice <- confusionMatrix(test_set$qPCR_HLBech_ct36_arbre, y_pred)
    
    # y_pred = predict(classifier, newdata = test_fold[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))])
    # test_fold$qPCR_HLBech_ct36_arbre <- factor(test_fold$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
    # matrice <- confusionMatrix(test_fold$qPCR_HLBech_ct36_arbre, y_pred)
    #return(matrice)
    accuracy <- accuracy + matrice$overall[1] #/10
    sensi <- sensi + matrice$byClass[1] #/10
    pospred <- pospred + matrice$byClass[3] #/10
    #}
    data_global[c,1] <- accuracy
    data_global[c,2] <- sensi
    data_global[c,3] <- pospred
  }
  
  sd_rep  <-   rbind(sd_rep,c(sd(data_global$accuracy), sd(data_global$pospred), sd(data_global$sensi)))   #rbind(sd_rep,sqrt((colMeans(data_global^2))-((colMeans(data_global))^2)))
  #var_rep <-  rbind(var_rep, (sd_rep)^2)  #rbind(var_rep,(colMeans(data_global^2))-((colMeans(data_global))^2))
  mean_rep <- rbind(mean_rep,colMeans(data_global))
  print(a) # ce qui est renvoye a l'ecran et correspond au nombre de rep_feuille
}



#names(var_rep)[1:3] <- c("a_var","s_var","p_var")

names(sd_rep)[1:3] <- c("a_sd","s_sd","p_sd")

names(mean_rep)[1:3] <- c("Accuracy","Sensitivity","Precision")


# data_global_rep <- cbind(mean_rep, var_rep,sd_rep)
# data_global_rep$nb <- seq(nrow(data_global_rep)) 
# 
# mean_rep_long <- pivot_longer( data =data_global_rep, cols = 1:(ncol(mean_rep)) ) #, values_to = "reflectance", names_to = "lambda"  )

mean_rep$nb_rep <- seq(nrow(mean_rep)) 

mean_rep_long <- melt(mean_rep,id.vars = "nb_rep", measure.var=colnames(mean_rep[1:3]))

sd_rep$ns <- seq(nrow(sd_rep))

sd_rep_long <- melt(sd_rep,id.vars = "ns", measure.var=colnames(sd_rep[1:3]))

# var_rep$nv <- seq(nrow(var_rep))
# 
# var_rep_long <- melt(var_rep,id.vars = "nv", measure.var=colnames(var_rep[1:3]))


data_global_rep <- cbind(mean_rep_long ,sd_rep_long )

data_global_rep<- data_global_rep[,-(which(colnames(data_global_rep) == "ns" ))]

names(data_global_rep)[2:5] <- c("Mean_values","Mean","Sd_values", "Sd")


#esquisse::esquisser() 

IC <- (1.96*data_global_rep$Sd)/sqrt(l_c)

# faire graph !!

names(data_global_rep)

graph_rep <- ggplot(data=data_global_rep, aes(nb_rep, Mean, colour= Mean_values))+
  scale_x_continuous(breaks=seq(1:6))+
  geom_errorbar(aes(ymin = Mean-IC , ymax = Mean+IC),width=0.05, lwd = 1.1 )+
  geom_point(size = 3)+
  geom_line(lwd = 1.1) +
  #scale_color_brewer(palette = "Dark2") +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction du statut HLB  des abres par SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM pour chaque répétition de feuille ", color = "Paramètres de robustesse en SVM ") +
  #theme_bw() +
  dark_theme_gray() +
  theme(legend.position = "right")+ # mettre "bottom" pour avoir titre en bas
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


graph_rep

#Faire vrai écart type IC = 1,96*sd/racine de n (boucle 10 ou 100) , voir fonction t.test


#save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_100byRep_Ed_ss_test_fold") 


#write.table(x =data_global_rep  , file = "Donnees/SVM_ct36_6rep_100byRep_Ed_ss_test_fold.csv")

#ggsave(filename = "Graphiques/Qualité_de_la_prédiction_SVM_ct36_6rep_100byRep_Ed_ss_test_fold.png", plot = last_plot())

#getwd()  # pour connaitre le cheminemant des fichiers


# Regarder la foinction progress() de library(mgcv)

# Faire un codage et analyse par lot


# Extraction des paramètre SVM du statut HLB à CT36 avec variation du nombre de feuille ####

#Essaie pour filtrer en fonction d'un nombre de feuille total choisi en piochant dans les listes le nombre de feuille voulu sachant que chaque tirroire = 1 feuille

tmb=split(data_SPIR_Ed,(data_SPIR_Ed$code_ech_feuille)) 



l_c = 100  # nombre de creation training_set et test_set 
l_a = 6   # nombdre de rep par feuille (sur un total de 6 repetition)
l_Z = 420   # nombre  de feuilles (sur un total de 420 feuilles)



data_global <- data.frame(accuracy= rep(NA,l_c), sensi=rep(NA,l_c), pospred=rep(NA,l_c))
var_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
sd_rep <-  data.frame(accuracy=c(), sensi=c(), pospred=c())
mean_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c())
#data_global_rep <- data.frame(accuracy=c(), sensi=c(), pospred=c(), nb=c())



progress = function(a, length_a){
  print(a/length_a)
}

r = 0
stock = 0



for (a in 1:l_a) {
  
  tmp <- sample(x = tmb,size = l_Z ,replace = FALSE,prob = NULL)
  
  maliste =lapply(1:length(tmp),function(b){ 
    tmp[[b]][sample(1:length(tmp[[b]]$code_ech_feuille),a),]   # on fait ça pour tte les feuilles mais tjrs en choisissant 2 rep tire aleatoirement
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
    
    
    
    # Désormais on prédit grâce au modèle classifier le statut qPCR_HLBech_ct36_arbre en utilisant le test_FOLD et on calcule la précision du modèle
    
    y_pred = predict(classifier, newdata = test_set[,c ((which(colnames(test_set) == "X350" )): (which(colnames(test_set) == "X2500" )))])
    test_set$qPCR_HLBech_ct36_arbre <- factor(test_set$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
    matrice <- confusionMatrix(test_set$qPCR_HLBech_ct36_arbre, y_pred)
    
    
    accuracy <- accuracy + matrice$overall[1] #/10
    sensi <- sensi + matrice$byClass[1] #/10
    pospred <- pospred + matrice$byClass[3] #/10
    #}
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


#esquisse::esquisser() 

# Graphique SVM ct36

# data_global_rep <- read.table(file = "Donnees/SVM_ct36_1lot_6rep_10byRep_Ed_ss_test_fold_420f.csv"
#                               
#                               , header = T
#                               , sep = ""
#                               , stringsAsFactors = T
#                               , row.names = 1
#                               , na.strings = c("","NA")
#                               , dec = "." )
# l_c = 100
# 
# IC <- (1.96*data_global_rep$Sd)/sqrt(l_c)
# 


names(data_global_rep)

graph_ct36_420f <- ggplot(data=data_global_rep, aes(nb_rep, Mean, colour= Mean_values))+
  scale_x_continuous(breaks=seq(1:6))+
  geom_errorbar(aes(ymin = Mean-IC , ymax = Mean+IC),width=0.05, lwd = 1.1 )+
  geom_point(size = 3)+
  geom_line(lwd = 1.1) +
  #scale_color_brewer(palette = "Dark2") +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction du statut HLB  des abres par SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM pour chaque répétition de feuille (420f)", color = "Paramètres de robustesse en SVM ") +
  #theme_bw() +
  dark_theme_gray() +
  theme(legend.position = "right")+ # mettre "bottom" pour avoir titre en bas
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


graph_ct36_420f



#save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_1lot_6rep_10byRep_Ed_ss_test_fold_420f") 


#write.table(x =data_global_rep  , file = "Donnees/SVM_ct36_1lot_6rep_10byRep_Ed_ss_test_fold_420f.csv")

#ggsave(filename = "Graphiques/SVM/Qualité_de_la_prédiction_SVM_ct36_1lot_6rep_10byRep_Ed_ss_test_fold_420f.png", plot = last_plot())


# Extraction des paramètre SVM du statut HLB à CT32 avec variation du nombre de feuille ####

tmb=split(data_SPIR_Ed,(data_SPIR_Ed$code_ech_feuille)) 



l_c = 100  # nombre de creation training_set et test_set 
l_a = 6   # nombdre de rep par feuille (sur un total de 6 repetition)
l_Z = 420   # nombre  de feuilles (sur un total de 420 feuilles)



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
    
    
    print_a = paste("Le programme en est à ", substr(stock + a*(100*(c-1)/l_c/(r-a+1)), start = 1, stop = 4) , "% de sa progression", sep = "")
    if (c == l_c){stock = 100*(c-1)/l_c/(r-a+1)+stock}
    print(print_a)
    
    #set.seed(123)
    
    split = sample.split(test_ed$qPCR_HLBech_ct32_arbre, SplitRatio = 0.75) 
    
    training_set <- subset(test_ed, split == TRUE) 
    test_set <- subset(test_ed, split == FALSE) 
    
    
    accuracy <- 0
    pospred <- 0
    sensi <- 0
    
    
    classifier = svm(y = training_set$qPCR_HLBech_ct32_arbre, x = training_set[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))],
                     data = training_set,
                     type = 'C-classification',
                     kernel = 'linear')
    
    
    
    y_pred = predict(classifier, newdata = test_set[,c ((which(colnames(test_set) == "X350" )): (which(colnames(test_set) == "X2500" )))])
    test_set$qPCR_HLBech_ct32_arbre <- factor(test_set$qPCR_HLBech_ct32_arbre)           
    matrice <- confusionMatrix(test_set$qPCR_HLBech_ct32_arbre, y_pred)
    
    
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



# Graphique des paramètre de SVM ct32
# 
# data_global_rep <- read.table(file = "Donnees/SVM_ct32_6rep_100byRep_Ed_ss_test_fold_420f.csv"
#                               
#                               , header = T
#                               , sep = ""
#                               , stringsAsFactors = T
#                               , row.names = 1
#                               , na.strings = c("","NA")
#                               , dec = "." )
# l_c = 100
# 
# IC <- (1.96*data_global_rep$Sd)/sqrt(l_c)

graph_ct32_420f <- ggplot(data=data_global_rep, aes(nb_rep, Mean, colour= Mean_values))+
  scale_x_continuous(breaks=seq(1:6))+
  geom_errorbar(aes(ymin = Mean-IC , ymax = Mean+IC),width=0.05, lwd = 1.1 )+
  geom_point(size = 3)+
  geom_line(lwd = 1.1) +
  labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction du statut HLB  des abres par SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<32), obtenu après avoir fait la moyenne de 100 SVM pour chaque répétition de feuille (420f)", color = "Paramètres de robustesse en SVM ") +
  
  dark_theme_gray() +
  theme(legend.position = "right")+
  scale_colour_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "grey20"))


graph_ct32_420f

#save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct32_6rep_100byRep_Ed_ss_test_fold_420f") 


#write.table(x =data_global_rep  , file = "Donnees/SVM_ct32_6rep_100byRep_Ed_ss_test_fold_420f.csv")

#ggsave(filename = "Graphiques/SVM/Qualité_de_la_prédiction_SVM_ct32_6rep_100byRep_Ed_ss_test_fold_420f.png", plot = last_plot())


#groupement des graphiques

#plot_grid_SVM <- plot_grid(graph_ct36_420f , graph_ct32_420f)

#plot_grid_SVM


# SVM : Graphique ####


l_a = 6   # nombdre de rep par feuille (sur un total de 6 repetition)
l_Z = 420   # nombre  de feuilles (sur un total de 420 feuilles)


tmb=split(data_SPIR_Ed,(data_SPIR_Ed$code_ech_feuille)) 

for (a in 1:l_a) {
  
  tmp <- sample(x = tmb,size = l_Z ,replace = FALSE,prob = NULL)
  
  maliste =lapply(1:length(tmp),function(b){ 
    tmp[[b]][sample(1:length(tmp[[b]]$code_ech_feuille),a),]   # on fait ça pour tte les feuilles mais tjrs en choisissant 2 rep tire aleatoirement
  })
  
  test_ed=do.call(rbind,maliste)
  
  split = sample.split(test_ed$qPCR_HLBech_ct36_arbre, SplitRatio = 0.75) 
  
  training_set <- subset(test_ed, split == TRUE) 
  test_set <- subset(test_ed, split == FALSE) 
  
}


classifier = svm(y = training_set$qPCR_HLBech_ct36_arbre, x = training_set[,c ((which(colnames(training_set) == "X350" )): (which(colnames(training_set) == "X2500" )))],
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'linear')

y_pred = predict(classifier, newdata = test_set[,c ((which(colnames(test_set) == "X350" )): (which(colnames(test_set) == "X2500" )))])
test_set$qPCR_HLBech_ct36_arbre <- factor(test_set$qPCR_HLBech_ct36_arbre)            # pour mettre test_fold$qPCR_HLBech_ct36_arbre en facteur pour la fonction cufusionMatrix
matrice <- confusionMatrix(test_set$qPCR_HLBech_ct36_arbre, y_pred)

plot(classifier,training_set  )

print(classifier)

summary(classifier)

ggplot(classifier)

# Test discrimination du jeu de donnee par Random Forest ####

ct36_arbre <- ((data_SPIR_Ed [, -(which(colnames(data_SPIR_Ed) == "code_ech_feuille" | colnames(data_SPIR_Ed) == "code_ech_lot" | colnames(data_SPIR_Ed) == "code_ech_arbre" | colnames(data_SPIR_Ed) == "rep_arbre" | colnames(data_SPIR_Ed) =="rep_lot" | colnames(data_SPIR_Ed) =="rep_feuille" | colnames(data_SPIR_Ed) =="HLB" | colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_arbre"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_lot"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct36_lot"))])) 

summary(ct36_arbre$qPCR_HLBech_ct36_arbre)

gtree <- ctree( qPCR_HLBech_ct36_arbre ~ . , data=ct36_arbre)

# plot 

plot(gtree ,)

gtree <- ctree( qPCR_HLBech_ct36_arbre ~ . , data=ct36_arbre)

plot(gtree, inner_panel = node_barplot,
     edge_panel = function(ctreeobj, ...) { function(...) invisible() },
     tnex = 1)

nodes(gtree, 1)

table(Predict(gtree), ct36_arbre$qPCR_HLBech_ct36_arbre)

# Prediction

tr_tree <- lmtree( qPCR_HLBech_ct36_arbre ~ . , data=ct36_arbre)


# Test discrimination du jeu de donnee par Partial Least Square ####

ct36_arbre <- ((data_SPIR_Ed [, -(which(colnames(data_SPIR_Ed) == "code_ech_feuille" | colnames(data_SPIR_Ed) == "code_ech_lot" | colnames(data_SPIR_Ed) == "code_ech_arbre" | colnames(data_SPIR_Ed) == "rep_arbre" | colnames(data_SPIR_Ed) =="rep_lot" | colnames(data_SPIR_Ed) =="rep_feuille" | colnames(data_SPIR_Ed) =="HLB" | colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_arbre"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_lot"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct36_lot"))])) 

ct36_arbre$qPCR_HLBech_ct36_arbre <- as.numeric(ct36_arbre$qPCR_HLBech_ct36_arbre)


split = sample.split(ct36_arbre$qPCR_HLBech_ct36_arbre, SplitRatio = 0.75) 

train.data <- subset(ct36_arbre, split == TRUE)
test.data <- subset(ct36_arbre, split == FALSE)

preproc.param <- train.data %>%  preProcess(method = c("center", "scale"))

train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)



pls1 = plsreg1(ct36_arbre[, 2:2152], ct36_arbre[, 1, drop = FALSE], comps = 3)

plot(pls1)

plot(ct36_arbre$qPCR_HLBech_ct36_arbre, pls1$y.pred, xlab="Original", ylab = "Predicted")

# Test discrimination du jeu de donnee par Principal Component Regression ####

ct36_arbre <- ((data_SPIR_Ed [, -(which(colnames(data_SPIR_Ed) == "code_ech_feuille" | colnames(data_SPIR_Ed) == "code_ech_lot" | colnames(data_SPIR_Ed) == "code_ech_arbre" | colnames(data_SPIR_Ed) == "rep_arbre" | colnames(data_SPIR_Ed) =="rep_lot" | colnames(data_SPIR_Ed) =="rep_feuille" | colnames(data_SPIR_Ed) =="HLB" | colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_arbre"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_lot"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct36_lot"))])) 

require(pls)

pcr_model <- pcr(qPCR_HLBech_ct36_arbre ~ . , data=ct36_arbre, scale = TRUE, validation = "CV")

plot(pcr_model)

validationplot(pcr_model)

predplot(pcr_model)

# Test discrimination du jeu de donnee par Discriminant Analysis Essentials ####

ct36_arbre <- ((data_SPIR_Ed [, -(which(colnames(data_SPIR_Ed) == "code_ech_feuille" | colnames(data_SPIR_Ed) == "code_ech_lot" | colnames(data_SPIR_Ed) == "code_ech_arbre" | colnames(data_SPIR_Ed) == "rep_arbre" | colnames(data_SPIR_Ed) =="rep_lot" | colnames(data_SPIR_Ed) =="rep_feuille" | colnames(data_SPIR_Ed) =="HLB" | colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_arbre"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct32_lot"| colnames(data_SPIR_Ed) =="qPCR_HLBech_ct36_lot"))])) 

split = sample.split(ct36_arbre$qPCR_HLBech_ct36_arbre, SplitRatio = 0.75) 

train.data <- subset(ct36_arbre, split == TRUE)
test.data <- subset(ct36_arbre, split == FALSE)

preproc.param <- train.data %>%  preProcess(method = c("center", "scale"))

# creuser si preproc à besoin de centrer reduire la colonne HLB

train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

## a) Linear discriminant analysis - LDA ####

# Fit the model

model_lda <- lda(qPCR_HLBech_ct36_arbre~., 
                 data = ct36_arbre) 
                 #prior = c(1,1)/2, 
                 #subset = train.transformed)

#print(model_lda)

names(model_lda$scaling)

#lda.data<-cbind(train.transformed, predict(model_lda)$x)

ggplot(model_lda, aes(LD1,LD2))

names(lda.data)

coef(model_lda) # corespond à LD1 
model_lda[["scaling"]][2,] 

#model

plot_model_lda <- plot(model_lda) # 2 groups pour 0 et 1
plot_model_lda

# Make predictions

predictions_lda <- model_lda %>% predict(test.transformed)


str(model_lda)

# Predicted classes
head(predictions_lda$class, 10)
# Predicted probabilities of class memebership.
head(predictions_lda$posterior, 10) 
# Linear discriminants
head(predictions_lda$x, 10) 
#names(predictions)

# Model accuracy

m_lda <- mean(predictions_lda$class==test.transformed$qPCR_HLBech_ct36_arbre) ; m_lda

lda.data_lda <- cbind(train.transformed, predict(model_lda)$x) # predict(model)$x = LD1 mais où est LD2 ?

RowMeans <- rowMeans(train.transformed[,-1])

lda.data_lda <- cbind(lda.data_lda, RowMeans)

View(summary(lda.data_lda ))

#esquisse::esquisser() 

ggplot(lda.data_lda, aes(LD1, RowMeans)) + 
  geom_point(aes(color = qPCR_HLBech_ct36_arbre)) +
  labs(x = "Linear discriminant analysis  ", y = "Valeurs moyennes", title = "Test discrimination du jeu de donnee par Linear discriminant analysis ", subtitle = "", color = "Qpcr à Ct36 ") +
  dark_theme_gray() +
  #scale_colour_viridis_d(option = "plasma") +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major.y = element_line(colour = "grey20"))



## b) Quadratic discriminant analysis - QDA #### 



model_qda <- qda(qPCR_HLBech_ct36_arbre~., data = train.transformed) 


plot(model_qda) 


predictions_qda <- model_qda %>% predict(test.transformed)

#names(predictions_qda)


m_qda <- mean(predictions_qda$class==test.transformed$qPCR_HLBech_ct36_arbre) ; m_qda

lda.data_qda <- cbind(train.transformed, predict(model_qda)$x)



## c) Mixture discriminant analysis - MDA ####



model_mda <- mda(qPCR_HLBech_ct36_arbre~., data = train.transformed) 


plot(model_mda) 


predictions_mda <- model_mda %>% predict(test.transformed)

#names(predictions_mda)

m_mda <- mean(predictions_mda$class==test.transformed$qPCR_HLBech_ct36_arbre) ; m_mda

lda.data_mda <- cbind(train.transformed, predict(model_mda)$x)



## d) Flexible discriminant analysis - FDA ####



model_fda <- fda(qPCR_HLBech_ct36_arbre~., data = train.transformed) 


plot(model_fda) 


predictions_fda <- model_fda %>% predict(test.transformed)

#names(predictions_fda)

m_fda <- mean(predictions_fda$class==test.transformed$qPCR_HLBech_ct36_arbre) ; m_fda

lda.data_fda <- cbind(train.transformed, predict(model_fda)$x)



## e) Regularized discriminant analysis - RDA ####



model_rda <- rda(qPCR_HLBech_ct36_arbre~., data = train.transformed) 


plot(model_rda) # 2 groups pour 0 et 1


predictions_rda <- model_rda %>% predict(test.transformed)

#names(predictions_rda)

m_rda <- mean(predictions_rda$class==test.transformed$qPCR_HLBech_ct36_arbre) ; m_rda

lda.data_rda <- cbind(train.transformed, predict(model_rda)$x)


