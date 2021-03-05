
rm(list=ls())

# LIBRARY ####

library(caret) # package Machine learning # confusion matrix
library(caTools) # sample.split
library(data.table)
library(e1071) # SVM
library(ggplot2)
library(ggparty)
library(ggdark)
library(ggpubr)
library(gt) 
library(plyr)
library(party)  # Random Forest
library(pls)
library(plsdepot) # Partial Least Square
library(Rmisc) #intervalle de confiance
library(tidyr) # pivot longer
library(tidyverse) 
library(snowfall)

library(randomForest)
library(cvAUC)

# Importation des donnees data_SPIR_Ed ####

data_SPIR_Ed <- read.table(file = "Donnees/Data_SPIR_Edouard_Hoarau/SPIR_Hoarau.csv"
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )

# Attention des qu'il y a du texte(ou des virgules) tt le jeu de donnee est compris comme etant non numerique , il faut donc le mettre en numerique !


# data_SPIR_Ed
code_labo <- rownames(data_SPIR_Ed)

data_SPIR_Ed$code_nbr_rep <- factor(substr(code_labo, 8, 8))
data_SPIR_Ed$code_ech_arbre <- factor(substr(code_labo, 1, 2))
data_SPIR_Ed$code_ech_feuille <- factor(substr(code_labo, 1, 3))
data_SPIR_Ed$code_rep_feuille <- factor(paste(data_SPIR_Ed$code_ech_feuille, data_SPIR_Ed$code_nbr_rep,sep="")) 

rm (code_labo)

data_SPIR_Ed <- data_SPIR_Ed[ !is.na(data_SPIR_Ed$X350),]

#length(levels(data_SPIR_Ed$code_ech_arbre)) # correspondant à 2x HLB 3x lot et 7x arbres soit 42

# Importation et preparation des resultat de la Qpcr ####

data_Qpcr_Ed <- read.table(file = "Donnees/Resultats_qPCR_Ed/Detect_HLB_Hoareau_Ed_19022021_data.csv"
                           
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )

# on stocke les noms d'échantillons positifs selon nos 2 seuils de Ct (cycle de qPCR ), a savoir moins de 32 cyclces et moins de 36 cycles qPCR
seuils <- c(32, 36)

trueP <- lapply(seuils, function(x) unique(data_Qpcr_Ed$Sample.Name[ which(data_Qpcr_Ed$C..Mean < x & data_Qpcr_Ed$C..SD < 1) ]))
names(trueP) <- paste("seuil", seuils, sep = ".")          
              

data_SPIR_Ed[paste("qPCR", seuils, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(data_SPIR_Ed$code_ech_arbre %in% x ) ) # on cherche quels sont les code_ech_arbre qui se trouvent dans le vecteur x, x reprenant automatiquement les noms des arbres positifs selon le seuil choisi

select.lambda <- grep("^X", names(data_SPIR_Ed))
data_SPIR_Ed <- data_SPIR_Ed[,c(names(data_SPIR_Ed)[-select.lambda], names(data_SPIR_Ed)[select.lambda] )]

rm (data_Qpcr_Ed,trueP,select.lambda,seuils)

data_SPIR_Ed[c("qPCR_32", "qPCR_36")] <- lapply(data_SPIR_Ed[c("qPCR_32", "qPCR_36")], factor)

# Exploration du jeu de donnees     ####

ftable (code_nbr_rep ~ code_ech_arbre , data = data_SPIR_Ed) # repartition global pour voir les erreur de prise de données

ftable ( code_ech_arbre ~ qPCR_32 , data = data_SPIR_Ed) 

ftable ( code_ech_arbre ~ qPCR_36 , data = data_SPIR_Ed) 

# Format_long ####

select.lambda <- grep("^X", names(data_SPIR_Ed))
data_long <- pivot_longer( data = data_SPIR_Ed, cols = select.lambda, values_to = "reflectance", names_to = "lambda"  ) 


data_long$lambda <- as.numeric(gsub("X", "", data_long$lambda))
summary(data_long)

data_long <- data_long[ !is.na(data_long$reflectance),]

# graph moyennés à CT 32 ####

ggplot(data_long) +
  aes(x = lambda, y = reflectance, group = qPCR_32, color = qPCR_32 ) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_brewer(palette = "Dark2", labels = c("Négatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif à Ct<32)", color = "Résultat du test HLB à Ct<32") +
  #dark_theme_gray() 
  theme_gray()

# graph zoomé sur les spectres 400 à 680 nm

ggplot(data_long[data_long$lambda >= 400 & data_long$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_32 )+
  # geom_line() +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  theme_gray()

# graph moyennés à CT 36 ####

ggplot(data_long) +
  aes(x = lambda, y = reflectance, group = qPCR_36, color = factor(qPCR_36) ) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif à Ct<36)", color = "Resultat du test HLB à Ct<36") +
  #dark_theme_gray() 
  theme_gray()

# graph zoomé sur les spectres 400 à 680 nm

ggplot(data_long[data_long$lambda >= 400 & data_long$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = factor(qPCR_36) )+
  # geom_line() +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  theme_gray()

# multiple graph feuilles moyennés à CT 32 ####

p1 <- ggplot(data_long[data_long$qPCR_32 == 1,]) +
  aes(x = lambda, y = reflectance, color = code_ech_feuille, group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  F) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyenn des feuilles par arbres", subtitle = "arbre positif à Ct<32)", color = "Resultat du test HLB à Ct<36") +
  #dark_theme_gray() 
  theme_gray()



p0 <- ggplot(data_long[data_long$qPCR_32 == 0,]) +
  aes(x = lambda, y = reflectance, color = code_ech_feuille, group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  F) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyenn des feuilles par arbres", subtitle = "arbre négatif à Ct<32)", color = "Resultat du test HLB à Ct<36") +
  #dark_theme_gray() 
  theme_gray()

ggarrange(p0, p1, ncol = 2)

#ggsave("Graphiques/exploration Hoareau spectres moyens feuille selon sante.pdf", units = "cm", width = 20, height = 15, scale = 2)

# test de Student ####
data_long.lambda <- split(data_long, data_long$lambda)

test.de.student <- lapply(data_long.lambda, function(x) t.test(x$reflectance[x$QPCR.32 == 0], x$reflectance[x$QPCR.32 == 1], var.equal = F ))

toto <- sapply(test.de.student, function(x) x$p.value)

longueurs <- 350:2500
plot(longueurs, toto, type= "l")
abline(v = longueurs[which.max(toto)], col = "red", lwd = 2)
text(longueurs[which.max(toto)], 4, longueurs[which.max(toto)], adj = -0.25, font = 3, cex = 2)

toto <- names(test.de.student)[ sapply(test.de.student, function(x) x$p.value < 0.001 )]


table(data_long$QPCR.32)/(2500 - 350 + 1) # nb d'echantillons (spectres) supposes malades et sains

# Faire modele mixt ( reflectance ~ QPCR. )
# car donnÃ©es pas idÃ©pendante :   Verger/arbre/feuille/Repetition

#->Support Vector Machine en kernel linear ####


source(file = "Scripts/Fonction_SVM.R")  




## a) SVM a Ct32 & Graph ####

nb.simu <- 5  # Minimu 1000 simu
rep.max <- 3


Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T) # drop = T pour enlever les tiroirs vides !!

#Blablabla

### a.1) Calcul parallele ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs
sfLibrary(caTools)             # la library des packages utilisés
sfLibrary(e1071)
sfLibrary(caret)
sfExport("fct_svm","Tirage","rep.max","nb.simu") # les éléments extérieur à la fonction
T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

 
res.svm.32 <- sfClusterApplySR(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage , restore = F, perUpdate = 6 ) # restore = T seulement si ça plante !

#res.svm.32 <- sfClusterApplyLB(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage ) # Le plus rapide

# Commande sans calcul parralele  # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

#res.svm.32 <- sfSapply(rep(1:rep.max, each = nb.simu), fct_svm , seuil.ct = 32 , list.feuilles= Tirage)

T2 <- Sys.time()

sfStop()        # stop l'utilisation du sfInit aux autres lignes de codes


difftime(T2,T1) # information sur le temps qu'à mis l'operation 

# Time difference of 4.272735 mins avec nb.simu = 5 ,rep.max = 3 & sfSapply

# Time difference of 4.264103 mins avec nb.simu = 5 ,rep.max = 3 & sfClusterApplyLB

# Time difference of 4.569168 mins avec nb.simu = 5 ,rep.max = 3 & sfClusterApplySR


### a.2) Enregistrement des criteres de precision ####

#load("Sauvegardes_objet_R.data/SVM_ct32_6rep_10simu_para")

data_global.32 <- pivot_longer(as.data.frame(t(res.svm.32)), cols = 1:3, names_to = "critere", values_to = "valeurs")
data_global.32$nb.rep <- rep(1:rep.max, each = rep.max*nb.simu)


tab.32 <-  aggregate(valeurs ~ nb.rep + critere, data_global.32, mean) 
names(tab.32)[3] <- "moyenne"
tab.32$et <- aggregate(valeurs ~ nb.rep + critere, data_global.32, sd)$valeurs
tab.32$mediane <- aggregate(valeurs ~ nb.rep + critere, data_global.32, median )$valeurs
tab.32$nb <-  aggregate(valeurs ~ nb.rep + critere, data_global.32, length )$valeurs


# save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct32_6rep_10simu_para") 

# write.table(x =data_global.32 , file = "Donnees/SVM_ct32_6rep_100simu.csv" , sep = ';')


## a.3) graph calculé avec stat_summary ####

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

## a.4) graph calculé avec IC ####

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

    
## b) SVM a Ct36 & Graph ####

sfInit(parallel = T, cpus = 4) # optimisation des processeurs sur les 4 coeurs

T1 <- Sys.time() # information sur le temps que met l'operation a se realiser

res.svm.36 <- lapply(1:rep.max, function(nomb.rep)
  sapply(1:nb.simu, function(x) fct_svm(nb.rep = nomb.rep, seuil.ct = 36, list.feuilles = Tirage)))  # fct_svm on remplsis les 3 arguments qui sont : nb.rep, seuil.ct, list.feuilles

T2 <- Sys.time()

sfStop()
difftime(T2,T1)

    data_global.36 <- lapply(res.svm.36, function(x) data.frame(critere = c("Accuracy","Sensitivity","Precision")
                                                                , moyenne = apply(x, 1, mean)
                                                                , et = apply(x, 1, sd)
                                                                , mediane = apply(x, 1, median)) )
    data_global.36 <- do.call(rbind, data_global.36)
    data_global.36$nb.rep <- rep(1:rep.max, each = 3)
    
    IC <- 1.96*data_global.36$et/sqrt(nb.simu)
    
    ggplot(data = data_global.36) +
      aes(x = nb.rep, y = moyenne, colour= critere, group = critere)+
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = moyenne - IC , ymax = moyenne + IC),width=0.05, lwd = 1.1 )+
      #stat_summary(fun.data  = function(x) mean_se(x, mult = qt(0.975, length(x) -1)), color = "pink") +
      scale_x_continuous(breaks=seq(1:6)) +  
      geom_line(lwd = 1.1) +
      labs(x = "Nombre de mesures SPIR effectuées sur chaque feuille ", y = "Valeurs moyennes", title = "Prédiction des paramètres de SVM en fonction du nombre de mesure SPIR par feuille ", subtitle = "Détection du statut HLB à (Ct<36), obtenu après avoir fait la moyenne de 100 SVM", color = "Paramètres de robustesse en SVM ") +
      dark_theme_gray() +
      theme(legend.position = "right")+
      scale_colour_viridis_d() +
      theme(panel.grid.major.y = element_line(colour = "grey20"))
  
    save(list = ls(), file = "Sauvegardes_objet_R.data/SVM_ct36_6rep_100simu") 
    
    write.table(x =data_global.36 , file = "Donnees/SVM_ct36_6rep_100simu.csv" , sep = ';')

    ggsave(filename = "Graphiques/SVM_ct36_6rep_100simu.png", plot = last_plot() ,width = 900 ,height = 680 )

rm (fct_svm,Tirage,nb.simu,rep.max,IC)


#->Random Forest ####

#nb.simu <- 2

nb.rep <- 3
seuil.ct <- 32

Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)


maliste <- lapply(Tirage,function(feuille){ 
  list_rf <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
  # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
  code <- grep("^code_ech", names(list_rf))
  qPCR <- grep("^qPCR_", names(list_rf))
  sortie <- cbind.data.frame(unique(list_rf[c(code,qPCR)])
                             , matrix(apply(list_rf[-(which(colnames(list_rf) == "code_nbr_rep"):which(colnames(list_rf) == "qPCR_36"))], 2, mean)
                                      , nr = 1, dimnames = list(NULL, names(list_rf)[-(which(colnames(list_rf) == "code_nbr_rep"):which(colnames(list_rf) == "qPCR_36"))])))
  sortie
})

test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 

decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
train_rf <- test_ed[decoup,]
test_rf <- test_ed[!decoup,] 
#}



## a) Essais avec le jeu de donne global ####
    
    RF_ct32 <- data_SPIR_Ed
    RF_ct32$qPCR_32 <- as.numeric(RF_ct32$qPCR_32)
    masque_numeric <- sapply(RF_ct32, is.numeric)
    RF_ct32 <- RF_ct32[,masque_numeric]

    gtree_glob <- ctree( qPCR_32 ~ . , data=RF_ct32)
    
    # plot 
    
    plot(gtree_glob ,)
    
    plot(gtree_glob, inner_panel = node_barplot,
         edge_panel = function(ctreeobj, ...) { function(...) invisible() },
         tnex = 1)
    
    nodes(gtree_glob, 1)
    
    
    
## b) Essais avec le training et test ####
    
    train_rf$qPCR_32 <- as.numeric(train_rf$qPCR_32)
    masque_numeric <- sapply(train_rf, is.numeric)
    train_rf_32 <- train_rf[,masque_numeric]
    train_rf_32$qPCR_32 <- train_rf$qPCR_32
    
    gtree_32 <- ctree( qPCR_32 ~ . , data=train_rf_32)
    
    plot(gtree_32 ,)
    
    plot(gtree_32, inner_panel = node_barplot,
         edge_panel = function(ctreeobj, ...) { function(...) invisible() },
         tnex = 1)
    
    # Prediction
    
    test_rf$qPCR_32 <- as.numeric(test_rf$qPCR_32)
    masque_numeric <- sapply(test_rf, is.numeric)
    test_rf_32 <- test_rf[,masque_numeric]
    
    tr_tree <- lmtree( qPCR_32 ~ . , data=test_rf_32)
    
    plot(tr_tree ,)
    
    table(Predict(gtree_32), test_rf_32$qPCR_32)
    
## c) Package randomForest ####

    
    
    
    
#fct_rf <- function(nb.rep, seuil.ct, list.feuilles) {
    
   
    # Identity the response column
    
    Resp <- which(colnames(train_rf) == "qPCR_32")
    
    # Identify the predictor columns
    
    Pred <- grep("^X", names(train_rf))
    

    #set.seed(1)  # For reproducibility
    
    system.time(
      model_rf <- randomForest(
        x = train_rf[,Pred], 
        y = train_rf[,Resp],
        xtest = test_rf[,Pred],
        ntree = 100
      )
    )
    
    # preds <- model$test_rf$qPCR_32[, 2]
    # labels <- test_rf[,Resp]
    # cvAUC::AUC(predictions = preds, labels = labels)
    
    plot(model_rf)
    
    rf_classifier = randomForest(qPCR_32 ~ ., data=train_rf, ntree=100, mtry=2, importance=TRUE)
    
    plot(rf_classifier)
    
## d) Package carret ####
    
    model_carret <- caret::train(
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
    
    
#->Partial Least Square ####
    
    nb.rep <- 3
    seuil.ct <- 32
    
    Tirage <- split(data_SPIR_Ed, data_SPIR_Ed$code_ech_feuille, drop = T)
    
    
    maliste_pls <- lapply(Tirage,function(feuille){ 
      list_pls <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
      # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
      code <- grep("^code_ech", names(list_pls))
      qPCR <- grep("^qPCR_", names(list_pls))
      sortie <- cbind.data.frame(unique(list_pls[c(code,qPCR)])
                                 , matrix(apply(list_pls[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))], 2, mean)
                                          , nr = 1, dimnames = list(NULL, names(list_pls)[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))])))
      sortie
    })
    
    test_ed <- do.call(rbind, maliste_pls) # on colle toute les listes créées précédement 
    
    decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
    train_pls <- test_ed[decoup,]
    test_pls <- test_ed[!decoup,] 


## a) Essais avec le jeu de donne global ####    
   
    Pls_ct32 <- data_SPIR_Ed
    Pls_ct32$qPCR_32 <- as.numeric(Pls_ct32$qPCR_32)
    masque_numeric <- sapply(Pls_ct32, is.numeric)
    Pls_ct32 <- Pls_ct32[,masque_numeric]
    
    pls.model = plsr(qPCR_32 ~ ., data = Pls_ct32, scale = TRUE, validation = "CV")
    
    # Prediction
    
    pls.pred = predict(pls.model, test_pls, ncomp=100)
    
    plot( pls.pred,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
    abline( col="red")
    
    
    
    
    # Graph
    
    #load("Sauvegardes_objet_R.data/Essais_pls_model_data_global")
    
    pls1 <- plot(pls.model)
    
    pls.RMSEP = RMSEP(pls.model, estimate="CV")
    pls2 <- plot(pls.RMSEP, main="RMSEP PLS Solubility", xlab="components")
    
    pls3 <- validationplot(pls.model)
    
    
    # min_comp = which.min(pls.model$val)
    # points(min_comp, min(pls.model$val), pch=1, col="red", cex=1.5)
    
    
    coefficients = coef(pls.model)
    sum.coef = sum(sapply(coefficients, abs))
    coefficients = coefficients * 100 / sum.coef
    coefficients = sort(coefficients[, 1 , 1])
    barplot(tail(coefficients, 10))
    barplot(head(coefficients, 10))
    
    coef_pls <- data.frame(coefficients)
    coef_pls$Lambda <- as.numeric(substr(rownames(coef_pls), 2, 5))
    
    
    pls4 <- ggplot(data = coef_pls ) +
      aes(x = Lambda, y = coefficients) +
      #geom_point(colour = "#0c4c8a") +
      geom_text(label=rownames(coef_pls))+
      labs(title = "Coefficient du model PLS sur le jeu de donné global", subtitle = "A Ct32") +
      theme_minimal()
    
    pls4
    
     
## b) Essais avec le training et test ####    
    
    train_pls$qPCR_32 <- as.numeric(train_pls$qPCR_32)
    masque_numeric <- sapply(train_pls, is.numeric)
    train_pls_32 <- train_pls[,masque_numeric]
    train_pls_32$qPCR_32 <- train_pls$qPCR_32
    
    pls.model = plsr(qPCR_32 ~ ., data = train_pls_32, scale = TRUE, validation = "CV")
    
    # Prediction
    
    test_pls$qPCR_32 <- as.numeric(test_pls$qPCR_32)
    masque_numeric <- sapply(test_pls, is.numeric)
    test_pls <- test_pls[,masque_numeric]
    test_pls_32$qPCR_32 <- test_pls$qPCR_32
    
    pls.pred = predict(pls.model, test_pls_32, ncomp=100)
    
    plot( pls.pred,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
    abline( col="red")
    
    
    
    
    # Graph
    
    #load("Sauvegardes_objet_R.data/Essais_pls_model_data_global")
    
    pls1 <- plot(pls.model)
    
    pls.RMSEP = RMSEP(pls.model, estimate="CV")
    pls2 <- plot(pls.RMSEP, main="RMSEP PLS Solubility", xlab="components")
    
    pls3 <- validationplot(pls.model)
    
    
    # min_comp = which.min(pls.model$val)
    # points(min_comp, min(pls.model$val), pch=1, col="red", cex=1.5)
    
    
    coefficients = coef(pls.model)
    sum.coef = sum(sapply(coefficients, abs))
    coefficients = coefficients * 100 / sum.coef
    coefficients = sort(coefficients[, 1 , 1])
    barplot(tail(coefficients, 10))
    barplot(head(coefficients, 10))
    
    coef_pls <- data.frame(coefficients)
    coef_pls$Lambda <- as.numeric(substr(rownames(coef_pls), 2, 5))
    
    
      pls4 <- ggplot(data = coef_pls ) +
      aes(x = Lambda, y = coefficients) +
      #geom_point(colour = "#0c4c8a") +
      geom_text(label=rownames(coef_pls))+
      labs(title = "Coefficient du model PLS sur le training_set", subtitle = "A Ct32") +
      theme_minimal()
    
      pls4
      
      #save(list = 'pls.model', file = "Sauvegardes_objet_R.data/Essais_pls_model_data_global") 
    
      #esquisse::esquisser() 
      
     
    
## c) Package pls ####
  
    # Build the model on training set
    set.seed(123)
    model_pls <- train(
      qPCR_32 ~. , data = train_pls, method = "pls",
      scale = TRUE,
      trControl = trainControl("cv", number = 10),
      tuneLength = 10
    )
    # Plot model RMSE vs different values of components
    plot(model_pls)
    # Print the best tuning parameter ncomp that
    # minimize the cross-validation error, RMSE
    model_pls$bestTune
    
    
    
## d) Principal Component Regression ####
    
    
    # Build the model on training set
    set.seed(123)
    model_pls <- train(
      qPCR_32 ~. , data = train_pls, method = "pcr",
      scale = TRUE,
      trControl = trainControl("cv", number = 10),
      tuneLength = 10
    )
    # Plot model RMSE vs different values of components
    plot(model_pcr)
    # Print the best tuning parameter ncomp that
    # minimize the cross-validation error, RMSE
    model_pcr$bestTune
