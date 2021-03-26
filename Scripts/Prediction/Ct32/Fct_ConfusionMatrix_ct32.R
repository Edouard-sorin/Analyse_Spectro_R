
fct_ConfusionMatrix <- function(nb.rep, seuil.ct, list.feuilles) {
  
  maliste <- lapply(list.feuilles,function(feuille){ 
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
  training_set <- test_ed[decoup,]
  test_set <- test_ed[!decoup,]  

  # Criteres choisi pour la separation des donnees
  
  crit.pos <- 0.4
  crit.neg <- 0.35
  
  trueP$seuil.36 <- NULL
  
  # Prediction RF ####
  
  model_rf_32 <- randomForest(
    x = training_set[, grep("^X", names(training_set))], 
    y = training_set[[paste0("qPCR_", seuil.ct)]],
    ntree = 100
  )
  
  # Prediction sur les feuilles de la base d'apprentissage
  
  rf_pred <- test_set
  rf_pred$rf_pred_32 <- predict(model_rf_32,newdata = test_set, decision.values = T)
  
  # Conversion des resultats de la prediction en numerique
  
  rf_pred$rf_pred_32 = as.numeric(as.character(rf_pred$rf_pred_32))
  
  # Moyennage des resultats de la prediction pour chaque arbres
  
  rf_pred_arbres_32 = aggregate(rf_pred_32 ~ code_ech_arbre, data = rf_pred, mean, na.rm = T)
  
  rf_pred_arbres_32$crit <- 0.5
  rf_pred_arbres_32$crit[rf_pred_arbres_32$rf_pred_32 >= crit.pos  ] <- 1
  rf_pred_arbres_32$crit[rf_pred_arbres_32$rf_pred_32 <= crit.neg  ] <- 0
  
  # Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite
  
  rf_pred_arbres_32[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
    as.numeric(rf_pred_arbres_32$code_ech_arbre %in% x ))
  
  # Résultats de la prédiction sous forme de matrice de confusion 
  
  rf_pred_arbres_32$crit <- factor(rf_pred_arbres_32$crit ,levels= c(0,1))
  
  rf_confusion_matrix_32 <- ftable( qPCR_32 ~ crit , data = rf_pred_arbres_32 )
  
  rf_confusion_matrix_32 <- as.matrix(rf_confusion_matrix_32)
  
  TP_rf <- rf_confusion_matrix_32[1,1]
  
  TN_rf <- rf_confusion_matrix_32[2,2]
  
  FN_rf <- rf_confusion_matrix_32[2,1]
  
  FP_rf <- rf_confusion_matrix_32[1,2]
  
  # Calcul des parametres rf de la matrice de confusion
  
  Accuracy_rf32 <- ((TP_rf+TN_rf) / (TP_rf+TN_rf+FN_rf+FP_rf)*100)
  
  Precision_rf32 <- ((TP_rf / (TP_rf+FP_rf)*100))
  
  Sensitivity_rf32 <- ((TP_rf / (TP_rf+FN_rf)*100))
  
  Parametre_rf_32 <- rbind(Accuracy_rf32,Precision_rf32,Sensitivity_rf32)
  
  # Prediction SVM ####
  
  model_SVM_32 <- svm(y = training_set[,paste0("qPCR_", seuil.ct)]              
                 , x = training_set[, grep("^X", names(training_set))]      
                 , type = 'C-classification'
                 , kernel = 'linear'
  ) 
  
  # Prediction sur le test_set
  
  svm_pred <- test_set
  svm_pred$svm_pred_32 <- predict(model_SVM_32 ,newdata = test_set[, grep("^X", names(test_set))] , decision.values = T)
  
  # Conversion des resultats de la prediction en numerique
  
  svm_pred$svm_pred_32 = as.numeric(as.character(svm_pred$svm_pred_32))
  
  # Moyennage des resultats de la prediction pour chaque arbres
  
  svm_pred_arbres_32 = aggregate(svm_pred_32 ~ code_ech_arbre, data = svm_pred, mean, na.rm = T)
  
  # Critere choisi sur les observations graphiques
  
  svm_pred_arbres_32$crit <- 0.5
  svm_pred_arbres_32$crit[svm_pred_arbres_32$svm_pred_32 >= crit.pos  ] <- 1
  svm_pred_arbres_32$crit[svm_pred_arbres_32$svm_pred_32 <= crit.neg  ] <- 0
  
  # Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite
  
  svm_pred_arbres_32[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
    as.numeric(svm_pred_arbres_32$code_ech_arbre %in% x ))
  
  # Résultats de la prédiction sous forme de matrice de confusion 
  
  svm_pred_arbres_32$crit <- factor(svm_pred_arbres_32$crit ,levels= c(0,1))
  
  svm_confusion_matrix_32 <- ftable( qPCR_32 ~ crit , data = svm_pred_arbres_32 )
  
  svm_confusion_matrix_32 <- as.matrix(svm_confusion_matrix_32)
  
  TP_svm <- svm_confusion_matrix_32[1,1]
  
  TN_svm <- svm_confusion_matrix_32[2,2]
  
  FN_svm <- svm_confusion_matrix_32[2,1]
  
  FP_svm <- svm_confusion_matrix_32[1,2]
  
  # Calcul des parametres svm de la matrice de confusion
  
  Accuracy_svm32 <- ((TP_svm+TN_svm) / (TP_svm+TN_svm+FN_svm+FP_svm)*100)
  
  Precision_svm32 <- ((TP_svm / (TP_svm+FP_svm)*100))
  
  Sensitivity_svm32 <- ((TP_svm / (TP_svm+FN_svm)*100))
  
  Parametre_svm_32 <- rbind(Accuracy_svm32,Precision_svm32,Sensitivity_svm32)

  # Prediction PLS ####
  
  
  model_pls_32 <-  plsr(
    
    training_set[[paste0("qPCR_", seuil.ct)]] ~ . ,
    data = training_set[, grep("^X", names(training_set))], 
    scale = TRUE, 
    validation = "CV"
    
  )
  
  # Prediction
  
  pls_pred <- test_set
  
  # Prediction sur les feuilles de la base d'apprentissage
  
  pls_pred$pls_pred_32 <- predict(model_pls_32, newdata = test_set[, grep("^X", names(test_set))], decision.values = T, ncomp=100)
  
  # Conversion des resultats de la prediction en numerique
  
  pls_pred$pls_pred_32 = as.numeric(as.character(pls_pred$pls_pred_32))
  
  # Moyennage des resultats de la prediction pour chaque arbres
  
  pls_pred_arbres_32 = aggregate(pls_pred_32 ~ code_ech_arbre, data = pls_pred, mean, na.rm = T)

  # Parametrage pour presentation graphique et la matrice de confusion
  
  pls_pred_arbres_32$crit <- 0.5
  pls_pred_arbres_32$crit[pls_pred_arbres_32$pls_pred_32 >= crit.pos  ] <- 1
  pls_pred_arbres_32$crit[pls_pred_arbres_32$pls_pred_32 <= crit.neg  ] <- 0
  
  # Creation d'une nouvelle colonne des resulats issu de la qPCR, pour comparer à la valeurs predite

  pls_pred_arbres_32[paste("qPCR", seuil.ct, sep = "_")] <- lapply(trueP, function(x)
    as.numeric(pls_pred_arbres_32$code_ech_arbre %in% x ))
  
  # Résultats de la prédiction sous forme de matrice de confusion 
  
  pls_pred_arbres_32$crit <- factor(pls_pred_arbres_32$crit ,levels= c(0,1))
  
  pls_confusion_matrix_32 <- ftable( qPCR_32 ~ crit , data = pls_pred_arbres_32 )
  
  pls_confusion_matrix_32 <- as.matrix(pls_confusion_matrix_32)
  
  TP_pls <- pls_confusion_matrix_32[1,1]
  
  TN_pls <- pls_confusion_matrix_32[2,2]
  
  FN_pls <- pls_confusion_matrix_32[2,1]
  
  FP_pls <- pls_confusion_matrix_32[1,2]
  
  # Calcul des parametres pls de la matrice de confusion
  
  Accuracy_pls32 <- ((TP_pls+TN_pls) / (TP_pls+TN_pls+FN_pls+FP_pls)*100)
  
  Precision_pls32 <- ((TP_pls / (TP_pls+FP_pls)*100))
  
  Sensitivity_pls32 <- ((TP_pls / (TP_pls+FN_pls)*100))
  
  Parametre_pls_32 <- rbind(Accuracy_pls32,Precision_pls32,Sensitivity_pls32)
  
  
  
  # Regroupement des predictions ####
  
#All_Confusion_matrix_32 <- rbind(rf_confusion_matrix_32,svm_confusion_matrix_32,pls_confusion_matrix_32)

#rownames(All_Confusion_matrix_32) <- c("Predicted Positive_RF","Predicted Negative_RF","Predicted Positive_SVM","Predicted Negative_SVM","Predicted Positive_PLS","Predicted Negative_PLS")

#colnames(All_Confusion_matrix_32) <- c("Actual Positive","Actual Negative")
  
parametre_32 <- rbind(Parametre_rf_32,Parametre_svm_32,Parametre_pls_32)
  
  All_parametre_32 <- as.vector(parametre_32)
  
  names(All_parametre_32) <- c(rownames(parametre_32))

  All_parametre_32

} 
