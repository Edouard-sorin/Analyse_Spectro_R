

fct_var <- function(nb.rep, list.feuilles) {
  
  maliste <- lapply(list.feuilles,function(feuille){ 
    list_svm <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
    # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
    code <- grep("^code_ech", names(list_svm))
    variete <- grep("code_var", names(list_svm))
    sortie <- cbind.data.frame(unique(list_svm[c(code,variete)])
                               , matrix(apply(list_svm[-(which(colnames(list_svm) == "code_variete"):which(colnames(list_svm) == "qPCR_36"))], 2, mean)
                                        , nr = 1, dimnames = list(NULL, names(list_svm)[-(which(colnames(list_svm) == "code_variete"):which(colnames(list_svm) == "qPCR_36"))])))
    sortie
  })
  
  test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 
  
  test_ed$var <- 1
  test_ed$var[test_ed$code_variete == "Tangor"  ] <- 2
  test_ed$var[test_ed$code_variete == "Zanzibar"] <- 3
  
  test_ed$var <- factor(test_ed$var)
  
  decoup <- sample.split(test_ed[,which(colnames(test_ed) == "var")], SplitRatio = 0.5) # on decoupe le jeu de donné en training set et test set
  training_set <- test_ed[decoup,]
  test_set <- test_ed[!decoup,] 

  # Prediction SVM ####
  
  model_SVM_var <- svm(y = training_set$var              
                      , x = training_set[, grep("^X", names(training_set))]      
                      , type = 'C-classification'
                      , kernel = 'linear'
  ) 
  
  # Prediction sur le test_set
  
  svm_pred <- predict(model_SVM_var ,newdata = test_set[, grep("^X", names(test_set))] , decision.values = T)
  
  svm.confusion <- confusionMatrix(data = svm_pred, reference = test_set$var)   
  
  sortie <- c(svm.confusion$overall[1] , svm.confusion$byClass[1] , svm.confusion$byClass[3] )           # combine les 3 parametres recherches 
  names(sortie) <- c("Accuracy", "Sensitivity", "Precision")                                             # renomer
  sortie
  
} 

