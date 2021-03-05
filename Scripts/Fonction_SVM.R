

fct_svm <- function(nb.rep, seuil.ct, list.feuilles) {
  
  maliste <- lapply(list.feuilles,function(feuille){ 
    if (nb.rep < nrow(feuille))  list_svm <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
    else list_svm <- feuille
    # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
    code <- grep("^code_ech", names(list_svm))
    qPCR <- grep("^qPCR_", names(list_svm))
    sortie <- cbind.data.frame(unique(list_svm[c(code,qPCR)])
                               , matrix(apply(list_svm[-(which(colnames(list_svm) == "code_nbr_rep"):which(colnames(list_svm) == "qPCR_36"))], 2, mean)
                                        , nr = 1, dimnames = list(NULL, names(list_svm)[-(which(colnames(list_svm) == "code_nbr_rep"):which(colnames(list_svm) == "qPCR_36"))])))
    sortie
  })
  
  test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 
  
  decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
  training_set <- test_ed[decoup,]
  test_set <- test_ed[!decoup,] 
  
  res.svm <- svm(y = training_set[,paste0("qPCR_", seuil.ct)]               # ici on prend les seuil donc sois 32 sois 36 , avec "paste0" colle sans separateurs
                 , x = training_set[, grep("^X", names(training_set))]      #  ici on prend ttes les longueurs d'ondes
                 , type = 'C-classification'
                 , kernel = 'linear'
  ) 
  
  svm_pred <-  predict(res.svm, newdata = test_set[,grep("^X", names(test_set))])
  
  svm.confusion <- confusionMatrix(data = svm_pred, reference = test_set[,paste0("qPCR_", seuil.ct)])   
  sortie <- c(svm.confusion$overall[1] , svm.confusion$byClass[1] , svm.confusion$byClass[3] )           # combine les 3 parametres recherches 
  names(sortie) <- c("Accuracy", "Sensitivity", "Precision")                                             # renomer
  sortie
} 
