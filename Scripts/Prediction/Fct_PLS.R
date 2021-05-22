

fct_pls <- function(nb.rep, seuil.ct, list.feuilles) {
  
  maliste <- lapply(list.feuilles,function(feuille){   # On cree une list qui piochera au hazard 1 feuille avec le nombre de rep correspondant
    if (nb.rep < nrow(feuille))  list_pls <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
    else list_pls <- feuille
    # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
    code <- grep("^code_ech", names(list_pls))
    qPCR <- grep("^qPCR_", names(list_pls))
    sortie <- cbind.data.frame(unique(list_pls[c(code,qPCR)])
                               , matrix(apply(list_pls[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))], 2, mean) # on moyenne la valeur des rep pour chaque feuille
                                        , nr = 1, dimnames = list(NULL, names(list_pls)[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))])))
    sortie
  })
  
  test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 
  
  
  test_ed <- test_ed[-(which(colnames(test_ed) == "qPCR_32"))]
  
  test_ed[[paste0("qPCR_", seuil.ct)]] <- as.numeric(as.character(test_ed[[paste0("qPCR_", seuil.ct)]] ))
  
  decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
  training_set <- test_ed[decoup,]
  test_set <- test_ed[!decoup,] 
  
  res.pls <-  plsr(
    
    training_set[,paste0("qPCR_", seuil.ct)]  ~ . ,
    data = training_set[, grep("^X", names(training_set))], 
    scale = TRUE, 
    validation = "CV"
    
  )
  
  test_set[[paste0("qPCR_", seuil.ct)]] <- as.factor(test_set[[paste0("qPCR_", seuil.ct)]] )
  
  pls_pred <-  predict(res.pls, newdata = test_set[,grep("^X", names(test_set))] ,decision.values = T, ncomp=100)
  
  pls_pred[pls_pred < 0.35] = 0 
  pls_pred[pls_pred >= 0.35] = 1 
  
  pls_pred <- as.data.frame(pls_pred)
  colnames(pls_pred) <- "pls_pred"
  
  pls_pred <- as.factor(pls_pred$pls_pred)
  
  pls.confusion <- confusionMatrix(data = pls_pred, reference = test_set[,paste0("qPCR_", seuil.ct)])   
  sortie <- c(pls.confusion$overall[1] , pls.confusion$byClass[1] , pls.confusion$byClass[3] )           # combine les 3 parametres recherches 
  names(sortie) <- c("Accuracy", "Sensitivity", "Precision")                                             # renomer
  sortie
} 
