
rm(list=ls())

# LIBRARY ####

library(caret) # package Machine learning # confusion matrix
library(caTools) # sample.split

library(ggplot2) # Package ggplot pour graphiques
library(ggdark) # Met un style de graphique ggplot en noir
library(ggpubr)# Utilisation de la fonction ggarrange qui permet de coller 2 graphiques

library(pls)
library(plsdepot) # Partial Least Square


# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/data_SPIR_Ed.Rdata")

# Importation du jeu de donnee de Hoarau

load("Sauvegardes_objet_R.data/SPIR_Ho.Rdata")

# Importation du jeu de donnee de Pothin

load("Sauvegardes_objet_R.data/SPIR_Po.Rdata")

# ------------------------------Partial Least Square a Ct<32------------------------------

# Preparation du training_pls et test_pls  ####

nb.rep <- 6
seuil.ct <- 32

Tirage <- split(data_SPIR_Ed[,-1], data_SPIR_Ed$code_ech_feuille, drop = T)


maliste <- lapply(Tirage,function(feuille){ 
  list_pls <- feuille[sample(1:length(feuille$code_ech_feuille), nb.rep),]  
  # on fait ça pour tte les feuilles mais tjrs en choisissant le nombre de rep tire aleatoirement
  code <- grep("^code_ech", names(list_pls))
  qPCR <- grep("^qPCR_", names(list_pls))
  sortie <- cbind.data.frame(unique(list_pls[c(code,qPCR)])
                             , matrix(apply(list_pls[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))], 2, mean)
                                      , nr = 1, dimnames = list(NULL, names(list_pls)[-(which(colnames(list_pls) == "code_nbr_rep"):which(colnames(list_pls) == "qPCR_36"))])))
  sortie
})

test_ed <- do.call(rbind, maliste) # on colle toute les listes créées précédement 

decoup <- sample.split(test_ed[,paste0("qPCR_", seuil.ct)], SplitRatio = 0.75) # on decoupe le jeu de donné en training set et test set
train_pls <- test_ed[decoup,]
test_pls <- test_ed[!decoup,] 

# Construction du Partial least square ####

# Identify the predictor columns

Pred <- grep("^X", names(train_rf))


# Identity the response column

Resp = "qPCR_32"

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
