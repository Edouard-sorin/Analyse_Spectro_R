rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(party)  # Plot l'arbre de decision en Random Forest

#library(pals) # large palette de couleur


library(ggplot2)
library(ggparty)
library(ggdark)

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")



# Discriminabilite de la relation statut/parcelle via l'anova ####

# masque logique pour ne prendre que les arbres négatifs a CT36 

data_long_agri <- data_long_Ed
#masque_qPCR <- data_long_agri$qPCR_36 == 0
#data_long_agri <- data_long_agri[masque_qPCR,]

data_long.lambda <- split(data_long_agri, data_long_agri$lambda)

model.lm <- lapply(data_long.lambda, function(x) lm(reflectance ~ qPCR_36 * code_agri ,data = x))

model.anov<-lapply(model.lm, function(x) anova(x))

# Pour annalyser l'anova en tant que tel, chercher dans l'anova pour la lougueur d'onde souhaité

intermed.anov <- as.data.frame(do.call(rbind, model.anov))

intermed.anov2 <- intermed.anov[,(which(colnames(intermed.anov) == "F value" ):which(colnames(intermed.anov) == "Pr(>F)"))]

Anov_Discrim <- intermed.anov2[ !is.na(intermed.anov2),]

nrow_anov <- nrow(Anov_Discrim)

rownames(Anov_Discrim)<- 0:(nrow_anov-1)

Anov_Discrim  <- Anov_Discrim [(which(rownames(Anov_Discrim) == 0 ):which(rownames(Anov_Discrim) == ((2151*3)-1))),]

colnames(Anov_Discrim) <- c("F_value","P_value")

Anov_Discrim$lambda <- as.numeric(rep (350:2500 , each = 3))

Anov_Discrim$Anov_by <- factor(rep(c("Influence statut","Influence parcelle","Influence statut*parcelle"),each=1))

#esquisse::esquisser() 

rm(data_long_agri,data_long.lambda,model.anov,model.lm,intermed.anov,intermed.anov2,trueP)

## Graphique Relation statut/parcelles ###

graph_discrim_agri <- ggplot(Anov_Discrim) +
  aes(x = lambda, y = F_value, colour = Anov_by, size = P_value) +
  geom_line() + 
  labs(x = "Longueur d'onde (en nm)", y = "F_value", title = "Relation statut/parcelle a qPCR_36 via l'anova", subtitle = "Mesure via la F_value", color = "") +
  
  dark_theme_gray() +
  #scale_colour_viridis_d( ) + #labels = c("Influence parcelle","Influence statut","Influence statut*parcelles")
  
  #scale_color_brewer(palette = "Paired") +
  
  scale_colour_manual(values = c("#D2B48C","#7EC0EE","#8B4500"))+
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
  annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 14+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 2+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) +
  
  theme(legend.position = "bottom") 

graph_discrim_agri


ggsave("Graphiques/Graph_relation_variables/Relation statut_parcelles a qPCR_36 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)

# Discriminabilite de la relation statut/variete via l'anova ####


data_long_agri <- data_long_Ed
#masque_qPCR <- data_long_agri$qPCR_36 == 0
#data_long_agri <- data_long_agri[masque_qPCR,]

data_long.lambda <- split(data_long_agri, data_long_agri$lambda)

model.lm <- lapply(data_long.lambda, function(x) lm(reflectance ~ qPCR_36 * code_variete ,data = x))

model.anov<-lapply(model.lm, function(x) anova(x))

# Pour annalyser l'anova en tant que tel, chercher dans l'anova pour la lougueur d'onde souhaité

intermed.anov <- as.data.frame(do.call(rbind, model.anov))

intermed.anov2 <- intermed.anov[,(which(colnames(intermed.anov) == "F value" ):which(colnames(intermed.anov) == "Pr(>F)"))]

Anov_Discrim <- intermed.anov2[ !is.na(intermed.anov2),]

nrow_anov <- nrow(Anov_Discrim)

rownames(Anov_Discrim)<- 0:(nrow_anov-1)

Anov_Discrim  <- Anov_Discrim [(which(rownames(Anov_Discrim) == 0 ):which(rownames(Anov_Discrim) == ((2151*3)-1))),]

colnames(Anov_Discrim) <- c("F_value","P_value")

Anov_Discrim$lambda <- as.numeric(rep (350:2500 , each = 3))

Anov_Discrim$Anov_by <- factor(rep(c("Influence statut","Influence variete","Influence statut*variete"),each=1))

#esquisse::esquisser() 

rm(data_long_agri,data_long.lambda,model.anov,model.lm,intermed.anov,intermed.anov2,trueP)

## Graphique Relation statut/variete ###
  
  graph_discrim_var <- ggplot(Anov_Discrim) +
  aes(x = lambda, y = F_value, colour = Anov_by, size = P_value) +
  geom_line() + 
    labs(x = "Longueur d'onde (en nm)", y = "F_value", title = "Relation statut/variete a qPCR_36 via l'anova", subtitle = "Mesure via la F_value", color = "") +
    
  dark_theme_gray() +
  #scale_colour_viridis_d( ) + #labels = c("Influence parcelle","Influence statut","Influence statut*parcelles")
  
  #scale_color_brewer(palette = "Paired") +
  
  scale_colour_manual(values = c("#7EC0EE","#191970","#CAFF70"))+
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
  annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 14+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 2+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) +
  
  theme(legend.position = "bottom") 

graph_discrim_var


#ggsave("Graphiques/Graph_relation_variables/Relation statut_varietes a qPCR_36 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)


# Discriminabilite de la relation statut/variete/parcelle via l'anova ####


data_long_agri <- data_long_Ed
#masque_qPCR <- data_long_agri$qPCR_36 == 0
#data_long_agri <- data_long_agri[masque_qPCR,]

data_long.lambda <- split(data_long_agri, data_long_agri$lambda)

model.lm <- lapply(data_long.lambda, function(x) lm(reflectance ~ qPCR_36 * code_variete * code_agri ,data = x))

model.anov<-lapply(model.lm, function(x) anova(x))

# Pour annalyser l'anova en tant que tel, chercher dans l'anova pour la lougueur d'onde souhaité

intermed.anov <- as.data.frame(do.call(rbind, model.anov))

intermed.anov2 <- intermed.anov[,(which(colnames(intermed.anov) == "F value" ):which(colnames(intermed.anov) == "Pr(>F)"))]

Anov_Discrim <- intermed.anov2[ !is.na(intermed.anov2),]

nrow_anov <- nrow(Anov_Discrim)

rownames(Anov_Discrim)<- 0:(nrow_anov-1)

Anov_Discrim  <- Anov_Discrim [(which(rownames(Anov_Discrim) == 0 ):which(rownames(Anov_Discrim) == ((2151*6)-1))),]

colnames(Anov_Discrim) <- c("F_value","P_value")

Anov_Discrim$lambda <- as.numeric(rep (350:2500 , each = 6))

Anov_Discrim$Anov_by <- factor(rep(c("Influence statut","Influence variete","Influence parcelle","Influence statut_variete","Influence statut_parcelle","Influence variete*parcelle"),each=1)) 

#esquisse::esquisser() 

rm(data_long_agri,data_long.lambda,model.anov,model.lm,intermed.anov,intermed.anov2,trueP)

select.influ_agri <- grep("Influence statut_parcelle", Anov_Discrim$Anov_by)

select.influ_var <- grep("Influence statut_v", Anov_Discrim$Anov_by)

Anov_Influ <- Anov_Discrim[c(select.influ_agri,select.influ_var),]

## Graphique Relation statut/parcelles ###

graph_discrim_var_agri <- ggplot(Anov_Discrim) +
  aes(x = lambda, y = F_value, colour = Anov_by, size = P_value) +
  geom_line() + 
  labs(x = "Longueur d'onde (en nm)", y = "F_value", title = "Relation statut/variete/parcelles a qPCR_36 via l'anova", subtitle = "Mesure via la F_value", color = "") +
  
  dark_theme_gray() +
  #scale_colour_viridis_d( ) + #labels = c("Influence parcelle","Influence statut","Influence statut*parcelles")
  
  #scale_color_brewer(palette = "Paired") +
  
  scale_colour_manual(values = c("#D2B48C","#7EC0EE","#8B4500","#191970","#CAFF70","#228B22"))+
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
  annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 14+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 2+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) +
  
  theme(legend.position = "bottom") 

graph_discrim_var_agri

#ggsave("Graphiques/Graph_relation_variables/Relation statut_variete_parcelles à qPCR_36 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)

# graph var&agri ####


g0 <- ggplot(data_long_Ed) +
  aes(x = lambda, y = reflectance, group = qPCR_36, color = qPCR_36 ) +
  stat_summary(fun = mean, geom = "line" , size = 0.5) +
  #scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4"))
  scale_color_brewer(palette = "Dark2", labels = c("Negatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif a Ct<36)", color = "Resultat du test HLB a Ct<36") +
  dark_theme_gray() +
  theme(legend.position = "bottom") +
  facet_wrap( facets = data_long_Ed$code_variete)

#theme_gray()

g0



g1 <- ggplot(data_long_Ed) +
  aes(x = lambda, y = reflectance, group = qPCR_36, color = qPCR_36 ) +
  stat_summary(fun = mean, geom = "line" , size = 0.5) +
  #scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4"))
  scale_color_brewer(palette = "Dark2", labels = c("Negatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif a Ct<36)", color = "Resultat du test HLB a Ct<36") +
  dark_theme_gray() +
  theme(legend.position = "bottom") +
  facet_wrap( facets = data_long_Ed$code_agri)

#theme_gray()

g1


#ggsave("Graphiques/GGraph_relation_variables/ Discriminabilite des variete positives à qPCR_36 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)
# Histogram sur une longueur d'onde ####

X703 <- data_SPIR_Ed[,c(which(colnames(data_SPIR_Ed) == "X703"),which(colnames(data_SPIR_Ed) == "code_variete"),which(colnames(data_SPIR_Ed) == "code_agri"),which(colnames(data_SPIR_Ed) == "qPCR_36"))]

ggplot(X703) +
  aes(x = X703, fill = qPCR_36, colour = code_variete) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(l = 50 , h = c(150,20)) + 
  xlim(0, 0.5) +
  scale_colour_manual(values = c("#00EEEE","#FFFF00","#0000CD"))+
  dark_theme_gray() +
  theme(legend.position = "bottom") +
  facet_wrap(vars(code_agri))

#esquisse::esquisser() 

#ggsave("Graphiques/Graph_relation_variables/Histogramme de la longueur d'onde 703.pdf", units = "cm", width = 20, height = 15, scale = 2)

