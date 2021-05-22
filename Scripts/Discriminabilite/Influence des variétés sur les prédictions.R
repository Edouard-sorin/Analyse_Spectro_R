rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(party)  # Plot l'arbre de decision en Random Forest

#library(pals) # large palette de couleur


library(ggplot2)
library(ggparty)
library(ggdark)
library(ggpubr) # pour la fonction ggarrange qui permet de coller 2 graphiques

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")


# Discriminabilite des varietes via l'anova ####

# masque logique pour ne prendre que les arbres postitif a CT36 

data_long_var <- data_long_Ed
masque_qPCR <- data_long_var$qPCR_36 == 1
data_long_var <- data_long_var[masque_qPCR,]

data_long.lambda <- split(data_long_var, data_long_var$lambda)

model.lm <- lapply(data_long.lambda, function(x) lm(reflectance ~ code_variete ,data = x))

model.anov<-lapply(model.lm, function(x) anova(x))

# Pour annalyser l'anova en tant que tel, chercher dans l'anova pour la lougueur d'onde souhaité

intermed.anov <- as.data.frame(do.call(rbind, model.anov))

intermed.anov2 <- intermed.anov[,(which(colnames(intermed.anov) == "F value" ):which(colnames(intermed.anov) == "Pr(>F)"))]

Anov_Discrim <- intermed.anov2[ !is.na(intermed.anov2),]

rownames(Anov_Discrim ) = 0:4301

Anov_Discrim  <- Anov_Discrim [(which(rownames(Anov_Discrim ) == 0 ):which(rownames(Anov_Discrim ) == 2150)),]

colnames(Anov_Discrim) <- c("F_value","P_value")

Anov_Discrim $lambda <- as.numeric(350:2500)

## Graphique de présentation ####

graph_discrim <- ggplot(Anov_Discrim) +
  
  geom_line(aes(x = lambda , y = F_value)) +
  
  labs(x = "Longueur d'onde (en nm)", y = "F_value", title = "Discriminabilite a qPCR_36 via l'anova", subtitle = "Mesure via la F_value", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
  annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 9+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 2+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) 

graph_discrim

ggsave("Graphiques/Graph_explo_donnees/ Discriminabilite des variete positives à qPCR_36 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)

## graphP_value ####

graph_Pvalue<- ggplot(Anov_Discrim) +
  
  geom_line(aes(x = lambda , y = P_value)) +
  
  labs(x = "Longueur d'onde (en nm)", y = "P_value", title = "Discriminabilite a qPCR_36 via l'anova", subtitle = "Mesure via la P_value", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
  annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$P_value))], y = 0.02+ Anov_Discrim$P_value[(which.max(Anov_Discrim$P_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$P_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$P_value))], y = 0.002+Anov_Discrim$P_value[(which.max(Anov_Discrim$P_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) 


graph_Pvalue



## Zoom sur une lougueur d'onde recherche ####

summary(Anov_Discrim)

a <- Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))] # Valeur maximal du F_value
a <-  which(Anov_Discrim$lambda == 804)  # longueur d'onde à chercher
a <- Anov_Discrim$lambda[(which.min(Anov_Discrim$F_value))]

plot(model.lm[[a]])
res<- residuals(model.lm[[a]])
hist(res)

## Representation graphique des varietes avec ggparty ####

# Changer le nom des colonnes de X350 a 350nm

T1 <- Sys.time()

data_SPIR_var <- data_SPIR_Ed
masque_qPCR <- data_SPIR_var$qPCR_36 == 1
data_SPIR_var <- data_SPIR_var[masque_qPCR,]

RF_ct <- data_SPIR_var
RF_ct <- round(RF_ct[, grep("^X", names(RF_ct))],3)
masque_numeric <- sapply(RF_ct, is.numeric)
RF_ct <- RF_ct[,masque_numeric]
RF_ct$code_variete <- data_SPIR_var[,(which(colnames(data_SPIR_var) == "code_variete" ))]

gtree_var <- ctree( code_variete ~ . , data=RF_ct)


#autoplot(gtree_glob)

is.ggplot(ggparty(gtree_var)) # Pour savoir si le package party rend en compte le randomForest 

# ,layout = data.frame(id = c(4, 5, 8, 10),x = c(0.35, 0.15, 0.7, 0.8), y = c(0.95, 0.6, 0.8, 0.55)))


pdf("Graphiques/Graph_RF/Arbre de decision par variete.pdf",  width = 30, height = 15)


rf.var <- ggparty(gtree_var, terminal_space = 0.2 ) +  # terminal_space =  taille aloué aux plots (terminal)
  
  geom_edge(colour = "grey77", size = 0.5 ,  ) +
  
  
  geom_edge_label(colour = "Black", size = 4, alpha = 0.1) +
  
  geom_node_label(# map color of complete label to splitvar
    mapping = aes(),
    # map content to label for each line
    line_list = list(aes(label = splitvar),
                     aes(label = paste("p =",
                                       formatC(p.value,
                                               format = "e",
                                               digits = 2))),
                     aes(label = ""),
                     aes(label = paste0("Noeud ", id," N = ", round(nodesize/6,2)))
    ),
    # set graphical parameters for each line in same order
    line_gpar = list(list(size = 12),
                     list(size = 8),
                     list(size = 6),
                     list(size = 7,
                          col = "black",
                          fontface = "bold",
                          alignment = "left")
    ),
    # only inner nodes
    ids = "inner") +
  
  geom_node_label(aes(label = paste0("Nd ", id," N = ", round(nodesize/6,2))),
                  fontface = "bold",
                  ids = "terminal",
                  size = 2.7,
                  nudge_y = 0.025)+
  
  geom_node_plot(gglist = list(geom_bar(aes(x = "" , fill = code_variete),
                                        position = position_fill(),
                                        alpha = 0.8),
                               theme_classic(base_size = 15) ,
                               scale_fill_brewer(palette = "Dark2"),xlab("")),
                 shared_axis_labels = TRUE,
                 legend_separator = FALSE,
                 size = 1.2,
                 height = 1) +
  
  #scale_color_brewer(palette = "Spectral) +
  
  labs(title = "    Arbre de decision sur la reflectance des données global par rapport au différents vairétés positives à ct36", subtitle = "      formule : code_variete ~ reflectance | taux d'erreu out of bag = 25 % ") 

rf.var


dev.off()



T2 <- Sys.time()

difftime(T2,T1)




g0 <- ggplot(data_long_Ed) +
  aes(x = lambda, y = reflectance, group = qPCR_32, color = qPCR_32 ) +
  stat_summary(fun = mean, geom = "line" , size = 0.5) +
  scale_color_brewer(palette = "Dark2", labels = c("Négatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif à Ct<32)", color = "Résultat du test HLB à Ct<32") +
  dark_theme_gray() 
#theme_gray()

g0

## graph moyennés à CT var tout spectres ####

data_long_var <- data_long_Ed
data_SPIR_var <- data_SPIR_Ed
masque_qPCR <- data_SPIR_var$qPCR_36 == 0
data_SPIR_var <- data_SPIR_var[masque_qPCR,]

g0 <- ggplot(data_long_var) +
  aes(x = lambda, y = reflectance, group = code_variete, color = code_variete ) +
  stat_summary(fun = mean, geom = "line" , size = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen en fonction de la variété ", subtitle = "arbre positif à Ct<36)", color = "Résultat du test HLB à Ct<var") +
  dark_theme_gray() +
  theme(legend.position = "bottom")
#theme_gray()

g0

## graph moyennés à CT var zoomé sur les spectres 400 à 680 nm ####

g1 <- ggplot(data_long_var[data_long_var$lambda >= 400 & data_long_var$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = code_variete )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_brewer(palette = "Dark2")+
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen de 400 à 680 nm en fonction de la variété ")+
       #, subtitle = "arbre positif à Ct<36)", color = "Résultat du test HLB à Ct<var") +
  dark_theme_gray() +
  theme(legend.position = "none")
#theme_gray()
g1

## graph moyennés à CT var zoomé sur les spectres 700 à 1400 nm ####

g2 <- ggplot(data_long_var[data_long_var$lambda >= 700 & data_long_var$lambda <= 1400,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = code_variete )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen de 700 à 1400 nm en fonction de la variété ")+
       #, subtitle = "arbre positif à Ct<36)", color = "Résultat du test HLB à Ct<var") +
  dark_theme_gray() +
  theme(legend.position = "none")
#theme_gray()
g2


T1 <- Sys.time()

g_var <- ggarrange( g1, g2,g0, ncol = 3)

T2 <- Sys.time()

TG1 <- difftime(T2,T1)

g_var

ggsave("Graphiques/Graph_explo_donnees/Variété/Spectres moyen des feuilles pour les variétés sur l'ensemble des arbres négatifs.pdf",plot = g_var, units = "cm", width = 20, height = 10, scale = 2)

