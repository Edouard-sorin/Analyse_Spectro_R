rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####

library(party)  # Plot l'arbre de decision en Random Forest

#library(pals) # large palette de couleur


library(ggplot2)
library(ggparty)
library(ggdark)

# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")


# Discriminabilite du statut HLB  via l'anova ####

data_long.lambda <- split(data_long_Ed, data_long_Ed$lambda)

model.lm <- lapply(data_long.lambda, function(x) lm(reflectance ~ qPCR_32 ,data = x))

model.anov<-lapply(model.lm, function(x) anova(x))

# Pour annalyser l'anova en tant que tel, chercher dans l'anova pour la lougueur d'onde souhaité

intermed.anov <- as.data.frame(do.call(rbind, model.anov))

intermed.anov2 <- intermed.anov[,(which(colnames(intermed.anov) == "F value" ):which(colnames(intermed.anov) == "Pr(>F)"))]

Anov_Discrim <- intermed.anov2[ !is.na(intermed.anov2),]

rownames(Anov_Discrim ) = 0:4301

Anov_Discrim  <- Anov_Discrim [(which(rownames(Anov_Discrim ) == 0 ):which(rownames(Anov_Discrim ) == 2150)),]

colnames(Anov_Discrim) <- c("F_value","P_value")

Anov_Discrim $lambda <- as.numeric(350:2500)

# Discriminabilite des viretes via l'anova ####

data_long.lambda <- split(data_long_Ed, data_long_Ed$lambda)

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
  
  labs(x = "Longueur d'onde (en nm)", y = "F_value", title = "Discriminabilite a qPCR_32 via l'anova", subtitle = "Mesure via la F_value", color = "") +
  
  dark_theme_gray() +
  
  #scale_colour_viridis_d(option = "plasma") +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme(panel.grid.major.y = element_line(colour = "grey20")) +
  
annotate(geom = "text", x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 9+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))], label = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], size = 4, colour = "red")+
  
  geom_point(aes(x = Anov_Discrim$lambda[(which.max(Anov_Discrim$F_value))], y = 2+Anov_Discrim$F_value[(which.max(Anov_Discrim$F_value))] ), shape = 124, fill="darkred" ,color = "red" , size = 2) 
  
graph_discrim

#ggsave("Graphiques/Graph_explo_donnees/ Discriminabilite a qPCR_32 via l'anova.pdf", units = "cm", width = 20, height = 15, scale = 2)

## graphP_value ####

graph_Pvalue<- ggplot(Anov_Discrim) +
  
  geom_line(aes(x = lambda , y = P_value)) +
  
  labs(x = "Longueur d'onde (en nm)", y = "P_value", title = "Discriminabilite a qPCR_32 via l'anova", subtitle = "Mesure via la P_value", color = "") +
  
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

# Discriminabilite des varietes via Random Forest ####

seuil.ct = 32

RF_ct <- data_SPIR_Ed
masque_numeric <- sapply(RF_ct, is.numeric)
RF_ct <- RF_ct[,masque_numeric]
RF_ct[,paste0("qPCR_", seuil.ct)] <- data_SPIR_Ed[,paste0("qPCR_", seuil.ct)] 

gtree_glob <- ctree( qPCR_32 ~ . , data=RF_ct)

# plot 

pdf("Graphiques/Graph_RF/Arbre de decision a ct32.pdf",  width = 30, height = 15)

plot(gtree_glob ,)

#on ferme le graphique

dev.off()

print(gtree_glob)

## Representation graphique du statut HLB avec ggparty ####

#load("Sauvegardes_objet_R.data/Jeux de donnee/SPIR_Ba.Rdata")

seuil.ct = 32

RF_ct <- data_SPIR_Ed
RF_ct <- round(RF_ct[, grep("^X", names(RF_ct))],3)
masque_numeric <- sapply(RF_ct, is.numeric)
RF_ct <- RF_ct[,masque_numeric]
RF_ct[,paste0("qPCR_", seuil.ct)] <- data_SPIR_Ed[,paste0("qPCR_", seuil.ct)] 

gtree_glob <- ctree( qPCR_32 ~ . , data=RF_ct)


#autoplot(gtree_glob)

is.ggplot(ggparty(gtree_glob)) # Pour savoir si le package party rend en compte le randomForest 

# ,layout = data.frame(id = c(4, 5, 8, 10),x = c(0.35, 0.15, 0.7, 0.8), y = c(0.95, 0.6, 0.8, 0.55)))

rf.ct <- ggparty(gtree_glob, terminal_space = 0.2 ) +  # terminal_space =  taille aloué aux plots (terminal)
  
geom_edge(colour = "grey77", size = 0.5) + # mettre reflectance
  
 geom_edge_label(colour = c("300", "600", "633", "600" ,"633", "900", "999" ,"063" ,"066", "78", "79" ,"80", "81" ,"182", "83", "84", "85", "86", "87" ,"88", "89" ,"90", "91", "92" ,"93", "94", "95", "96" ,"97", "98", "99") , size = 4 , alpha = 0.1) +
  
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
                     list(size = 10),
                     list(size = 6),
                     list(size = 8,
                          col = "black",
                          fontface = "bold",
                          alignment = "left")
    ),
    # only inner nodes
    ids = "inner") +
  
  geom_node_label(aes(label = paste0("Nd ", id," N = ", round(nodesize/6,2)), size = 2),
                  fontface = "bold",
                  ids = "terminal",
                  size = 2.7,
                  nudge_y = 0.013)+

  geom_node_plot(gglist = list(geom_bar(aes(x = "" , fill = qPCR_32),
                                        position = position_fill(),
                                        alpha = 0.8),
                               theme_classic(base_size = 10) ,
                               scale_fill_brewer(palette = "Dark2"),xlab("")),
                 shared_axis_labels = TRUE, # Met 1 legende pour tous les plots
                 legend_separator = TRUE, # trait entre la legende et les plots
                 size = 1.2,
                 height = 1) +
  
  theme_void() +
  
  #scale_color_brewer(palette = "Spectral) +

labs(title = "    Arbre de decision sur la reflectance des données global par rapport au statut HLB", subtitle = "     Arbre positif à Ct<32 | formule : qPCR36 ~ reflectance | taux d'erreu out of bag = 25 % ") 
    
rf.ct

#ggsave("Graphiques/Graph_RF/Arbre de decision sur le jeu de donnée global a ct32 (party).pdf",plot = last_plot(), units = "cm", width = 50, height = 10, scale = 2)

# palette de couleur du package pals :

pal.bands(coolwarm, parula, ocean.haline, brewer.blues, cubicl, kovesi.rainbow, ocean.phase,alphabet, alphabet2, glasbey, kelly, polychrome, stepped, stepped2, stepped3, tol, watlington, stepped, main="Colormap suggestions")

# plus d'info sur https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html

# Plus d'info sur ggparty : https://cran.r-project.org/web/packages/ggparty/vignettes/ggparty-graphic-partying.html

#https://github.com/kwstat/pals/issues/3

ggplot(mtcars) + 
  geom_bar(aes(x=factor(hp), fill=factor(hp))) +
  scale_fill_manual(values=as.vector(polychrome(26)))


## Representation graphique des varietes avec ggparty ####

# Changer le nom des colonnes de X350 a 350nm

RF_ct <- data_SPIR_Ed
RF_ct <- round(RF_ct[, grep("^X", names(RF_ct))],3)
masque_numeric <- sapply(RF_ct, is.numeric)
RF_ct <- RF_ct[,masque_numeric]
RF_ct$code_variete <- data_SPIR_Ed[,(which(colnames(data_SPIR_Ed) == "code_variete" ))]

gtree_glob <- ctree( code_variete ~ . , data=RF_ct)


#autoplot(gtree_glob)

is.ggplot(ggparty(gtree_glob)) # Pour savoir si le package party rend en compte le randomForest 

# ,layout = data.frame(id = c(4, 5, 8, 10),x = c(0.35, 0.15, 0.7, 0.8), y = c(0.95, 0.6, 0.8, 0.55)))

rf.var <- ggparty(gtree_glob, terminal_space = 0.2 ) +  # terminal_space =  taille aloué aux plots (terminal)
  
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
                               theme_classic(base_size = 5) ,
                               scale_fill_brewer(palette = "Dark2"),xlab("")),
                 shared_axis_labels = TRUE,
                 legend_separator = FALSE,
                 size = 1.2,
                 height = 1) +
  
  #scale_color_brewer(palette = "Spectral) +
  
  labs(title = "    Arbre de decision sur la reflectance des données global par rapport au différents vairétés", subtitle = "      formule : code_variete ~ reflectance | taux d'erreu out of bag = 25 % ") 

rf.var
