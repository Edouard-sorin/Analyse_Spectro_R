rm(list=ls())  # nettoyage des listes de l'environnement de travail

library(tidyr) # pivot longer

library(ggplot2) #ggplot pour graphiques
library(ggpubr) # pour la fonction ggarrange qui permet de coller 2 graphiques

# Importation des donnees SPIR Hoarau ####

SPIR_Po <- read.table(file = "Donnees/Data_SPIR_Edouard/SPIR_Pothin_var.csv"
                      , header = T
                      , sep = ";"
                      , stringsAsFactors = T
                      , row.names = 1
                      , na.strings = c("","NA")
                      , dec = "," )


# SPIR_Po

code_labo <- rownames(SPIR_Po)

names(SPIR_Po) [1] = c("code_variete")
SPIR_Po$code_nbr_rep <- factor(substr(code_labo, 8, 8))
SPIR_Po$code_ech_arbre <- factor(substr(code_labo, 1, 2))
SPIR_Po$code_ech_feuille <- factor(substr(code_labo, 1, 3))
SPIR_Po$code_rep_feuille <- factor(paste(SPIR_Po$code_ech_feuille, SPIR_Po$code_nbr_rep,sep="")) 


rm (code_labo)

SPIR_Po <- SPIR_Po[ !is.na(SPIR_Po$X350),]

# Importation et preparation des resultats de la Qpcr Hoarau ####

Qpcr_Po <- read.table(file = "Donnees/Resultats_qPCR_Ed/qPCR_Pothin_19022021.csv"
                      
                      , header = T
                      , sep = ";"
                      , stringsAsFactors = T
                      , row.names = 1
                      , na.strings = c("","NA")
                      , dec = "," )

# on stocke les noms d'échantillons positifs selon nos 2 seuils de Ct (cycle de qPCR ), a savoir moins de 32 cyclces et moins de 36 cycles qPCR
seuils <- c(32, 36)

truep_Po <- lapply(seuils, function(x) unique(Qpcr_Po$Sample.Name[ which(Qpcr_Po$C..Mean < x & Qpcr_Po$C..SD < 1) ]))
names(truep_Po) <- paste("seuil", seuils, sep = ".")          


SPIR_Po[paste("qPCR", seuils, sep = "_")] <- lapply(truep_Po, function(x)
  as.numeric(SPIR_Po$code_ech_arbre %in% x ) ) # on cherche quels sont les code_ech_arbre qui se trouvent dans le vecteur x, x reprenant automatiquement les noms des arbres positifs selon le seuil choisi

select.lambda <- grep("^X", names(SPIR_Po))
SPIR_Po <- SPIR_Po[,c(names(SPIR_Po)[-select.lambda], names(SPIR_Po)[select.lambda] )]

rm (Qpcr_Po,select.lambda,seuils)

SPIR_Po[c("qPCR_32", "qPCR_36")] <- lapply(SPIR_Po[c("qPCR_32", "qPCR_36")], factor)


# Format_long ####

select.lambda <- grep("^X", names(SPIR_Po))
data_long_Po <- pivot_longer( data = SPIR_Po, cols = select.lambda, values_to = "reflectance", names_to = "lambda"  ) 



data_long_Po$lambda <- as.numeric(gsub("X", "", data_long_Po$lambda))

summary(data_long_Po)

data_long_Po <- data_long_Po[ !is.na(data_long_Po$reflectance),]

# Sauvegarde des donnees Hoarau en Rdata ####

save(data_long_Po, SPIR_Po, truep_Po,  file = "Sauvegardes_objet_R.data/SPIR_Po.Rdata")

# Exploration du jeu de donnees     ####

length(levels(SPIR_Po$code_ech_arbre)) # correspondant à 2x HLB 3x lot et 7x arbres soit 42

ftable (code_nbr_rep ~ code_ech_arbre , data = SPIR_Po) # repartition global pour voir les erreur de prise de données

ftable ( qPCR_32 ~code_ech_arbre, data = SPIR_Po) 

ftable ( qPCR_36 ~ code_ech_arbre, data = SPIR_Po) 

ftable ( code_variete ~ code_ech_arbre, data = SPIR_Po)

# graphiques d'exploration ####
## graph moyennés à CT 32 tout spectres ####

ggplot(data_long_Po) +
  aes(x = lambda, y = reflectance, group = qPCR_32, color = qPCR_32 ) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_brewer(palette = "Dark2", labels = c("Négatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Réflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif à Ct<32)", color = "Résultat du test HLB à Ct<32") +
  #dark_theme_gray() 
  theme_gray()

## graph moyennés à CT 32 zoomé sur les spectres 400 à 680 nm ####

ggplot(data_long_Po[data_long_Po$lambda >= 400 & data_long_Po$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_32 )+
  # geom_line() +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  theme_gray()

## graph moyennés à CT 36 tout spectres ####

ggplot(data_long_Po) +
  aes(x = lambda, y = reflectance, group = qPCR_36, color = factor(qPCR_36) ) +
  stat_summary(fun = mean, geom = "line") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif à Ct<36)", color = "Resultat du test HLB à Ct<36") +
  #dark_theme_gray() 
  theme_gray()

# graph moyennés à CT 36 zoomé sur les spectres 400 à 680 nm #####

ggplot(data_long_Po[data_long_Po$lambda >= 400 & data_long_Po$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = factor(qPCR_36) )+
  # geom_line() +
  stat_summary(fun = mean, geom = "line") +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  theme_gray()

# multiple graph feuilles moyennés à CT 32 ####

p1 <- ggplot(data_long_Po[data_long_Po$qPCR_32 == 1,]) +
  aes(x = lambda, y = reflectance, color = code_variete , group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  T) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyen des feuilles par arbres et par variété", subtitle = "arbre positif à Ct<32)", color = "Variété") +
  #dark_theme_gray() 
  theme_gray()



p0 <- ggplot(data_long_Po[data_long_Po$qPCR_32 == 0,]) +
  aes(x = lambda, y = reflectance, color = code_variete, group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  F) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyen des feuilles par arbres et par variété", subtitle = "arbre négatif à Ct<32)") +
  #dark_theme_gray() 
  theme_gray()

px <- ggarrange(p0, p1, ncol = 2)


#ggsave("Graphiques/Spectres moyenn des feuilles par arbres, par variété et par statut à ct32.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)
