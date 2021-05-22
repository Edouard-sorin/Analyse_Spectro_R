
rm(list=ls())  # nettoyage des listes de l'environnement de travail

library(tidyr) # pivot longer & pivot_wider

library(ggplot2) #ggplot pour graphiques
library(ggdark) # graphique noir
library(ggpubr) # pour la fonction ggarrange qui permet de coller 2 graphiques

# Importation des donnees SPIR Global ####

data_SPIR_Ed <- read.table(file = "Donnees/Data_SPIR_Edouard/SPIR_Global.csv"
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )


# Creation des colonne "code" de data_SPIR_Ed

code_labo <- rownames(data_SPIR_Ed)

names(data_SPIR_Ed) [1] = c("code_variete")
names(data_SPIR_Ed) [2] = c("code_agri")
data_SPIR_Ed$code_nbr_rep <- factor(substr(code_labo, 8, 8))
data_SPIR_Ed$code_ech_arbre <- factor(substr(code_labo, 1, 2))
data_SPIR_Ed$code_ech_feuille <- factor(substr(code_labo, 1, 3))
data_SPIR_Ed$code_rep_feuille <- factor(paste(data_SPIR_Ed$code_ech_feuille, data_SPIR_Ed$code_nbr_rep,sep="")) 


rm (code_labo)

data_SPIR_Ed <- data_SPIR_Ed[ !is.na(data_SPIR_Ed$X350),]

# Importation et preparation des resultats de la Qpcr Globaux ####

data_Qpcr_Ed <- read.table(file = "Donnees/Resultats_qPCR_Ed/qPCR_global_Ed.csv"
                           
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )

# on stocke les noms d'echantillons positifs selon nos 2 seuils de Ct (cycle de qPCR ), a savoir moins de 32 cyclces et moins de 36 cycles qPCR
seuils <- c(32, 36)

trueP <- lapply(seuils, function(x) unique(data_Qpcr_Ed$Sample.Name[ which(data_Qpcr_Ed$C..Mean < x & data_Qpcr_Ed$C..SD < 1) ]))
names(trueP) <- paste("seuil", seuils, sep = ".")          


data_SPIR_Ed[paste("qPCR", seuils, sep = "_")] <- lapply(trueP, function(x)
  as.numeric(data_SPIR_Ed$code_ech_arbre %in% x ) ) # on cherche quels sont les code_ech_arbre qui se trouvent dans le vecteur x, x reprenant automatiquement les noms des arbres positifs selon le seuil choisi

select.lambda <- grep("^X", names(data_SPIR_Ed))
data_SPIR_Ed <- data_SPIR_Ed[,c(names(data_SPIR_Ed)[-select.lambda], names(data_SPIR_Ed)[select.lambda] )]

rm (data_Qpcr_Ed,select.lambda,seuils)

data_SPIR_Ed[c("qPCR_32", "qPCR_36")] <- lapply(data_SPIR_Ed[c("qPCR_32", "qPCR_36")], factor)

#save(trueP,  file = "Sauvegardes_objet_R.data/Jeux de donnee/trueP.Rdata")

# Format_long ####

select.lambda <- grep("^X", names(data_SPIR_Ed))
data_long_Ed <- pivot_longer( data = data_SPIR_Ed, cols = select.lambda, values_to = "reflectance", names_to = "lambda"  ) 



data_long_Ed$lambda <- as.numeric(gsub("X", "", data_long_Ed$lambda))

summary(data_long_Ed)

data_long_Ed <- data_long_Ed[ !is.na(data_long_Ed$reflectance),]

# Sauvegarde des donnees en Rdata ####

#save(data_long_Ed, data_SPIR_Ed, trueP,  file = "Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Exploration du jeu de donnees     ####

length(levels(data_SPIR_Ed$code_ech_arbre)) # correspondant a 2x HLB 3x lot et 7x arbres soit 42

ftable (code_nbr_rep ~ code_ech_arbre , data = data_SPIR_Ed) # repartition global pour voir les erreur de prise de donnees

ftable ( qPCR_32 ~ code_ech_arbre, data = data_SPIR_Ed) 

ftable ( qPCR_36 ~ code_ech_arbre, data = data_SPIR_Ed) 

ftable ( code_variete ~ code_ech_arbre, data = data_SPIR_Ed)

sum_var_36 <- ftable ( code_variete + qPCR_36 ~ code_ech_arbre, data = data_SPIR_Ed)

sum_var_36 <- as.matrix(sum_var_36)

#write.table(x = var_36  , file = "Donnees/var_36.csv" , sep = ';')

totalt_Citron_0 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),1])/60

totalt_Citron_1 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),2])/60

totalt_Tangor_0 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),3])/60

totalt_Tangor_1 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),4])/60

totalt_Zanzibar_0 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),5])/60

totalt_Zanzibar_1 <- sum(sum_var_36[1:length(levels(data_SPIR_Ed$code_ech_arbre)),6])/60

Total_all_Var <- rbind(totalt_Citron_0,totalt_Citron_1,totalt_Tangor_0,totalt_Tangor_1,totalt_Zanzibar_0,totalt_Zanzibar_1 )

Total_all_Var

# graphiques d'exploration ####
## graph moyennes a CT 32 tout spectres ####


T1 <- Sys.time()

g0 <- ggplot(data_long_Ed) +
  aes(x = lambda, y = reflectance, group = qPCR_32, color = qPCR_32 ) +
  stat_summary(fun = mean, geom = "line" , size = 0.5) +
  #scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4"))
  scale_color_brewer(palette = "Dark2", labels = c("Negatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif a Ct<32)", color = "Resultat du test HLB a Ct<32") +
  dark_theme_gray() +
  theme(legend.position = "bottom")
  #theme_gray()

g0

## graph moyennes a CT 32 zoome sur les spectres 400 a 680 nm ####

g1 <- ggplot(data_long_Ed[data_long_Ed$lambda >= 400 & data_long_Ed$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_32 )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen de 400 a 680 nm en fonction du statut HLB des arbres ") +
  # , subtitle = "arbre positif a Ct<32", color = "Resultat du test HLB a Ct<32"
  dark_theme_gray()  +
  theme(legend.position = "none")
  #theme_gray()
g1

## graph moyennes a CT 32 zoome sur les spectres 700 a 1400 nm ####

g2 <- ggplot(data_long_Ed[data_long_Ed$lambda >= 700 & data_long_Ed$lambda <= 1400,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_32 )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen de 700 a 1400 nm en fonction du statut HLB des arbres ") +
  #, subtitle = "arbre positif a Ct<32", color = "Resultat du test HLB a Ct<32"
  dark_theme_gray()  +
  theme(legend.position = "none")
#theme_gray()
g2

#ggsave("Graphiques/Graph_explo_donnees/Spectres_moyen_de 700 a 1400 nm_ct32.pdf",plot = g2, units = "cm", width = 30, height = 20, scale = 2)


g_32 <- ggarrange( g1, g2,g0, ncol = 3)

T2 <- Sys.time()

TG1 <- difftime(T2,T1)
TG1

## graph moyennes a CT 36 tout spectres ####

T1 <- Sys.time()

g3 <- ggplot(data_long_Ed) +
  aes(x = lambda, y = reflectance, group = qPCR_36, color = factor(qPCR_36) ) +
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  #scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) 
  scale_color_brewer(palette = "Dark2", labels = c("Negatif","Positif")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction du statut HLB des arbres ", subtitle = "arbre positif a Ct<36", color = "Resultat du test HLB a Ct<36") +
  dark_theme_gray() +
    theme(legend.position = "bottom")
  #theme_gray()
g3



# graph moyennes a CT 36 zoome sur les spectres 400 a 680 nm #####

g4 <- ggplot(data_long_Ed[data_long_Ed$lambda >= 400 & data_long_Ed$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_36 )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen de 400 a 680 nm en fonction du statut HLB des arbres " )+
  dark_theme_gray() +
theme(legend.position = "none")
  #theme_gray()
g4





## graph moyennes a CT 36 zoome sur les spectres 700 a 1400 nm  ####

g5 <- ggplot(data_long_Ed[data_long_Ed$lambda >= 700 & data_long_Ed$lambda <= 1400,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = qPCR_36 )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen de 700 a 1400 nm en fonction du statut HLB des arbres ") +
  #, subtitle = "arbre positif a Ct<36)", color = "Resultat du test HLB a Ct<36"
  dark_theme_gray()  +
  theme(legend.position = "none")
#theme_gray()
g5


g_36 <- ggarrange( g4, g5,g3, ncol = 3)
 
T2 <- Sys.time()

TG2 <- difftime(T2,T1)
TG2

T1 <- Sys.time()

g_glob <- ggarrange(g_32, g_36, nrow = 2)

T2 <- Sys.time()

TG3 <- difftime(T2,T1)
TG3



ggsave("Graphiques/Graph_explo_donnees/Spectres_moyen_feuille_ct32&ct36.pdf",plot = g_glob, units = "cm", width = 30, height = 20, scale = 2)


# multiple graph feuilles moyennes a CT 32 ####

p1 <- ggplot(data_long_Ed[data_long_Ed$qPCR_36 == 1,]) +
  aes(x = lambda, y = reflectance, color = code_variete , group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  T) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyen des feuilles par arbres et par variete", subtitle = "arbre positif a Ct<36)", color = "Variete") +
  #dark_theme_gray() 
  theme_gray()



p0 <- ggplot(data_long_Ed[data_long_Ed$qPCR_36 == 0,]) +
  aes(x = lambda, y = reflectance, color = code_variete, group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  F) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyen des feuilles par arbres et par variete", subtitle = "arbre negatif a Ct<36)") +
  #dark_theme_gray() 
  theme_gray()

T1 <- Sys.time()

px <- ggarrange(p0, p1, ncol = 2)

T2 <- Sys.time()

difftime(T2,T1)


#ggsave("Graphiques/Graph_explo_donnees/Spectres_moyen_feuille_var_negatif_ct36.pdf",plot = last_plot(), units = "cm", width = 20, height = 15, scale = 2)

#save(p0,p1  , file = "Graphiques/graphGeneral_data_SPIR_Ed.Rdata")

# Calcul des valeurs represente graphiquement

# Graph par variete ####

Resum.spect<-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Ed, mean) 
names(Resum.spect)[5] <- "mean_ref"
Resum.spect$Ref.sd <-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Ed, sd)$reflectance

#Resum.spect$Ref.median <-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Ed, median)$reflectance

  ggplot(Resum.spect) +
  aes(x = lambda, y = mean_ref) +
  geom_point() +
    geom_ribbon(aes(x = lambda, y = mean_ref,
                    ymin = mean_ref - Ref.sd,
                    ymax = mean_ref + Ref.sd))+
  scale_color_viridis_d(option = "viridis") +
  theme_minimal() +
    facet_wrap(vars(code_ech_arbre))
  
# Spectre par Variete
  
  Spec.by.var <-  aggregate(reflectance ~ code_variete  + qPCR_36 + lambda, data_long_Ed, mean) 
  names(Spec.by.var)[4] <- "mean_ref"
  Spec.by.var$Ref.sd <-  aggregate(reflectance ~ code_variete  + qPCR_36 + lambda, data_long_Ed, sd)$reflectance
  
  ggplot(Spec.by.var) +
    aes(x = lambda, y = mean_ref) +
    geom_point() +
    geom_ribbon(aes(x = lambda, y = mean_ref,
                    ymin = mean_ref - Ref.sd,
                    ymax = mean_ref + Ref.sd))+
    scale_color_viridis_d(option = "viridis") +
    theme_minimal() +
    facet_wrap(vars(code_variete))



