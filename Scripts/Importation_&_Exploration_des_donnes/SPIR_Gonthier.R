rm(list=ls())  # nettoyage des listes de l'environnement de travail

library(tidyr) # pivot longer & pivot_wider

library(ggplot2) #ggplot pour graphiques
library(ggdark) # graphique noir
library(ggpubr) # pour la fonction ggarrange qui permet de coller 2 graphiques

# Importation des donnees SPIR Global ####

SPIR_Go <- read.table(file = "Donnees/Data_SPIR_Edouard/SPIR_Gonthier_var.csv"
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )


# Creation des colonne "code" de SPIR_Go

code_labo <- rownames(SPIR_Go)

names(SPIR_Go) [1] = c("code_variete")
names(SPIR_Go) [2] = c("code_agri")
SPIR_Go$code_nbr_rep <- factor(substr(code_labo, 8, 8))
SPIR_Go$code_ech_arbre <- factor(substr(code_labo, 1, 2))
SPIR_Go$code_ech_feuille <- factor(substr(code_labo, 1, 3))
SPIR_Go$code_rep_feuille <- factor(paste(SPIR_Go$code_ech_feuille, SPIR_Go$code_nbr_rep,sep="")) 


rm (code_labo)

SPIR_Go <- SPIR_Go[ !is.na(SPIR_Go$X350),]

# Importation et preparation des resultats de la Qpcr Globaux ####

Qpcr_Go <- read.table(file = "Donnees/Resultats_qPCR_Ed/qPCR_Gonthier_31032021.csv"
                           
                           , header = T
                           , sep = ";"
                           , stringsAsFactors = T
                           , row.names = 1
                           , na.strings = c("","NA")
                           , dec = "," )

# on stocke les noms d'echantillons positifs selon nos 2 seuils de Ct (cycle de qPCR ), a savoir moins de 32 cyclces et moins de 36 cycles qPCR
seuils <- c(32, 36)

truep_Go <- lapply(seuils, function(x) unique(Qpcr_Go$Sample.Name[ which(Qpcr_Go$C..Mean < x & Qpcr_Go$C..SD < 1) ]))
names(truep_Go) <- paste("seuil", seuils, sep = ".")          


SPIR_Go[paste("qPCR", seuils, sep = "_")] <- lapply(truep_Go, function(x)
  as.numeric(SPIR_Go$code_ech_arbre %in% x ) ) # on cherche quels sont les code_ech_arbre qui se trouvent dans le vecteur x, x reprenant automatiquement les noms des arbres positifs selon le seuil choisi

select.lambda <- grep("^X", names(SPIR_Go))
SPIR_Go <- SPIR_Go[,c(names(SPIR_Go)[-select.lambda], names(SPIR_Go)[select.lambda] )]

rm (Qpcr_Go,select.lambda,seuils)

SPIR_Go[c("qPCR_32", "qPCR_36")] <- lapply(SPIR_Go[c("qPCR_32", "qPCR_36")], factor)


# Format_long ####

select.lambda <- grep("^X", names(SPIR_Go))
data_long_Go <- pivot_longer( data = SPIR_Go, cols = select.lambda, values_to = "reflectance", names_to = "lambda"  ) 



data_long_Go$lambda <- as.numeric(gsub("X", "", data_long_Go$lambda))

summary(data_long_Go)

data_long_Go <- data_long_Go[ !is.na(data_long_Go$reflectance),]

# Sauvegarde des donnees en Rdata ####

save(data_long_Go, SPIR_Go, truep_Go,  file = "Sauvegardes_objet_R.data/Jeux de donnee/SPIR_Go.Rdata")

# Exploration du jeu de donnees     ####

length(levels(SPIR_Go$code_ech_arbre)) # correspondant a 2x HLB 3x lot et 7x arbres soit 42

ftable (code_nbr_rep ~ code_ech_arbre , data = SPIR_Go) # repartition global pour voir les erreur de prise de donnees

ftable ( qPCR_32 ~ code_ech_arbre, data = SPIR_Go) 

ftable ( qPCR_36 ~ code_ech_arbre, data = SPIR_Go) 

ftable ( code_variete ~ code_ech_arbre, data = SPIR_Go)

sum_var_36 <- ftable ( code_variete + qPCR_36 ~ code_ech_arbre, data = SPIR_Go)

sum_var_36 <- as.matrix(sum_var_36)

#write.table(x = var_36  , file = "Donnees/var_36.csv" , sep = ';')

totalt_Citron_0 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),1])/60

totalt_Citron_1 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),2])/60

totalt_Tangor_0 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),3])/60

totalt_Tangor_1 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),4])/60

totalt_Zanzibar_0 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),5])/60

totalt_Zanzibar_1 <- sum(sum_var_36[1:length(levels(SPIR_Go$code_ech_arbre)),6])/60

Total_all_Var <- rbind(totalt_Citron_0,totalt_Citron_1,totalt_Tangor_0,totalt_Tangor_1,totalt_Zanzibar_0,totalt_Zanzibar_1 )

Total_all_Var

# graphiques d'exploration ####
## graph moyennes des varietes ####

g1 <- ggplot(data_long_Go) +
  aes(x = lambda, y = reflectance, group = code_variete, color = factor(code_variete) ) +
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  #scale_color_manual(values = c('0' = "darkgreen", '1' = "brown4")) 
  scale_colour_viridis_d()+
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction des varietes", color = "Varietes") +
  dark_theme_gray() +
  theme(legend.position = "bottom")
#theme_gray()
g1



# graph moyennes des varietes zoome sur les spectres 400 a 680 nm #####

g2 <- ggplot(data_long_Go[data_long_Go$lambda >= 400 & data_long_Go$lambda <= 680,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = code_variete )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_colour_viridis_d()+
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction des varietes" )+
  dark_theme_gray() +
  theme(legend.position = "none")
#theme_gray()
g2





## graph moyennes des varietes zoome sur les spectres 700 a 1400 nm  ####

g3 <- ggplot(data_long_Go[data_long_Go$lambda >= 700 & data_long_Go$lambda <= 1400,] ) +
  aes(x = lambda, y = reflectance, group = code_ech_arbre, color = code_variete )+
  stat_summary(fun = mean, geom = "line", size = 0.5) +
  scale_colour_viridis_d()+
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectre moyen en fonction des varietes") +
  #, subtitle = "arbre positif a Ct<36)", color = "Resultat du test HLB a Ct<36"
  dark_theme_gray()  +
  theme(legend.position = "none")
#theme_gray()
g3

T1 <- Sys.time()

g_var <- ggarrange( g1, g2 ,g3, ncol = 3)

T2 <- Sys.time()

TG1 <- difftime(T2,T1)






ggsave("Graphiques/Graph_explo_donnees/Spectres_moyen_feuille_var_Gonthier.pdf",plot = g_var, units = "cm", width = 30, height = 20, scale = 2)


# multiple graph feuilles moyennes a CT 32 ####

p1 <- ggplot(data_long_Go[data_long_Go$qPCR_36 == 1,]) +
  aes(x = lambda, y = reflectance, color = code_variete , group = code_ech_feuille) +
  stat_summary(fun = mean, geom = "line", show.legend =  T) +
  scale_y_continuous(limits = c(0,1.2)) +
  #scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(code_ech_arbre)) +
  labs(x = "Longueur d'onde (en nm)", y = "Reflectance moyenne", title = "Spectres moyen des feuilles par arbres et par variete", subtitle = "arbre positif a Ct<36)", color = "Variete") +
  #dark_theme_gray() 
  theme_gray()



p0 <- ggplot(data_long_Go[data_long_Go$qPCR_36 == 0,]) +
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

#save(p0,p1  , file = "Graphiques/graphGeneral_SPIR_Go.Rdata")

# Calcul des valeurs represente graphiquement

# Graph par variete ####

Resum.spect<-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Go, mean) 
names(Resum.spect)[5] <- "mean_ref"
Resum.spect$Ref.sd <-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Go, sd)$reflectance

#Resum.spect$Ref.median <-  aggregate(reflectance ~ code_variete + code_ech_arbre + qPCR_36 + lambda, data_long_Go, median)$reflectance

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

Spec.by.var <-  aggregate(reflectance ~ code_variete  + qPCR_36 + lambda, data_long_Go, mean) 
names(Spec.by.var)[4] <- "mean_ref"
Spec.by.var$Ref.sd <-  aggregate(reflectance ~ code_variete  + qPCR_36 + lambda, data_long_Go, sd)$reflectance

ggplot(Spec.by.var) +
  aes(x = lambda, y = mean_ref) +
  geom_point() +
  geom_ribbon(aes(x = lambda, y = mean_ref,
                  ymin = mean_ref - Ref.sd,
                  ymax = mean_ref + Ref.sd))+
  scale_color_viridis_d(option = "viridis") +
  theme_minimal() +
  facet_wrap(vars(code_variete))



