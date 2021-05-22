rm(list=ls())  # nettoyage des listes de l'environnement de travail

# Library ####


library(effects) # Pour utiliser la fonction allEffects

library(forestmodel) # Pour utiliser la fonction forest_model


library(gtsummary) # Pour utiliser la fonction


library(ggplot2)
library(ggdark)
library(GGally) # Pour utiliser la fonction ggcoef_model
library(ggeffects) # Pour la fonction ggeffect


library(nnet) # Pour utiliser la fonction multinom


# Importation du jeu de donnee Global

load("Sauvegardes_objet_R.data/Jeux de donnee/data_SPIR_Ed.Rdata")

# Importation du jeu Barret

#load("Sauvegardes_objet_R.data/Jeux de donnee/SPIR_Ba.Rdata")

# Regression Logistique Binaire pour qPCR36 ~ tout les lambdas #### 

# voir http://larmarange.github.io/analyse-R/regression-logistique.html#r%C3%A9gression-logistique-binaire

data_long_Ed <- data_long_Ba


reg <- glm( qPCR_36 ~ 
              
              X680	+
              X681	+
              X682	+
              X683	+
              X684	+
              X685	+
              X686	+
              X687	+
              X688	+
              X689	+
              X690	+
              X691	+
              X692	+
              X693	+
              X694	+
              X695	+
              X696	+
              X697	+
              X698	+
              X699	+
              X700	+
              X701	+
              X702	+
              X703	+
              X704	+
              X705	+
              X706	+
              X707	+
              X708	+
              X709	+
              X710	+
              X711	+
              X712	+
              X713	+
              X714	+
              X715	+
              X716	+
              X717	+
              X718	+
              X719	+
              X720	+
              X721	+
              X722	+
              X723	+
              X724	+
              X725	+
              X726	+
              X727	+
              X728	+
              X729	+
              X730	+
              X731	+
              X732	+
              X733	+
              X734	+
              X735	+
              X736	+
              X737	+
              X738	+
              X739	+
              X740	+
              X741	+
              X742	+
              X743	+
              X744	+
              X745	+
              X746	+
              X747	+
              X748	+
              X749	+
              X750	+
              X751	+
              X752	+
              X753	+
              X754	+
              X755	+
              X756	+
              X757	+
              X758	+
              X759	+
              X760	+
              X761	+
              X762	+
              X763	+
              X764	+
              X765	+
              X766	+
              X767	+
              X768	+
              X769	+
              X770	+
              X771	+
              X772	+
              X773	+
              X774	+
              X775	+
              X776	+
              X777	+
              X778	+
              X779	+
              X780	+
              X781	+
              X782	+
              X783	+
              X784	+
              X785	+
              X786	+
              X787	+
              X788	+
              X789	+
              X790	+
              X791	+
              X792	+
              X793	+
              X794	+
              X795	+
              X796	+
              X797	+
              X798	+
              X799	+
              X800	
              
              
             
            
            , data = data_SPIR_Ed,  family = binomial(logit)) # La variable d’interet doit être binaire (0,1) pour utiliser la condition family = binomial(logit)

reg

#summary(reg)

# graph forest

gx <- forest_model(reg)

pdf("Graphiques/Graph_RF/forest_model_qPCR36_680-800nm.pdf",  width = 30, height = 70 )

gx

dev.off()

#save(reg_all, gx,  file = "Sauvegardes_objet_R.data/glm_all.Rdata")


# tableau avec les parametres de la regression

t1 <- tbl_regression(reg_700, exponentiate = TRUE)
t1

# Regression Logistique Binaire pour qPCR36 ~ autres variables #### 

reg <- glm( qPCR_36 ~ code_variete + code_agri , data = data_SPIR_Ed,  family = binomial(logit)) # La variable d’interet doit être binaire (0,1) pour utiliser la condition family = binomial(logit)

reg

summary(reg)

t0 <- exp(cbind(coef(reg), confint(reg)))

#g0 <- plot(reg)
#g0

# tableau avec les parametres de la regression

t1 <- tbl_regression(reg, exponentiate = TRUE)
t1

# graph coef

#g1 <- ggcoef_model(reg, exponentiate = TRUE)
#g1

# graph forest

#pdf("Graphiques/Graph_RF/forest_model_qPCR36.pdf",  width = 30, height = 15)



g2 <- forest_model(reg)
g2

dev.off()

# graph ggeffect

g3 <- cowplot::plot_grid(plotlist = plot(ggeffect(reg)))
g3

# Regression Logistique Multinomiale ####

data_long_Ed <- data_long_Ba

intermed.arbre <- aggregate(reflectance ~ code_variete +  qPCR_32 + qPCR_36 + lambda, data =  data_long_Ed, mean)

#mean.arbre <- pivot_wider(intermed.arbre, names_from = "lambda", values_from = "reflectance", names_prefix = "X")


regm <- multinom(code_variete  ~ qPCR_36 + qPCR_32  + code_ech_arbre , data = SPIR_Ba) # , family = binomial(logit)


regm2 <- step(regm)

summary(regm2)

t0 <- odds.ratio(regm2)
t0 

t1 <- tbl_regression(regm2, exponentiate = TRUE)
t1

g0 <- ggcoef_multinom(
  regm2,
  exponentiate = TRUE
)
gà

g1 <- plot(allEffects(regm2))
g1

g2 <- plot(ggeffect(regm2, "code_agri"))
g2

