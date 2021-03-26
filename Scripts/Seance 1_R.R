rm(list=ls())


#####################################################################
#                                                     Edouard sorin #
#                     Initiation R                                  #
#                                                     02/11/2020    #
#####################################################################

# Objectifs: 

# Prise en mains de R ####

# pour connaitre la version de R : R.Version() ou  R.version

# Les différente entités contenu dans R

" - Fonction

  - Vecteur : suite d'léléments de même nature
     
            - numérique
            - text
            - logique (True/False) R transforme automatiquement T ne 1 et F en 0 quand on trnaspose en numérique
  
  Quand on index ces vecteurs :
  
                              - Vecteur de positions (numérique)
                              - Vecteur de noms (texte)
                              - Vecteur de sélection (logique)
  
  - Facteur : est composé de 2 vecteurs , un vecteur d'observations (au format numérique correspondant au numéro de la modalité) (ex : espece) et d'un vecteur de modalités ( au format texte mais ou chacun a un numéro attribué) ( levels)
           
  - Matrice : vecteur à 2 dimensions. ATTENTION ! 1 Matrice = 1 Nature ( numérique, texte, logique)
  
  - Liste   : suite d'éléments de nature différentes
  
  - Data_frame: cas particulier de liste qui peut se gérer comme une matrice. Dans data_frame un tirroire = une variable = colonne. ATTENTION ! Contrainte sur le nombre d'observation dans chaque colonne. Ils doivent faire le même nombre pour être compatible et s'assembler dans un data frame.

"

# CITATION ####
# Pour Cit? :
# 
# citation("nom_package")

# Exemple
# citation("cluster")

#To cite the R package 'cluster' in publications use:
  
#  Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2019).  cluster: Cluster Analysis Basics and Extensions. R package version 2.1.0.

#A BibTeX entry for LaTeX users is

#@Manual{,
#  title = {cluster: Cluster Analysis Basics and Extensions},
#  author = {Martin Maechler and Peter Rousseeuw and Anja Struyf and Mia Hubert and Kurt Hornik},
#  year = {2019},
#  note = {R package version 2.1.0 --- For new features, see the 'Changelog' file (in the package source)},
#}




# #Pour afficher titre mettre 4 # -> ####
# #ou ctrl + shift + R dans "code, insert section


# FONCTION UTILE ####

# pour RUN  Ctrl+ Entrée

# # pour retour en arriere "ctrl + Z "

# # pour mettre en commentaire "ctrl + shift + c"

#   pour trouver & remplacer Ctrl + shift + J ou aller dans Edit

# pour afficher aide de la fonction cliquer sur  F1 qd on est sur une fonction

# pour lancer une entité "Ctrl + Entre" ou "; "nom de l'entité"  ex : compil <- c(5:1, 1:5) ;compil , R va afficher "compil"

# sur une fonction cliquer sur tabulation qd curser entre les 2 paranthèse pour avoir les proposition d'arguments de la fonction concerné


# pour la fonction sequence , la limite (pour cet ordi ) est : seq.default(from = 0,to = 1,by = 10^-9)

# Pour libérer mémoir RAM : gc ()

# Pour mettre " <- " cliquer sur "ALt" + "-"

# Tirage de DE de 100 :
D<-sample(1:100,1)
D

#DE truqué sur le 1 :

s<-sample(x = 100, size = 300, replace = T, prob = c(25,rep(1,99)))
s

table(s)

# code loi NORMAL

# 
# p <- seq(-3,3,0.1)
# p
# 
# toto <- dnorm (p)
# toto
# 
# g <- ggplot()+ aes(x = p, y = toto) + geom_line()
# g

# Symbol ####

#  ET = &

# OU = | 

# RECUPERE (indexation) = []  on met des vecteurs de position, de nom et de sélection

# "-" que sur des vecteurs numériques ( donne l'opposé )

# "!" que sur des vecteurs logiques (donne l'opposé )

# is.  pose une question sur la nature de l'objet

# as.  transforme la nature de l'objet

# avant une fonction catégorie de la fonction : : la fonction, car plusieurs fonctions ont les mêmes nom ex :
# base ::length 

# MODULES(PACKAGES) ####

# pour connaitre les modules actifs faire : search()

# pour aller chercher un module Package/Install liste déroulante des packages ou install.packages() pour installer depuis R

# pour aller voir les modules existant -> https://ftp.igh.cnrs.fr/pub/CRAN/         + Task views/


"# -x = L'opposé de x (numérique)
 # 1/x = Linverse de x  ()
 # !x = Le complémentaire de x  (logique)"

# Vous créez un nouveau répertoire de stockage des modules par défaut (solution intelligente).
# Voici comment le signifier à R  :
#   1. Créez un nouveau dossier à la racine de votre ordinateur ( ://C ou ://D) que vous pouvez
# appeler “Bibliothèque R” et copiez son chemin absolu (par exemple : “ ://C/Bibliothèque
#                                                       R”)
# 2. Dans Panneau de configuration > Système et sécurité > Système > Paramètres système
# avancés > Variable d’environnement, cliquez sur “Nouvelle... ” qui se trouve dans l’onglet
# Variable système, entrez comme nom : “R_LIBS” et comme valeur : le chemin de votre
# dossier.
# 3. Redémarrez RStudio. Dorénavant, lorsque vous installerez des nouveaux modules, ils seront stockés dans ce dossier quelle que soit votre version de R. Vous n’y penserez même
# plus !
  

# CREE UN PROJET R dans un répertoire de travail ####

# touche en haut à gauche, (petit carré bleu avec R dedans) 

# prendre la 2eme option, créer dans un répertoir de travail existant et organisé au préalable

# l'interet, qd on importe des données via "read.table" , pour la premier argument "file" qui permet d'ouvrir le fichier, on peut aller chercher le fichier en cliquant sur tabulation pour aller le chercher


# LA FONCTION APPLY ####

# apply(x,M,Fq,...)  # x = argument de apply correspondant à la matrice, M margin (1 pour ligne, 2 pour colonne) ,Fq = une fonction, ... argument de la fonction. -> en sortie soit (list, vecteur ou matrice) R choisi automatiquement

# lapply (S, Fq, ...) # S = une suite d'élément ne meme nature ou non, Fq = une fonction, ... argument de la fonction.
#la sortie est une liste

# sapply (S, Fq, ...) # S = une suite d'élément ne meme nature ou non,Fq = une fonction, ... argument de la fonction.
#la sortie est un (vecteur, matrice ou une liste)

"pour transfomer des variables "text" en facteurs dans les données "data"qui sont data.frame et donc une liste"
#
# a <- sapply(data,is.character)           # sapply génère un vecteur et non une liste

# pour transformer les characters selectionné en facteurs

# data[a] <- lapply (data[a], factor)     # factor = as.factor # lapply car notre donné d'entré est une liste donc on veut une liste ! d'ou lapply

# data[a] au début et data[a] dans l'attribution car on envoie ce qu'on à modifié ou on l'a pris !
#


# UTILISATION D'UN MASQUE LOGIQUE ####

# masque_logique <- Commode$tirroire == "ce_que_je_cherche"   # vecteur logique 
#
# Commode$tirroire[masque_logique] <-  "ce_qui_remplace"

# Instalation d'une nouvelle version de R et re instalation des packages 

# toto <- installed.packages(lib.loc = "c:/Bibliothèque R/")
# rownames(toto)

# pour connaitre cheminement des fichier faire :  getwd()



# pour convertir un fichier HTML R Markdown en PDF : open in browser / Clic droit / imprimer / imprimante/ Enregistrer en format pdf

# Conseil Rmarkdown ####

# pour ouvrir une zone script sur R Markdown : Ctrl + Alt + i

# pour le mode ouvrage :

# - le script markdown "index" doit toujours s'appeler index, il est le script qui regroupe tous les paramétrages effectués sur l'ouvrage

# pour pouvoir utiliser le "build" (en haut à droite avec "environment" et "Git" par exemple) qui permet de mettre en forme l'ouvrage il faut creer un projet R avec uniquement des Rmarkdown

# l

## a) conseils ####

# ne pas mettre les 4 #### sur Rmd 

# commenter en "pourquoi j'ai écrit ça" et pas "comment j'ai écrit ça" 

#  pas de points d'exlamation

# virgule à la fin (commme + dans ggplot) sans espace

# remettre les codes propres , aller dans Code/reformat code ou ctrl + shift + A et/ou ctrl + I pour les lignes. ATTENTION NE PAS FAIRE SUR LE SCRIPT PRICIPALE . 

# eval=FALSE pour ne pas executer le code

# cache=TRUE met les resultat dans un cache

# width.cutoff=50 pour avoir Max 50 caracteres

# install.packages("bookdown")

# install.packages("remotes") pour lire remotes::install_github("EcoFoG/EcoFoG", build_vignettes = TRUE)

# pour pouvoir lire le package "EcoFog"

# Aller ensuite dans "cree_nouveau script/RMarkdown/From_Template/Choisir_Ouvrage pour ecrir un rapport.

# Etape 1 Installer et maitriser GitHub

# Etape 2 Installer et maitriser Package EcoFog + bookdown

# Etape 3 Ecrire rapport sur ce support + pour bibliographie
#
# Aller sur Zotero/Slectionner biblio/Exporter bibliographie/ ofrmat Bibtex / Exporter les notes/ Enregistrer les notes dans le projet R/ Lire la bibliographie avec le @article{,
#  title = {cluster: Cluster Analysis Basics and Extensions},
#  author = {Martin Maechler and Peter Rousseeuw and Anja Struyf and Mia Hubert and Kurt Hornik},
#  year = {2019},
#  note = {R package version 2.1.0 --- For new features, see the 'Changelog' file (in the package source)},
#}
#
#
# Voir 

# https://ecofog.github.io/EcoFoG/memo/docs/memo.html
# 
# https://ericmarcon.github.io/Rochebrune2018/
#   
# https://ericmarcon.github.io/travailleR/chap-rediger.html
# 
# https://ecofog.github.io/EcoFoG/





## b) EcoFoG ####

# installer le package EcoFoG: remotes::install_github("EcoFoG/EcoFoG")

# # File/New Files/Rmarkdown/From template/Memo EcoFoG ou aller voir sur  :
#
# (https://ericmarcon.github.io/travailleR/chap-rediger.html#memo-article-bookdown) 

# dans un dossier vide. Commence par un Memo (section 4.3.1) à titre
# d'exercice. Déplacer les fichiers dans le dossier du projet synchronisé
# avec GitHub et le tour est joué : il reste à rédiger
# 
# 
## c) Git Hub ####

  # - pour utiliser GitHub en tant que développeur, il faut :
  #     * ouvrir un compte sur GitHub
  #     * installer git sur son ordinateur
  # Ensuite, pour créer un projet, la méthode la plus simple :

# https://ericmarcon.github.io/travailleR/chap-git.html#cr%C3%A9er-un-nouve%20#%20au-d%C3%A9p%C3%B4t%20,%20section%203.2.4

# pour push sur Github 

# git remote add origin https://github.com/Edouard-sorin/Analyse_Spectro_R.git
# 
# git branch -M master
# 
# git push -u origin master

# pas de modif sur github que sur Rstudio  ! 



# Depuis Gitignore on peut signifier a git ce qu'on ne veut pas rendre sous cette notation :  Documents/




# Essaie de créer un premier projet de mémo et le rendre lisible sur GitHub
# à la façon de https://github.com/EricMarcon/Inference-bayesienne 

# Pour ne pas avoir de fichier "unstaged" dans la fenetre git en haut à droite il  faut ouvrir .gitignore dans les dossier ou se trouve le projet R et noté tel quel le nom des dossier de la fenetre git avec / à la fin



