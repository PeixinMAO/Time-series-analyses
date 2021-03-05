#*** Définition de l'environnement ***

# On met notre répertoire de travail dans "wd"
wd <- "C:/Users/Asus/Desktop/ENSAE/2A/Séries Temporelles Linéaires/Projet/kafangnibo_pmao"
setwd(wd)
getwd()

#*** Librairies nécessaires pour le projet ***
require(readr)             # pour lire le fichier csv
require(dplyr)             # pour manipuler les données (trier, réordonner etc.)
require(zoo)               # pour le format séries temp
require(tseries)           # pour des fonctions diverses en séries temporelles
require(fUnitRoots)        # pour des tests de racine unitaire plus modulables
require(forecast)          # pour la prévision
require(car)               # pour tracer des ellipses

#*** Importation des données ***
datafile <- "Donnees_The_Cafe.csv"
data <- read.csv(datafile,sep=";")
View(data) # Pour regarder les données

#*** Tansformations préliminaires ***

data <- arrange(data,Période) # pour mettre les données dans l'ordre 1990 à 2020 et non 2020 à 1990
# comme c'est le cas initialement
View(data) # pour observer le changement
dim(data) # pour regarder les dimensions de notre base: Il y a 362 observations

plot(data) # Premier graphe : pas tout à fait bon. On va mettre sous format zoo
dates <- as.yearmon(seq(from=1990+0/12,to=2020+1/12,by=1/12)) #index des dates pour les données
donnees <- zoo(data$Indice, order.by = dates)
plot(donnees, main = "IPI-Transformation du Thé et du Café",
     xlab = "Dates", ylab = "Valeurs de l'indice") # C'est mieux

# Transformations
# La serie ne semble pas varier autour d'une tendance linéaire déterministe
# Elle ne semble pas stationnaire

#*** Transformations pour stationnarité ***

diffdonnees <- diff(donnees,1) # difference premiere

plot(cbind(donnees,diffdonnees), xlab= "Dates",
     ylab= c("Originale", "Différenciée"), main = NA)
# La série différenciée une fois semble stable autour d'une constante nulle: elle semble stationnaire

#* Tests de racines unitaires pour détermination/confirmation de l'ordre d'intégration *

# Vérification de l'existence d'une constante et/ou une tendance linéaire non nulle/ Vérification
# nécéssaire avant la conduite des tests de racines unitaires pour la spécification
summary(lm(donnees~dates)) # régression de la série sur les dates pour voir la présence d'une tendance
# linéaire

# Constante et coefficient(négatif) significatifs donc on fait le test ADF avec constante et 
# tendance non nulles

adf <- adfTest(donnees, lag=0, type="ct") # test ADF dans le cas avec constante et tendance

# Il faut vérifier la validité du test i.e l'absence d'autocorrélation des résidus (bruit blanc)
# Fonction pour le test de Ljung-BoX sur les résidus du test ADF
Qtests <- function(series, k, fitdf=0){
  pvals <- apply(matrix(1:k), 1, FUN=function(l){
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients)) # Tests jusqu'à l'ordre 24 (2 ans)

# l'absence d'autorrélation est rejetée pour Q(4) à Q(12), les résidus sont donc autocorrélés
# et le test ADF sans aucun retard est non valide. Il faut donc ajouter des retards.

#* Tests ADF jusqu'à l'obtention des résidus non autocorrélés
# Fonction pour déterminer le nombre de retards nécessaires pour obtenir des résidus non autocorrélés
adfTest_valid <- function(series,kmax,type){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF avec ",k, " retards: résidus OK ? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0){
      noautocorr <- 1; cat("OK \n")
    }
    else cat("Non \n")
    k <- k + 1
  }
  return(adf)
}

adf <- adfTest_valid(donnees,24,"ct")
# 11 retards ont été nécéssaires pour les résidus soient non autocorrélés

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients)) # On vérifie bien que
# le test de validité prouve que les résidus sont non-autocorrélés
# Ce qui est le cas

adf # pour commenter le résultat du test
# racine unitaire n'est pas rejetée: serie au moins integrée une fois

# Maintenant test de la présence de racine unitaire  sur la série différenciée une fois
summary(lm(diffdonnees~dates[-1])) # sans la première date car on a différencié la série
# Constante et coefficient non significatifs donc faire test ADF sans constante et sans tendance

adf <- adfTest_valid(diffdonnees,24, type="nc") # pour le nombre de retards à introduire
# 2 retards à introduire pour que les résidus soient non-corrélés. On vérifie cela :

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients)) # Toutes les pvalues sont > 5%

adf
# On rejette la racine unité : pvalue = 0,01<5% 
# La série "diffdonnees" est stationnaire et la série "donnees" est donc I(1)

#*** Modélisation ARMA ***

#* Etude des autocorrélations totales et partielles *
x <- diffdonnees-mean(diffdonnees) # centrer la série

# Choix des ordres maximaux de l'ARMA
par(mfrow=c(1,1))
acf(x, main = "Auto-corrélogramme 
    de la série transformée" ) # AcF

pacf(x, main = "Auto-corrélogrammne partiel
     de la série transformée") # PACF

# on regarde jusqu'à 2 ans de retard

par(mfrow=c(1,2))
acf(x, 24, main = "ACF de la série
    transformée" ) # AcF
pacf(x,24, main = "PACF de la série
     transformée") # PACF
# pmax = 11 avec PACF
# qmax = 12 avec l'ACF

pmax = 11 ; qmax = 12

#* Choix du modèle
# On estime tous les modèles possibles, on choisit ceux bien ajustés et valides 
# Puis on prend celui/ceux qui minimise un critère d'information

# Ecriture d'une fonction présentant les coefficients, leurs écarts-types et leurs p-values 
# pour un modèle donné
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

# Fonction estimant tous les modèles 
modelchoice <- function(p,q,BD=x, k=24){
  estim <- try(arima(BD, c(p,0,q)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

# Fonction affichant tous les modèles 
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax); t(apply(matrix(1:dim(pqs)[1]),1,function(row){
    p <- pqs[row,1]
    q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax,qmax) #estime tous les arimas

selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] # sélectionne les modèles bien ajustés et valides
selec

# On a 20 modèles bien ajustés et valides

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #crée
#une liste des ordres p et q des modeles candidats

names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")")                 # renomme les éléments de la liste

models <- lapply(pqs, function(pq) arima(x,c(pq[["p"]],0,pq[["q"]])))     # crée une liste des modeles candidats
BICs <- vapply(models, FUN.VALUE=numeric(1), function(m) c("BIC"=BIC(m))) # calcule les BIC des modèles candidats
AICs <- vapply(models, FUN.VALUE=numeric(1), function(m) c("AIC"=AIC(m))) # calcule les AIC des modèles candidats


BICs
min(BICs)

# Le modèle avec le plus petit BIC parmi les modèles bien ajustés et valides
# est l'ARIMA(1,1,2)
arma12 <- arima(x,c(1,0,2),include.mean=F)

# Le modèle avec le plus petit AIC parmi les modèles bien ajustés et valides
# est l'ARIMA(4,1,2)
arma42 <- arima(x,c(4,0,2),include.mean=F)
# Les deux modèles sont bien ajustés et valides, et minimisent chacun un des critères d'informations.
# Ils sont les meilleurs modèles pour l'instant.
# On va utiliser le critère du R2 ajusté pour déterminer le meilleur d'entre eux

# Fonction pour le calcul du R2 ajusté

adj_r2 <- function(model){
  p <- model$arma[1]                                    # recupere l'ordre AR
  q <- model$arma[2]                                    # recupere l'ordre MA
  n <- model$nobs-max(p,q)                              # taille de l'echantillon
  ss_res <- sum(model$residuals^2)                      # somme des residus au carree
  ss_tot <- sum(x[-c(1:max(p,q))]^2)                    # somme des observations de l'echantillon au carre
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))         # R2 ajusté
  return(adj_r2)
}

# Calcul des R2 ajustés
adj_r2(arma42)
adj_r2(arma12)

# L'ARIMA (4,1,2) a le R2 ajusté le plus élevé. C'est le meilleur modèle.

#*** Prévision
#** Meilleure prévision pour X_(t+1) et X_(t+2)
SeriePrédite<-forecast(arma42, h=2,level=c(95))

SeriePrédite # Le vecteur prédit
X1= -0.4776403 # valeur prédite de X_(T+1)
X2= 0.5604930 # valeur prédite de X_(T+2)
par(mfrow = c(1,1))
plot(SeriePrédite, main = "Prévision pour les dates t+1 et t+2")

# Représentation de la région de confiance de X_(t+1) et X_(t+2)


Variance_residus = var(residuals(arma42))   # la variance des résidus
Matrix_cov = matrix(nrow=2,ncol=2)          # la matrice de variance covariance du vecteur
Matrix_cov[1,1]= Variance_residus
Matrix_cov[2,1]= Variance_residus *(- arma42$coef[5])
Matrix_cov[1,2]= Variance_residus *(- arma42$coef[5])
Matrix_cov[2,2]= Variance_residus*(1+(arma42$coef[5])^2)

alpha=0.95                      # le niveau du confiance
Quantile = qchisq(alpha , df=2) # le quantile d'une Chi(2) à 2 degrés de liberté

Center = c(X1,X2)                                                     # le centre de l'ellipse
conf = ellipse(center=Center,shape=Matrix_cov, radius=sqrt(Quantile)) # détermine la région de confiance
plot(conf, type='l',col="red" , xlab="X_(T+1)", ylab = "X_(T+2)")     # plot  de la région de confiance
points(X1, X2, col="green")                                           # le vecteur prédit
title("Région de confiance ")

