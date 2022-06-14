##Tarea6 Ej 2 y 4
#####################Paquetería
rm(list =ls())
library(ggplot2)   
library(tseries)      
library(GGally)
library(latex2exp)
library(plot3D)
library(dplyr)
library(Matrix)
########################Ejercicios###############################


#######################Funciones necesarias######################
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

plotR<-function(x,y,GRID = 1,...){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}
seqR<-function(x,length)seq(x[1],x[2],length =length) 

histR<-function(x,y,GRID = 1,freq = 0,ylab = "Densidad",...){
hist(x,freq = freq,ylab = ylab,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}

Ces1<-function(A){
    A<-apply(A,2,function(x)(x - mean(x)))
    apply(A,2,function(x)x/sum(x^2)^(1/2))
    } 

Es<-function(A)apply(A,2,function(x)x/sum(x^2)^(1/2))

##Función creada por Imanol Nuñez, gracias por la ayuda amigo :)
VarDes<-function(A){
     svdX <- svd(x = A)
     V<- svdX$v
     Sigma<-svdX$d
     eta <- max(Sigma)/Sigma
     Phi_kj <-t(V^2)/Sigma^2
     Phi_k <- colSums(Phi_kj)
     Pi_jk <- t(t(Phi_kj)/Phi_k)
     RR<-cbind(Sigma,Pi_jk,eta)
     colnames(RR) <- c("Valor Singular", paste("var",1:ncol(X)),"Indice de Cond")
     RR<-as_tibble(RR,.name_repair = NULL)
     return(RR)
     }
###########################################################################

###Inciso 1)
##Corrlación

dat2<-read.table(file.choose(),header = TRUE, sep = ",")
X<-as.matrix(dat2[2:5])
(R<-cor(X))
ggpairs(dat2[,2:5],lower = list(continuous = my_fn))


dat2 <- data.frame(Ces1(as.matrix(dat2[,2:5])), y = dat2[ ,6] - mean(dat2[,6]))

lm2 <- lm(y ~ x1 + x2 + x3 + x4-1,data = dat2)
summary(lm2)
AIC(lm2)
#        Estimate Std. Error t value Pr(>|t|)  
#\alpha_1   31.607     14.308   2.209   0.0545 
#\alpha_2  27.500     36.784   0.748    0.4738  
#\alpha_3    2.261     15.788   0.143   0.8893  
#\alpha_4   -8.353     38.762  -0.215   0.8342 
#R-squared:  0.9824
#AIC: 63.837

DDD<-data.frame(dat2[1:4])
names(DDD) <- c("X1","X2","X3","X4") 
X<-as.matrix(DDD) 
round(cor(X),3)
#       X1     X2     X3     X4
#X1  1.000  0.229 -0.824 -0.245
#X2  0.229  1.000 -0.139 -0.973
#X3 -0.824 -0.139  1.000  0.030
#X4 -0.245 -0.973  0.030  1.000
heatmap(cor(X))

###Inciso 2) 
##VIFS 
round(diag(solve(R)),3)
# VIF_1  VIF_2    VIF_3  VIF_4 
# 38.496 254.423  46.868 282.513 

###Número de condición 

X1<-as.matrix(cbind(dat2[,1:4]))

#EigenValores
(vals<-sort(eigen(t(X1)%*%X1)$values))
##1.623\cdot10^{-3} 0.187 1.576 2.236

#Número de condición 
max(vals)/min(vals)

#Raíz cuadrada del número de cond
(max(vals)/min(vals))^(1/2)
#37.106

round(sort((max(vals)/vals)^(1/2)),3)
#VarianceDescomposition
(A<-VarDes(X1))
#   Valores Singulares V(\alpha_1)  V(\alpha_2) V(\alpha_3) V(\alpha_4) `Indice de Cond`
#           1.50   0.00263 0.000559 0.00148 0.000475             1.000   
#           1.26   0.00427 0.000427 0.00495 0.000457             1.191 
#           0.432  0.0635  0.00208  0.0465  0.000724             3.461
#           0.0403 0.930   0.997    0.947   0.998               37.106 


##Propuestas

##Reg con x3 y x4
lm22 <- lm(y ~ x4 + x3-1,data = dat2)
summary(lm22)
AIC(lm22)
#   Estimate Std. Error t value Pr(>|t|)    
#  -26.622      3.999  -6.658 3.57e-05 ***
#  -42.014      3.999 -10.507 4.50e-07 ***
#R^2: 0.935
#AIC: 76.745

##Reg sin x3 y x4 
lm24 <- lm(y ~ x1 + x2-1,data = dat2)
summary(lm24)
AIC(lm24)
#   Estimate Std. Error t value Pr(>|t|)    
#   29.920      2.357   12.70 6.51e-08 ***
#   35.698      2.357   15.15 1.03e-08 ***
#0.979
#62.312


####Mejor
##Reg lm24

X<-as.matrix(dat2[,1:2])
round(cor(X),3)
#      X1    X2
#X1 1.000 0.229
#X2 0.229 1.000

round(diag(solve(cor(X))),3)
#   VIF_1    VIF_2 
#   1.055 1.055 
(A<-VarDes(X))
#Numero de cond 
#1.26


(beta<-lm24$coef)

dat2<-read.table(file.choose(),header = TRUE, sep = ",")


head(dat2)
X4<-dat2[,2:3]
W<-dat2[,6]

means <- apply(X4,2,mean)
(Dest  <- apply(X4,2,function(x) x -mean(x)))
(Dest  <- apply(Dest,2,function(x) (sum(x^2))^(1/2)))
(Dy<-sum((W - mean(W))^2)^(1/2))

#Coeficientes para un modelo con las escalas originales para las X's
(beta1 <- beta*(1/Dest))
(beta1 <- c(mean(W) - sum(beta1*means),beta1))
round(beta1,3)
#           x1     x2 
#52.577  1.468  0.662


lm24 <- lm(y ~ x1 + x2,data = dat2)
summary(lm24)


