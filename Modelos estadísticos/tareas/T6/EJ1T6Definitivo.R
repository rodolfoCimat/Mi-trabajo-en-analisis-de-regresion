##Tarea6 Ej 1 y 3
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
    apply(A,2,function(x)x/(sum(x^2))^(1/2))
    } 

Es<-function(A)apply(A,2,function(x)x/(sum(x^2))^(1/2))

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
     colnames(RR) <- c("Valor Singular", paste("var",1:ncol(A)),"Indice de Cond")
     RR<-as_tibble(RR,.name_repair = NULL)
     return(RR)
     }
###########################################################################
###Diagramas de dispersión 

dat1<-read.table(file.choose(),header = TRUE, sep = ",")

ggpairs(dat1,lower = list(continuous = my_fn))

dat1 <- data.frame(Ces1(as.matrix(dat1[,1:5])), W = dat1[ ,6] - mean(dat1[,6]))
head(dat1)

par(mfrow=c(2,3))
####V1 vs W
regXY<-lm(W ~ V1, data = dat1)
summary(regXY)
newdat<-data.frame(V1 = seq(min(dat1$V1), max(dat1$V1), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V1,conf[,1],type = "l",col="red",
      xlab = TeX("V1"), ylab =TeX("W"), 
      main = TeX("V1 vs W"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V1,dat1$W,pch = 20,cex = 1.3)
matlines(newdat$V1, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V1, pred[,2:3], col = "purple", lty=2)  
####

####V2 vs W
regXY<-lm(W ~ V2, data = dat1)
summary(regXY)
newdat<-data.frame(V2 = seq(min(dat1$V2), max(dat1$V2), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V2,conf[,1],type = "l",col="red",
      xlab = TeX("V2"), ylab =TeX("W"), 
      main = TeX("V2 vs W"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V2,dat1$W,pch = 20,cex = 1.3)
matlines(newdat$V2, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V2, pred[,2:3], col = "purple", lty=2)  
####

####V3 vs W
regXY<-lm(W ~ V3, data = dat1)
summary(regXY)
newdat<-data.frame(V3 = seq(min(dat1$V3), max(dat1$V3), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V3,conf[,1],type = "l",col="red",
      xlab = TeX("V3"), ylab =TeX("W"), 
      main = TeX("V3 vs W"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V3,dat1$W,pch = 20,cex = 1.3)
matlines(newdat$V3, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V3, pred[,2:3], col = "purple", lty=2)  
####

####V4 vs W
regXY<-lm(W ~ V4, data = dat1)
summary(regXY)
newdat<-data.frame(V4 = seq(min(dat1$V1), max(dat1$V4), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V4,conf[,1],type = "l",col="red",
      xlab = TeX("V4"), ylab =TeX("W"), 
      main = TeX("V4 vs W"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V4,dat1$W,pch = 20,cex = 1.3)
matlines(newdat$V4, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V4, pred[,2:3], col = "purple", lty=2)  
####

####V5 vs W
regXY<-lm(W ~ V5, data = dat1)
summary(regXY)
newdat<-data.frame(V5 = seq(min(dat1$V5), max(dat1$V5), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V5,conf[,1],type = "l",col="red",
      xlab = TeX("V5"), ylab =TeX("W"), 
      main = TeX("V5 vs W"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V5,dat1$W,pch = 20,cex = 1.3)
matlines(newdat$V5, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V5, pred[,2:3], col = "purple", lty=2)  
####

par(mfrow=c(2,2))
####V1 vs V2
regXY<-lm(V1 ~ V2, data = dat1)
summary(regXY)
newdat<-data.frame(V2 = seq(min(dat1$V2), max(dat1$V2), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V2,conf[,1],type = "l",col="red",
      xlab = TeX("V2"), ylab =TeX("V1"), 
      main = TeX("V1 vs V2"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V2,dat1$V1,pch = 20,cex = 1.3)
matlines(newdat$V2, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V2, pred[,2:3], col = "purple", lty=2)  
####

####V1 vs V3
regXY<-lm(V1 ~ V3, data = dat1)
summary(regXY)
newdat<-data.frame(V3 = seq(min(dat1$V3), max(dat1$V3), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V3,conf[,1],type = "l",col="red",
      xlab = TeX("V3"), ylab =TeX("V1"), 
      main = TeX("V1 vs V3"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V3,dat1$V1,pch = 20,cex = 1.3)
matlines(newdat$V3, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V3, pred[,2:3], col = "purple", lty=2)  
####

####V1 vs V4
regXY<-lm(V1 ~ V4, data = dat1)
summary(regXY)
newdat<-data.frame(V4 = seq(min(dat1$V4), max(dat1$V4), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V4,conf[,1],type = "l",col="red",
      xlab = TeX("V4"), ylab =TeX("V1"), 
      main = TeX("V1 vs V4"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V4,dat1$V1,pch = 20,cex = 1.3)
matlines(newdat$V4, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V4, pred[,2:3], col = "purple", lty=2)    
####

####V1 vs V5
regXY<-lm(V1 ~ V5, data = dat1)
summary(regXY)
newdat<-data.frame(V5 = seq(min(dat1$V5), max(dat1$V5), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$V5,conf[,1],type = "l",col="red",
      xlab = TeX("V5"), ylab =TeX("V1"), 
      main = TeX("V1 vs V5"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V5,dat1$V1,pch = 20,cex = 1.3)
matlines(newdat$V5, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V5, pred[,2:3], col = "purple", lty=2)   
####

par(mfrow=c(1,3))
####V2 vs V3
regXY<-lm(V2 ~ V3, data = dat1)
summary(regXY)
newdat<-data.frame(V3 = seq(min(dat1$V3), max(dat1$V3), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V3,conf[,1],type = "l",col="red",
      xlab = TeX("V3"), ylab =TeX("V2"), 
      main = TeX("V2 vs V3"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V3,dat1$V2,pch = 20,cex = 1.3)
matlines(newdat$V3, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V3, pred[,2:3], col = "purple", lty=2)  
####

####V2 vs V4
regXY<-lm(V2 ~ V4, data = dat1)
summary(regXY)
newdat<-data.frame(V4 = seq(min(dat1$V4), max(dat1$V4), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V4,conf[,1],type = "l",col="red",
      xlab = TeX("V4"), ylab =TeX("V2"), 
      main = TeX("V2 vs V4"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V4,dat1$V2,pch = 20,cex = 1.3)
matlines(newdat$V4, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V4, pred[,2:3], col = "purple", lty=2)  
####

####V2 vs V5
regXY<-lm(V2 ~ V5, data = dat1)
summary(regXY)
newdat<-data.frame(V5 = seq(min(dat1$V5), max(dat1$V5), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V5,conf[,1],type = "l",col="red",
      xlab = TeX("V5"), ylab =TeX("V2"), 
      main = TeX("V2 vs V3"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V5,dat1$V2,pch = 20,cex = 1.3)
matlines(newdat$V5, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V5, pred[,2:3], col = "purple", lty=2)     
####

par(mfrow=c(2,1))
####V3 vs V4
regXY<-lm(V3 ~ V4, data = dat1)
summary(regXY)
newdat<-data.frame(V4 = seq(min(dat1$V4), max(dat1$V4), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V4,conf[,1],type = "l",col="red",
      xlab = TeX("V4"), ylab =TeX("V3"), 
      main = TeX("V3 vs V4"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V4,dat1$V3,pch = 20,cex = 1.3)
matlines(newdat$V4, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V4, pred[,2:3], col = "purple", lty=2)  
####

####V3 vs V5
regXY<-lm(V3 ~ V5, data = dat1)
summary(regXY)
newdat<-data.frame(V5 = seq(min(dat1$V5), max(dat1$V5), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V5,conf[,1],type = "l",col="red",
      xlab = TeX("V5"), ylab =TeX("V3"), 
      main = TeX("V3 vs V5"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V5,dat1$V3,pch = 20,cex = 1.3)
matlines(newdat$V5, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V5, pred[,2:3], col = "purple", lty=2)     
####

par(mfrow = c(1,1))
####V4 vs V5
regXY<-lm(V4 ~ V5, data = dat1)
summary(regXY)
newdat<-data.frame(V5 = seq(min(dat1$V5), max(dat1$V5), length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)

plotR(newdat$V5,conf[,1],type = "l",col="red",
      xlab = TeX("V5"), ylab =TeX("V4"), 
      main = TeX("V4 vs V5"), lwd = 2,ylim =c(min(pred[,2]),max(pred[,3])))
points(dat1$V5,dat1$V4,pch = 20,cex = 1.3)
matlines(newdat$V5, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$V5, pred[,2:3], col = "purple", lty=2)     
####

###Calculos multicolinealidad 
round(cor(dat1),3)


lm1<-lm(W~V1+V2+V3+V4+V5-1, data = dat1)
summary(lm1)
AIC(lm1)
#   Estimate Std. Error t value Pr(>|t|)  
#V1   -10175      60255  -0.169   0.8687  
#V2     4762       1733   2.748   0.0177 *
#V3    31159      58132   0.536   0.6018  
#V4    -1818       2968  -0.613   0.5516  
#V5    -2495       1272  -1.961   0.0735 .
#Rsquared 0.9908
#AIC: 272.6599


###Factores de inflación 
X1<-data.frame(dat1[,1:5])
names(X1)
names(X1)<-c("X1","X2","X3","X4","X5")
X1<-as.matrix(X1)

#VIFS
round(diag(solve(cor(X1))),3)
round(cor(X1),3)
heatmap(cor(X1))

###Valores y Vectores Propios 

(vals<-eigen(t(X1)%*%X1)$values) 

#Número de condición y matriz de descomposición de varianzas. 
round((max(vals)/vals),2)

#Número de condición elevado a la 1/2
round((max(vals)/vals)^(1/2),2)
(A<-VarDes(X1))

#`Valor Singular`  `var 1` `var 2`  `var 3` `var 4` `var 5` `Indice de Cond`
#     2.05    5.85e- 6 0.00616 6.27e- 6 0.00217 0.00620             1   
#     0.817   6.42e-10 0.0213  1.22e-10 0.00621 0.277               2.51
#     0.308   3.04e- 5 0.861   2.80e- 5 0.131   0.0328              6.66
#     0.202   5.61e- 4 0.109   7.13e- 4 0.424   0.485              10.2 
#    0.00735  9.99e- 1 0.00313 9.99e- 1 0.437   0.199             279. 
####Modelo Con Menos Regresoras

### V2 + V5
lm1<-lm(W~V2+V5-1, data = dat1)
summary(lm1)
AIC(lm1)
#R^2: 0.9239
#AIC:302.6191 

### V4 + V5
lm1<-lm(W~V4+V5-1, data = dat1)
summary(lm1)
AIC(lm1)
#R^2: 0.9104
#AIC: 305.3942


### V3 + V5
lm1<-lm(W~V3+V5-1, data = dat1)
summary(lm1)
AIC(lm1)
#AIC: 275.3023
#R^2: 0.9847
#   Estimate Std. Error t value Pr(>|t|)    
#V3  24189.2      956.8  25.282 1.03e-13 ***
#V5  -3362.5      956.8  -3.514  0.00313 ** 



##Eliminando V1, V2 y V4

X2<-as.matrix(dat1[,c(3,5)])


vals<-eigen(t(X2)%*%X2)$values
max(vals)/min(vals)
(max(vals)/min(vals))^(1/2)

diag(solve(cor(X2)))
cor(X2)
heatmap(cor(X2))

dat1<-read.table(file.choose(),header = TRUE, sep = ",")

X4<-dat1[,c(3,5)]
W<-dat1[,6]
beta<-lm1$coeff

means <- apply(X4,2,mean)
(Dest  <- apply(X4,2,function(x) x -mean(x)))
(Dest  <- apply(Dest,2,function(x) (sum(x^2))^(1/2)))

(beta1 <- beta*(1/Dest))
(beta1 <- c(mean(W) - sum(beta1*means),beta1))


### V3 + V5 comprobación.  
lm1<-lm(W~V3 + V5, data = dat1)
summary(lm1)