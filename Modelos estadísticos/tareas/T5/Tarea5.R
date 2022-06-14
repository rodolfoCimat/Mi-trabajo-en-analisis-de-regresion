#####################Paquetería
rm(list =ls())
library(ggplot2)   
library(tseries)      
library(GGally)
library(latex2exp)
library(plot3D)
library(dplyr)
########################Ejercicios###############################
dat1<-read.table(file.choose(),header = TRUE, sep = ",")
dat2<-read.table(file.choose(),header = TRUE, sep = ",")
dat4<-read.table(file.choose(),header = TRUE, sep = ",")

#######################Funciones auxiliares######################

DWest<-function(e)sum((e[2:length(e)]-e[1:(length(e)-1)])^2)/sum(e^2)

StaStu<-function(e,h,s,n,p){
  r <- e/(s*sqrt(1 - h))
  r2<- r*((n - p - 1)/(n - p -r^2))^(1/2)
  return(list(sta= r,stu= r2))  
}

EstErrs<-function(n,p,h,e,s)sqrt(((n-p)/(n-p-1))*s^2 - e^2/((n-p-1)*(1 - h)) )
D<-function(r,p,h)((r^2)/p)*(h/(1 - h))
DFF<-function(t,h)t*sqrt(h/(1 - h))


DFB<-function(X,t,h){
   C<-solve(t(X)%*%X)
   R<-C%*%t(X)
   A<-matrix( , nr = nrow(X),nc = ncol(X))
   sapply(1:ncol(X),function(j)sapply(1:nrow(X),function(i)A[i,j]<<- (R[j,i]/sqrt(C[j,j]))*(t[i]/sqrt(1 - h[i]))))  
   return(A) 
}
COVR<-function(t,h,n,p)((((n - p - 1)/(n - p) + (t^2)/(n - p))^(p))*(1 - h))^(-1)



#######################Ejercicio 1###############################
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

###### Ejercicio 1 
names(dat1)
int<-lm(Y ~ X1 + X2, data = dat1)
summary(int)
AIC(int)
sinint<-lm(Y ~ X1 + X2 - 1, data = dat1)
summary(sinint)
AIC(lm(Y ~ X1 + X2-1, data = dat1))
ggpairs(dat1[,2:4])


regXY<-lm(Y ~ X2, data = dat1)
summary(regXY)
newdat<-data.frame(X2 = seq(min(dat1$X2), max(dat1$X2)+5, length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$X2,conf[,1],type = "l",col="red",
      xlab = TeX("X_2"), ylab =TeX("Y"), 
      main = TeX("X_2 vs Y"), lwd = 2, ylim = c(min(pred[,2]),85))
points(dat1$X2,dat1$Y,pch = 20)
matlines(newdat$X2, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$X2, pred[,2:3], col = "purple", lty=2)  

e<-int$resid
n<-nrow(dat1)
p<-ncol(dat1)-1
h<-unname(lm.influence(int)$hat)

###1 


###2 Correlación simple entre cajas X1 y X2
 cor(dat1$X1,dat1$X2)
###3 

head(dat1)
XX<-as.matrix(dat1[,2:3])
RXX <- cor(XX)
round(RXX,3) 
round(solve(RXX),3)
diag(solve(RXX)) > 10


##4 Conditioning number  < 30 
X<-unname(as.matrix(cbind(rep(1,nrow(dat1)),dat1[,2:3])))

svd(t(X)%*%X)
(6.728369e+06 )/(8.818801e+00)

####Belsley y Kuh
X[,1] <- X[,1]/sqrt(sum(X[,1]^2))   
X[,2] <- X[,2]/sqrt(sum(X[,2]^2))
X[,3] <- X[,3]/sqrt(sum(X[,3]^2)) 

svd(X)
1.6367604/0.2566265 


##Sugerencia Profesor
X<-unname(as.matrix(cbind(dat1[,2:3])))
X[,1] <- (X[,1] -mean(X[,1]))/sqrt(sum(X[,1]^2))   
X[,2] <- (X[,2] -mean(X[,2]))/sqrt(sum(X[,2]^2))

eigen(t(X)%*%X)
eigen(t(X)%*%X)$values[1]/eigen(t(X)%*%X)$values[2]

X<-unname(as.matrix(cbind(dat1[,2:3])))
X[,1] <- (X[,1] )/sqrt(sum(X[,1]^2))   
X[,2] <- (X[,2] )/sqrt(sum(X[,2]^2))


############# Ejercicio 2
##1 
head(dat2)
regXY<-lm(y ~ X, data = dat2)
summary(regXY)
newdat<-data.frame(X = seq(min(dat2$X), max(dat2$X)+5, length = 500))
conf <- predict(regXY,interval = "confidence",newdata = newdat)
pred <- predict(regXY,interval = "prediction",newdata = newdat)
head(conf)
head(pred)
AIC(regXY)

plotR(newdat$X,conf[,1],type = "l",col="red",
      xlab = TeX("X"), ylab =TeX("y"), 
      main = TeX("X vs y"), lwd = 2, ylim = c(min(pred[,2]),max(pred[,3])))
points(dat2$X,dat2$y,pch = 20)
matlines(newdat$X, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$X, pred[,2:3], col = "purple", lty=2)  

e<-regXY$resid
(n<-nrow(dat2))
p<-ncol(dat2)-1
h<-unname(lm.influence(regXY)$hat)
##2 
plotR(e[1:(length(e)-1)],e[2:length(e)],
      xlab = TeX("e_{i-1}"),
      ylab = TeX("e_{i}"),
      main = TeX("e_{i} vs e_{i-1}"),
      pch = 20) 

##3 

DWest(e)
4 - DWest(e)
n 
p
### parametros 3 no de datos 25, 1% dL=1.20 dU=1.40, 5% dL = 1.39, dU = 1.60

acf(e,type = "partial",main = "PACF Residuales",lag.max =32
    ,ylab = TeX("$\\phi_{k,k}$"),xlab = "k")
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2)

##4 Cochrane & Orcutt

et<-e[2:length(e)]
et_1<-e[1:(length(e)-1)]
rho <- unname(coef(lm(et ~ et_1 - 1)))
round(rho,3)
Dat2<-data.frame(X = dat2$X[2:n] - rho*dat2$X[1:(n-1)],
                 y = dat2$y[2:n] - rho*dat2$y[1:(n-1)]) 

EEE<-rep(1-rho,n-1)
CORc<-lm(y ~ EEE+ X - 1, data = Dat2)
summary(CORc)
summary(regXY)
AIC(CORc)
AIC(regXY)
(n<-nrow(Dat2)) 
(p<-ncol(Dat2))

newdat<-data.frame(X = seq(min(Dat2$X), max(Dat2$X)+5, length = 500),
                   EEE= rep(1-rho),500)
conf <- predict(CORc,interval = "confidence",newdata = newdat)
pred <- predict(CORc,interval = "prediction",newdata = newdat)
head(conf)
head(pred)

plotR(newdat$X,conf[,1],type = "l",col="red",
      xlab = TeX("X_{i} - \\rho X_{i-1}"), ylab =TeX("y_{i} - \\rho y_{i-1}"), 
      main = TeX("Modelo corregido"), lwd = 2, ylim = c(min(pred[,2]),max(pred[,3])))
points(Dat2$X,Dat2$y,pch = 20)
matlines(newdat$X, conf[,2:3], col = "#0082FF", lty=2)
matlines(newdat$X, pred[,2:3], col = "purple", lty=2)  

(beta<-CORc$coeff)
newx<-seq(min(dat2$X),max(dat2$X) + 5,length = 500)
y0<-beta[1] + newx*beta[2]

plotR(newdat$X,conf[,1],type = "l",col="red",
      xlab = TeX("X"), ylab =TeX("y"), 
      main = TeX("X vs y"), lwd = 2, ylim = c(min(y0),max(pred[,3])))
lines(newx,y0,col= "blue",lwd = 2)


names(CORc)
e<-CORc$resid


plotR(e[1:(length(e)-1)],e[2:length(e)],
      xlab = TeX("f_{i-1}"),
      ylab = TeX("f_{i}"),
      main = TeX("f_{i} vs f_{i-1}"),
      pch = 20) 

################### 1.38 1.60 5%,  1.19 1.39 1%. 
round(DWest(e),3)
round(4 - DWest(e),3)
n 
p

acf(e,type = "partial",main = "PACF Residuales",lag.max =32
    ,ylab = TeX("$\\phi_{k,k}$"),xlab = "k")
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2)



###################Ejercicio 3




####function(e,h,s,n,p)
(s<-summary(int)$sigma)
(e<-int$resid)
(n<-nrow(dat1))
(p<-ncol(dat1)-1)
(h<-unname(lm.influence(int)$hat))
RES<-StaStu(e,h,s,n,p)
(r<-unname(RES[[1]]))
(t<-unname(RES[[2]]))
X<-unname(cbind(rep(1,length(n)),dat1[,2:3]))
(X<-as.matrix(X))

#####################1 Cuales son potencialmente 
which( h > 2*p/n)
2*p/n
plotR(1:n,h,col = "black",pch= 20, 
      main = "Diagonal Matriz de Proyección",
      xlab = "Obs",
      ylab = "h")
abline(h = 2*p/n, col = "red")
round(h[which( h > 2*p/n)],3)

#EstErrs<-function(n,p,h,e,s)sqrt(((n-p)/(n-p-1))*s^2 - e^2/((n-p-1)*(1 - h)) )
#D<-function(r,p,h)((r^2)/p)*(h/(1 - h))
#DFF<-function(t,h)t*(h/(1 - h))
#DFB<-function(X,t,h){
#   R<-solve(t(X)%*%X)*t(X)
#   A<-matrix( , nr = nrow(X),nc = ncol(X))
#   sapply(1:ncol(X),function(j)sapply(1:nrow(X),function(i)A[i,j]<<- (R[j,i]/sum(R[ ,j]^2))*(t[i]/sqrt(1 - h[i]))))  
#   return(A) 
#}
#COVR<-function(t,h,n,p)(((n - p -1)/(n - p) + (t^2)/(n - p))*(1 - h))^(-1)
2*sqrt(p/n)
s<-EstErrs(n,p,h,e,s)
(dis<-unname(DFF(t,h)*(s*sqrt(h)) ))

which( D(r,p,h) >  qf(0.5,p,n-p)) 
D(r,p,h)[which( D(r,p,h) >  qf(0.5,p,n-p))]
which(abs(DFF(t,h)) > 2*sqrt(p/n)) 
DFF(t,h)[which(abs(DFF(t,h)) > 2*sqrt(p/n))] 
which(abs(DFB(X,t,h)) > 2/sqrt(n),arr.ind = TRUE) 
DFB(X,t,h)[which(abs(DFB(X,t,h)) > 2/sqrt(n),arr.ind = TRUE)]
which( COVR(t,h,n,p) > 1 + 3*p/n)
COVR(t,h,n,p)[which( COVR(t,h,n,p) > 1 + 3*p/n)]
which( COVR(t,h,n,p) < 1 - 3*p/n) 
COVR(t,h,n,p)[which( COVR(t,h,n,p) < 1 - 3*p/n)]

plotR(1:n,abs(dis),col = "black",pch= 20, 
      main = TeX("Desplazamiento en Y_i"),
      xlab = "Obs",
      ylab = TeX("|Y_{i} - Y_{i(i)}|"))
points(c(9,22),abs(dis)[c(9,22)],col = "red",pch = 20)
abs(dis)[c(9,22)]

plotR(1:n,D(r,p,h),col = "black",pch= 20, 
      main = "Cook's D",
      xlab = "Obs",
      ylab = TeX("D_i"))
points(c(9,22),D(r,p,h)[c(9,22)],col = "red",pch = 20)
abline(h = qf(0.5,p,n-p),col = "red")

plotR(1:n,abs(DFF(t,h)),col = "black",pch= 20, 
      main = "DFFITS",
      xlab = "Obs",
      ylab = TeX("|DFFITS_i|"))
points(c(9,22),abs(DFF(t,h))[c(9,22)],col = "red",pch = 20)
abline(h = 2*sqrt(p/n),col = "red")


plotR(1:n,abs(DFB(X,t,h)[,1]),col = "black",pch= 20, 
      main = TeX("DFFBETAS_{(0)}"),
      xlab = "Obs",
      ylab = TeX("|DFFBETAS_{(0)i}|"))
points(c(9,22),abs(DFB(X,t,h))[c(9,22),1],col = "red",pch = 20)
abline(h = 2/sqrt(n),col = "red")

plotR(1:n,abs(DFB(X,t,h)[,2]),col = "black",pch= 20, 
      main = TeX("DFFBETAS_{(1)}"),
      xlab = "Obs",
      ylab = TeX("|DFFBETAS_{(1)i}|"))
points(c(9,22),abs(DFB(X,t,h))[c(9,22),2],col = "red",pch = 20)
abline(h = 2/sqrt(n),col = "red")

plotR(1:n,abs(DFB(X,t,h)[,3]),col = "black",pch= 20, 
      main = TeX("DFFBETAS_{(2)}"),
      xlab = "Obs",
      ylab = TeX("|DFFBETAS_{(2)i}|"))
points(c(9,22),abs(DFB(X,t,h))[c(9,22),3],col = "red",pch = 20)
abline(h = 2/sqrt(n),col = "red")

plotR(1:n,COVR(t,h,n,p),col = "black",pch= 20, 
      main = "COVRATIO",
      xlab = "Obs",
      ylab = TeX("COVRATIO_i"),
      ylim = c(0.5 - 3*p/n,1.5 + 3*p/n))
abline(h = 1 + 3*p/n,col = "red")
abline(h = 1 - 3*p/n,col = "red")






#influence.measures(int)

##Plot 3d 

X1 <- seq(min(dat1$X1),max(dat1$X1),length = n)
X2 <- seq(min(dat1$X2),max(dat1$X2),length = n)
newx<-expand.grid(X1 = X1, X2 = X2)
z<-matrix(predict(int,newdata = newx),nc=n,nr=n)
fited<-predict(int)
scatter3D(dat1$X1,dat1$X2,dat1$Y,bty = "b2",
         pch = 20,cex = 2, theta = 0,phi=0,
         main = "Observaciones Entrega de Refrescos",
         xlab = "Cajas de Producto Abastecido",
         #ylab = "Distancia Recorrida",
         zlab = "Tiempo de Entrega")
text3D(dat1$X1,dat1$X2,dat1$Y, labels = 1:n,
       add = TRUE, colkey = FALSE, cex = 0.8) 

scatter3D(dat1$X1,dat1$X2,dat1$Y,bty = "b2",
         pch = 20,cex = 2, theta = 20,phi=20,
         surf = list(x = X1, y = X2, z = z,  
         facets = NA, fit = fited),
         main = "Observaciones Entrega de Refrescos + Regresión",
         xlab = "Cajas de Producto Abastecido",
         ylab = "Distancia Recorrida",
         zlab = "Tiempo de Entrega")

scatter3D(dat1$X1,dat1$X2,dat1$Y,bty = "b2",
         pch = 20,cex = 2, theta = 20,phi=-5,
         surf = list(x = X1, y = X2, z = z,  
         facets = NA, fit = fited),
         main = "Observaciones Entrega de Refrescos + Regresión",
         xlab = "Cajas de Producto Abastecido",
         ylab = "Distancia Recorrida",
         zlab = "Tiempo de Entrega")


dat1 <- dat1[which(dat1$X2 != 1460),] 
int2 <- lm(Y ~ X1 + X2, data = dat1)
summary(int2)
AIC(int2)

###Estudentizados
alpha <- 0.05 
cV <- qt(1 - alpha/2,df = n-p-1)
which(abs(t) > cV)
abs(t)[which(abs(t) > cV)]
names(int)

plotR(int$fitted.values,t,ylim = c(-cV-2.5,cV + 2.5),
      pch = 20, main = "Residuales Studentizados vs Valores Ajustados",
      xlab = TeX("\\hat{Y}"),ylab = "t")
points(int$fitted.values[c(9,22)],t[c(9,22)],col = "red",pch = 20)
abline(h = cV, col ="red")
abline(h = -cV, col ="red")

##Ejercicio 4 
dat4
dat4<-as_tibble(dat4)
dat4 <-mutate(dat4, Localizacion = factor(Localizacion, 
           levels = c("Colonia","CentroCom.","Centro"),
           labels = 0:2))
dat40<-filter(dat4,Localizacion ==0) 
dat41<-filter(dat4,Localizacion ==1) 
dat42<-filter(dat4,Localizacion ==2) 
lm0<-lm(Ventas ~ NumCasas, data = dat40)
lm1<-lm(Ventas ~ NumCasas, data = dat41)
lm2<-lm(Ventas ~ NumCasas, data = dat42)

summary(lm0)
summary(lm1)
summary(lm2)
newx0 <- data.frame(NumCasas = seq(min(dat4$NumCasas),
                    max(dat4$NumCasas),length =500))  
newx1 <- data.frame(NumCasas = seq(min(dat4$NumCasas),
                    max(dat4$NumCasas),length =500)) 
newx2 <- data.frame(NumCasas = seq(min(dat4$NumCasas),
                    max(dat4$NumCasas),length =500)) 
y0<- predict(lm0, newdata = newx0)
y1<- predict(lm1, newdata = newx1)
y2<- predict(lm2, newdata = newx2)

plotR(dat4$NumCasas,dat4$Ventas,type = "n", 
      main = "Numero de Casas vs Ventas", 
      xlab = "Número de Casas",
      ylab = "Ventas")
 
points(dat40$NumCasas,dat40$Ventas,pch = 20, col ="red") 
points(dat41$NumCasas,dat41$Ventas,pch = 15, col = "blue") 
points(dat42$NumCasas,dat42$Ventas,pch = 17, col = "black") 
lines(newx0$NumCasas,y0, col = "red")
lines(newx1$NumCasas,y1, col = "blue")
lines(newx2$NumCasas,y2, col = "black")

##Gráficas individuales

#1
newdat<-data.frame(NumCasas = seq(min(dat40$NumCasas),
                   max(dat40$NumCasas),length =500)) 
conf <- predict(lm0,interval = "confidence",newdata = newdat)
pred <- predict(lm0,interval = "prediction",newdata = newdat)


plotR(newdat$NumCasas,conf[,1],type = "l",col="red",
      main = TeX("Numero de Casas vs Ventas, Colonia."), ylab =TeX("Ventas"), 
      xlab = TeX("Numero de Casas"), lwd = 2, ylim = c(min(pred[,2]),max(pred[,3])))
points(dat40$NumCasas,dat40$Ventas,pch = 20)
matlines(newdat$NumCasas, conf[,2:3], col = "blue", lty=2)
matlines(newdat$NumCasas, pred[,2:3], col = "purple", lty=2)

#2
newdat<-data.frame(NumCasas = seq(min(dat41$NumCasas),
                   max(dat41$NumCasas),length =500)) 
conf <- predict(lm1,interval = "confidence",newdata = newdat)
pred <- predict(lm1,interval = "prediction",newdata = newdat)


plotR(newdat$NumCasas,conf[,1],type = "l",col="red",
      main = TeX("Numero de Casas vs Ventas, CentroCom."), ylab =TeX("Ventas"), 
      xlab = TeX("Numero de Casas"), lwd = 2, ylim = c(min(pred[,2]),max(pred[,3])))
points(dat41$NumCasas,dat41$Ventas,pch = 20)
matlines(newdat$NumCasas, conf[,2:3], col = "blue", lty=2)
matlines(newdat$NumCasas, pred[,2:3], col = "purple", lty=2)

#3
newdat<-data.frame(NumCasas = seq(min(dat42$NumCasas),
                   max(dat42$NumCasas),length =500)) 

conf <- predict(lm2,interval = "confidence",newdata = newdat)
pred <- predict(lm2,interval = "prediction",newdata = newdat)


plotR(newdat$NumCasas,conf[,1],type = "l",col="red",
      main = TeX("Numero de Casas vs Ventas, Centro."), ylab =TeX("Ventas"), 
      xlab = TeX("Numero de Casas"), lwd = 2, ylim = c(min(pred[,2])
      ,max(pred[,3])))
points(dat42$NumCasas,dat42$Ventas,pch = 20)
matlines(newdat$NumCasas, conf[,2:3], col = "blue", lty=2)
matlines(newdat$NumCasas, pred[,2:3], col = "purple", lty=2) 

dat44<-mutate(dat4,dum1 = 1*(Localizacion == 0),dum2 = 1*(Localizacion == 1))
dat44<-select(dat44,-c("Tienda","Localizacion"))
lm4<-lm(Ventas ~ NumCasas + dum1 + dum2, data = dat44)
summary(lm4)
AIC(lm4)

##Dif estimadas
a<-unname(lm4$coeff)
 round(-a[3],3)
 round(-a[4],3)
 round(a[3]-a[4],3)
### Errores estandar

(R<-unname(vcov(lm4)))

a1<-matrix(c(0,0,-1,0),nc=4,nr=1)
a2<-matrix(c(0,0,0,-1),nc=4,nr=1)
a3<-matrix(c(0,0,1,-1),nc=4,nr=1)
##Estimación para el error estandar de la dif 1 en la tabla
ee1<-sqrt(a1%*%R%*%t(a1))
(round(sqrt(a1%*%R%*%t(a1)),3))
##Estimación para el error estandar de la dif 2 en la tabla
ee2<-sqrt(a2%*%R%*%t(a2))
(round(sqrt(a2%*%R%*%t(a2)),3))
##Estimación para el error estandar de la dif 3 en la tabal
ee3<-sqrt(a3%*%R%*%t(a3))
round(sqrt(a3%*%R%*%t(a3)),3)
###t values
ee1
a[1]
(t1<--a[3]/ee1)
(t2<--a[4]/ee2)
(t3<-(a[3]-a[4])/ee3)
### p values
n<-15
p<-4

1-(pt(abs(t1),n-p) - pt(-abs(t1),n-p))
1-(pt(abs(t2),n-p) - pt(-abs(t2),n-p))
1-(pt(abs(t3),n-p) - pt(-abs(t3),n-p))




dat45<-mutate(dat4,dum2 = 1*(Localizacion == 1))
dat45<-select(dat45,-c("Tienda","Localizacion"))
lm5<-lm(Ventas ~ NumCasas + dum2, data = dat45)
summary(lm5)
AIC(lm5)

beta <-lm5$coeff
a1<-beta[1]
a2<-beta[1] + beta[3]
b<- beta[2]

newx0 <- seq(min(dat4$NumCasas),max(dat4$NumCasas),length =500)  

y0<- a1 + b*newx0
y1<- a2 + b*newx0


plotR(dat4$NumCasas,dat4$Ventas,type = "n", 
      main = "Numero de Casas vs Ventas", 
      xlab = "Número de Casas",
      ylab = "Ventas", ylim = c(min(y0),max(y1)))
 
points(dat40$NumCasas,dat40$Ventas,pch = 17, col ="black") 
points(dat41$NumCasas,dat41$Ventas,pch = 20, col = "blue") 
points(dat42$NumCasas,dat42$Ventas,pch = 17, col = "black") 
lines(newx0,y0, col = "black")
lines(newx0,y1, col = "blue")

##Diferencias en pendientes
dat46<-mutate(dat4,dum1 = 1*(Localizacion == 1),zx1 = NumCasas*(Localizacion ==0),zx2 = NumCasas*(Localizacion ==1) )
dat46 <- select(dat46,-c("Tienda","Localizacion"))
dat46
lm6<-lm(Ventas ~ NumCasas + dum1 + zx1 + zx2, data = dat46)
summary(lm6)
AIC(lm6)

(beta <-lm6$coeff)
a1<-beta[1]
a2<-beta[1] + beta[3]
b1<- beta[2]+ beta[4]
b2<- beta[2]+ beta[5] 
b3<- beta[2] 

newx0 <- seq(min(dat4$NumCasas),max(dat4$NumCasas),length =500)  

y0<- a1 + b1*newx0
y1<- a2 + b2*newx0
y2<- a1 + b3*newx0

points(dat40$NumCasas,dat40$Ventas,pch = 17, col ="black") 
points(dat41$NumCasas,dat41$Ventas,pch = 20, col = "blue") 
points(dat42$NumCasas,dat42$Ventas,pch = 15, col = "red") 
lines(newx0,y0, col = "black")
lines(newx0,y1, col = "blue")
lines(newx0,y2, col = "red")

#######
dat47<-mutate(dat4,dum1 =1*(Localizacion == 0),dum2 = 1*(Localizacion == 1),zx1 = NumCasas*(Localizacion ==0),zx2 = NumCasas*(Localizacion ==1) )
dat47 <- select(dat47,-c("Tienda","Localizacion"))
dat47
lm7<-lm(Ventas ~ NumCasas + dum1 + dum2 + zx1 + zx2, data = dat47)
summary(lm7)
AIC(lm7)

(beta <-lm7$coeff)
a1<-beta[1] + beta[3]
a2<-beta[1] + beta[4]
a3<-beta[1]
b1<- beta[2]+beta[5]
b2<- beta[2]+beta[6] 
b3<- beta[2] 

newx0 <- seq(min(dat4$NumCasas),max(dat4$NumCasas),length =500)  

y0<- a1 + b1*newx0
y1<- a2 + b2*newx0
y2<- a3 + b3*newx0

plotR(dat4$NumCasas,dat4$Ventas,type = "n", 
      main = "Numero de Casas vs Ventas", 
      xlab = "Número de Casas",
      ylab = "Ventas")
 
points(dat40$NumCasas,dat40$Ventas,pch = 20, col ="red") 
points(dat41$NumCasas,dat41$Ventas,pch = 15, col = "blue") 
points(dat42$NumCasas,dat42$Ventas,pch = 17, col = "black") 
lines(newx0,y0, col = "red")
lines(newx0,y1, col = "blue")
lines(newx0,y2, col = "black")


####Diferencias 
(R<-unname(vcov(lm7)))
(a <-unname(lm7$coeff))

round(a[5]-a[6],3)
a1<-matrix(c(0,0,0,0,-1,0),nc=6,nr=1)
a2<-matrix(c(0,0,0,0,0,-1),nc=6,nr=1)
a3<-matrix(c(0,0,0,0,1,-1),nc=6,nr=1)
##Estimación para el error estandar de la dif 1 en la tabla
(ee1<-sqrt(a1%*%R%*%t(a1)))
round(ee1,3)
##Estimación para el error estandar de la dif 2 en la tabla
(ee2<-sqrt(a2%*%R%*%t(a2)))
round(ee2,3)
##Estimación para el error estandar de la dif 3 en la tabal
(ee3<-sqrt(a3%*%R%*%t(a3)))
round(ee2,3)
###t values

(t1<--a[5]/ee1)
(t2<--a[6]/ee2)
(t3<-(a[5]-a[6])/ee3)
### p values
n<-15
p<-6

1-(pt(abs(t1),n-p) - pt(-abs(t1),n-p))
1-(pt(abs(t2),n-p) - pt(-abs(t2),n-p))
1-(pt(abs(t3),n-p) - pt(-abs(t3),n-p))

s<-summary(lm7)$sigma 
X
##P-valor final
K <-cbind(c(0,0,1,0,0,0),c(0,0,0,0,1,0))
(t(K)%*%a)%*%solve()
pf( )

