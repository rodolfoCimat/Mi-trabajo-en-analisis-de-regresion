rm(list = ls())

library(statmod) ##qresid 
library(MASS)
library(latex2exp)


plotR<-function(x,y,GRID = 1,...){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}

polip1 <-c(1,1,2,3,3,4,17,25,33)
x1<-c(22,23,16,23,23,42,22,17,23)
polip2<-c(7,10,15,28,28,40,44,46,50,61,63)
x2<-c(34,30,50,18,22,27,19,22,34,13,20)
x<-c(x1,x2)
polip<-c(polip1,polip2)
un<-c(rep(0,length(x1)),rep(1,length(x2))) 


par(mfrow = c(1,1))
boxplot(polip ~ as.factor(un), xlab = "No. Grupo", 
     ylab = "No. de polipos en los individuos por grupo",
     main = "Boxplot de los datos", col = "#229954")

##Datos problema
Dat<- data.frame(pol = polip,edad = x,grup= un )
Dat

Dat0<-Dat[Dat$grup == 0, ] 
Dat1<-Dat[Dat$grup == 1, ]
(edades0<- unique(Dat0$edad))
(edades1<- unique(Dat1$edad)) 

A0 <- sapply(edades0,function(x)c(mean(Dat0[Dat0$edad ==x, 1]),x)) 
A1 <- sapply(edades1,function(x)c(mean(Dat1[Dat1$edad ==x, 1]),x))
ylim <-c(min(A0[1,]), max(A1[1,]))

plotR(A1[2,],A1[1,], pch = 20,col = "red"
      ,ylim = ylim, main = "Polipos Promedio por Edad por Grupo"
      ,ylab = "Polipos Promedio", xlab = "Edad")
points(A0[2,],A0[1,], pch = 18, col = "blue")

mean(polip1);var(polip1)
mean(polip2);var(polip2)


###Ajuste Poisson. 

fit<-glm(pol~.,family ="poisson",data = Dat )
summary(fit)
anova(fit, test = "Chisq")
deviance(fit);df.residual(fit)

##Modelo con más parámetros posible;
fitt<-glm(pol~edad*grup,family ="poisson",data = Dat )
summary(fitt)
deviance(fitt);df.residual(fitt) 

qres <- qresid(fitt); qqnorm(qres, las=1, main ="Normal Q-Q Plot: Standardized Residuals"); abline(0,1)
scatter.smooth(qres~fitted(fit), las=1, main="Residuals vs fitted",
               xlab="Fitted value", ylab="Standardized residual")
fittet<-fitted(fit)

im <- influence.measures(fit); names(im)
im$infmat <- round(im$infmat, 3 ); head( im$infmat )
colSums( im$is.inf )
### c(8,9,10,18)

Ind<-c(8,9,12,19,10,18)
#boxplot(polip ~ as.factor(un), xlab = "No. Grupo", 
#     ylab = "No. de polipos en los individuos por grupo",
#     main = "Boxplot de los datos", col = "#229954")

points(fittet[Ind],qres[Ind], pch = 20, col = "red")



###QP
fitQ<-glm(pol~.,family ="quasipoisson",data = Dat )
summary(fitQ)
summary(fitQ)$dispersion
sqrt(summary(fitQ)$dispersion)
d<-cooks.distance(fitQ)
im <- influence.measures(fitQ); names(im)
im$infmat <- round(im$infmat, 3 ); head( im$infmat )
colSums( im$is.inf )


par(mfrow = c(1,3))
qres <- rstandard(fitQ); qqnorm(qres, las=1, main ="Normal Q-Q Plot: Standardized Residuals"); abline(0,1)
scatter.smooth(qres~fitted(fitQ), las=1, main="Residuals vs fitted",
               xlab="Fitted value", ylab="Standardized residuals")
plot(d,type = "h",main = "Cook's Distance", ylab="Cook's distance", las=1) 

##NB
fitNB<-glm.nb(pol~.,data = Dat)
fitNB$theta
fitNB <- glm.convert(fitNB)
printCoefmat(coef(summary(fitNB,dispersion=1)))
deviance(fitNB);df.residual(fitNB)
d<-cooks.distance(fitNB)
im <- influence.measures(fitNB); names(im)
im$infmat <- round(im$infmat, 3 ); head( im$infmat )
colSums( im$is.inf )

par(mfrow = c(1,3))
(qres <- qresid(fitNB,dispersion = 1))
qqnorm(qres,las=1,ylim = c(min(qres),max(qres)),main ="Normal Q-Q Plot: Quantile Residuals" ); abline(0,1)
scatter.smooth(qres~fitted(fitNB), las=1, main="Residuals vs fitted",
               xlab="Fitted value", ylab="Quantile residual")
plot(d,type = "h",main = "Cook's Distance", ylab="Cook's distance", las=1) 
qf(0.5, 3 , 17) 


coef(fit)
se.m1 <- coef(summary(fit))[, "Std. Error"]
se.nb <- coef(summary(fitNB,dispersion=1))[, "Std. Error"]
se.qp <- coef(summary(fitQ))[, "Std. Error"]

round(data.frame(SE.Pois=se.m1, SE.Quasi=se.qp, ratio=se.qp/se.m1,SE.NB=se.nb),3)

coef.mat <- rbind( coef(fit), coef(fitQ), coef(fitNB) )
rownames(coef.mat) <- c("Poisson glm", "Quasi-Poisson", "Neg bin glm")
coef.mat

#######Gráfica del Modelo CuasiPoisson Con intervalos de Confianza
Dat0<-Dat0[ ,1:2]
Dat1<-Dat1[ ,1:2]

neWx0<-data.frame(edad=seq(10,max(x),length.out= 1000), grup= rep(0,1000))
neWx1<-data.frame(edad=seq(10,max(x),length.out= 1000), grup= rep(1,1000))

preds0 <- predict(fitQ, newdata = neWx0, type = "link", se.fit = TRUE)
preds1 <- predict(fitQ, newdata = neWx1, type = "link", se.fit = TRUE)

critval <- qt(0.975,df = 17) ## approx 95% CI
upr0 <- preds0$fit + (critval * preds0$se.fit)
lwr0 <- preds0$fit - (critval * preds0$se.fit)
fit0 <- preds0$fit
upr1 <- preds1$fit + (critval * preds1$se.fit)
lwr1 <- preds1$fit - (critval * preds1$se.fit)
fit1 <- preds1$fit

fit02 <- fitQ$family$linkinv(fit0)
upr02 <- fitQ$family$linkinv(upr0)
lwr02 <- fitQ$family$linkinv(lwr0)
int0<-cbind(lwr02,upr02)
ylim0<-c(min(lwr02),max(upr02))
fit12 <- fitQ$family$linkinv(fit1)
upr12 <- fitQ$family$linkinv(upr1)
lwr12 <- fitQ$family$linkinv(lwr1)
int1<-cbind(lwr12,upr12)
ylim1<-c(min(lwr12),max(upr12))
ylim<-c( min( ylim0,ylim1), max(ylim0,ylim1))

min(x)
par(mfrow = c(1,2))
plotR(neWx0$edad,fit02,type = "l",ylim = ylim, 
       xlab = "Edad", ylab = "No. Polipos",
       main= "Ajuste Quasi-poisson", lwd = 2, col = "blue")
matlines(neWx0$edad, int0, col = "#0082FF", lty=2,type = "l")
points(Dat0$edad,Dat0$pol,pch = 20,col="blue")

lines(neWx1$edad,fit12,type = "l", lwd = 2, col = "red")
matlines(neWx1$edad, int1, col = "#C0392B", lty=2,type = "l")
points(Dat1$edad,Dat1$pol,pch = 20,col="red")

############Gráfica del Modelo binomial Negativo 
preds0 <- predict(fitNB, newdata = neWx0, type = "link"
                  ,se.fit = TRUE,dispersion = 1)
preds1 <- predict(fitNB, newdata = neWx1, type = "link"
                  ,se.fit = TRUE,dispersion = 1)

critval <- qnorm(0.975) ## approx 95% CI
upr0 <- preds0$fit + (critval * preds0$se.fit)
lwr0 <- preds0$fit - (critval * preds0$se.fit)
fit0 <- preds0$fit
upr1 <- preds1$fit + (critval * preds1$se.fit)
lwr1 <- preds1$fit - (critval * preds1$se.fit)
fit1 <- preds1$fit

fit02 <- fitNB$family$linkinv(fit0)
upr02 <- fitNB$family$linkinv(upr0)
lwr02 <- fitNB$family$linkinv(lwr0)
int0<-cbind(lwr02,upr02)
ylim0<-c(min(lwr02),max(upr02))
fit12 <- fitNB$family$linkinv(fit1)
upr12 <- fitNB$family$linkinv(upr1)
lwr12 <- fitNB$family$linkinv(lwr1)
int1<-cbind(lwr12,upr12)
ylim1<-c(min(lwr12),max(upr12))
ylim<-c( min( ylim0,ylim1), max(ylim0,ylim1))

plotR(neWx0$edad,fit02,type = "l",ylim = ylim, 
       xlab = "Edad",ylab = "",
       main= "Ajuste Binomial Negativo", lwd = 2, col = "blue")
matlines(neWx0$edad, int0, col = "#0082FF", lty=2,type = "l")
points(Dat0$edad,Dat0$pol,pch = 20,col="blue")

lines(neWx1$edad,fit12,type = "l", lwd = 2, col = "red")
matlines(neWx1$edad, int1, col = "#C0392B", lty=2,type = "l")
points(Dat1$edad,Dat1$pol,pch = 20,col="red")

###########################################################################




##NB Sin la edad
fitNB2<-glm.nb(pol~grup,data = Dat)
fitNB2$theta
fitNB2 <- glm.convert(fitNB2)
printCoefmat(coef(summary(fitNB2,dispersion=1)))
AIC(fitNB2);AIC(fitNB)
deviance(fitNB2);df.residual(fitNB2)
d<-cooks.distance(fitNB)


printCoefmat(coef(summary(fitNB,dispersion=1)))
printCoefmat(coef(summary(fitNB2,dispersion=1)))

summary(fitNB,dispersion=1)
summary(fitNB2,dispersion=1)

predict(fitNB2, type = "response", dispersion = 1,
        se.fit = TRUE, newdata = data.frame(grup =c(0,1)))


