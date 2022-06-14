library(latex2exp)
library(statmod) ##qresid 
library(MASS)

rm(list = ls())
###Gráficas
plotR<-function(x,y,GRID = 1,...){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}


x<-1981:1993 
length(x)
###Codificando los años como 0 = 1981, 1 = 1982, ..., 11 = 1992, 12 = 1993
x<-0:12
y<-c(12,14,33,50,67,74,123,141,165,204,253,246,240) 

plotR(x,log(y), pch = 20, 
     main = "Años vs Casos de AIDS",
     ylab = "No. de Casos", xlab = "Codigo año")

#####Proceso Iterativo Con Intercepto
n<-length(y)
mu<-sapply(1:n,function(i)ifelse(y[i] > 0, y[i], y[i] + 0.1))
eta<-log(mu)
z<-log(mu) + (y - mu)/mu
w<-mu 
D0<-2*sum(y*log(y/mu) - (y-mu))
D<-1
i<-1

while(D >= 10^{-8}){
lmod<-lm(z~x,weights=w)
eta<-lmod$fit 
mu<-exp(eta)
w<-mu
z<-log(mu) + (y - mu)/mu
D1<- 2*sum(y*log(y/mu) - (y-mu))
cat(i,round(c(coef(lmod),D1),3),"\n")
i<-i+1
D<-abs(D0 - D1)/(abs(D1) + 0.1)
D0<-D1
}


#####Proceso Iterativo con sqrt 
n<-length(y)
mu<-sapply(1:n,function(i)ifelse(y[i] > 0, y[i], y[i] + 0.01))
eta<-log(mu)
z<-log(mu) + (y - mu)/mu
w<-mu 
D0<-2*sum(y*log(y/mu) - (y-mu))
D<-1
i<-1

while(D >= 10^{-8}){
lmod<-lm(z~x + sqrt(x),weights=w)
eta<-lmod$fit 
mu<-exp(eta)
w<-mu
z<-log(mu) + (y - mu)/mu
D1<- 2*sum(y*log(y/mu) - (y-mu))
cat(i,c(coef(lmod),D1),"\n")
i<-i+1
D<-abs(D0 - D1)/(D1 + 0.1)
D0<-D1
}

####GLM R, con intercepto
fit1<-glm(y~x, family = poisson)
summary(fit1)
####GLM R, con sqrt
fit21<-glm(y~ x + sqrt(x) , family = poisson)
summary(fit21)
anova(fit21,test= "Chisq")




neWx<-data.frame(x=seq(min(x),max(x),length.out= 20))

(preds <- predict(fit1, newdata = neWx, type = "link", se.fit = TRUE))

critval <- qnorm(0.975) ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- fit1$family$linkinv(fit)
upr2 <- fit1$family$linkinv(upr)
lwr2 <- fit1$family$linkinv(lwr)
int<-cbind(lwr2,upr2)
ylim<-c(min(lwr2),max(upr2))

plotR(neWx$x,fit2,type = "l",ylim = ylim)
matlines(neWx$x, int, col = "#0082FF", lty=2,type = "l")
points(x,y,pch = 20)

#####
(preds <- predict(fit21, newdata = neWx, type = "link", se.fit = TRUE))

critval <- qnorm(0.975) ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- fit1$family$linkinv(fit)
upr2 <- fit1$family$linkinv(upr)
lwr2 <- fit1$family$linkinv(lwr)
int<-cbind(lwr2,upr2)
ylim<-c(min(lwr2),max(upr2))


plotR(neWx$x,fit2,type = "l",ylim = ylim, 
      main = "Regresión Poisson Casos de AIDS",
      xlab = "No. Año", ylab = "Casos")
matlines(neWx$x, int, col = "#0082FF", lty=2,type = "l")
points(x,y, pch = 20)






