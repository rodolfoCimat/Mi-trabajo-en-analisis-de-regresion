library(statmod) ## qresid ##glm.scoretest 
library(MASS) ##glm.nb , glm.convert
library(latex2exp)


#####Gráficas chidas 

plotR<-function(x,y,GRID = 1,...){
plot(x,y,...)
if(GRID == 1){
grid()
}
axis(side = 1,lwd = 2)
axis(side = 2,lwd = 2)
box(lwd=2) 
}

####Grafica QQ residuales de cuantiles
plotQQqres<-function(fit,...){ 
           qres <- qresid(fit,...) 
           qqnorm(qres, las=1, main ="Normal Q-Q Plot: Standardized Residuals")
           abline(0,1)
           }
####Gráfica QQ residuales de la deviance standarizados
plotQQqsd<-function(fit,...){ 
           qstandard <- rstandard(fit,...) 
           qqnorm(qres, las=1, main ="Normal Q-Q Plot: Standardized Residuals")
           abline(0,1)
           }
###Conteo observaciones influyentes
InfluenDiag<-function(fit,...){
          im <- influence.measures(fit); names(im)
          im$infmat <- round(im$infmat, 3 ); head( im$infmat )
          print(colSums(im$is.inf))
          }
###Intervalos de confianza pa la respuesta caso phi conocida

cofInt<-function(fit,alpha = 0.05,...){
        preds <- predict(fit, type = "link", se.fit = TRUE,...)

        critval <- qnorm(1 - alpha/2) ## approx 95% CI
        upr <- preds$fit + (critval * preds$se.fit)
        lwr <- preds$fit - (critval * preds$se.fit)
        fitt <- preds$fit

        fit2 <- fit$family$linkinv(fitt)
        upr2 <- fit$family$linkinv(upr)
        lwr2 <- fit$family$linkinv(lwr)
        int<-cbind(lwr2,upr2)
        ylim<-c(min(lwr2),max(upr2))
        return(list(fit = fit2,int = int, ylim = ylim) ) 
       }
###Intervalos de confianza pa la respuesta caso phi desconocida
cofIntT<-function(fit,df,alpha = 0.05,...){
        preds <- predict(fit, type = "link", se.fit = TRUE,...)

        critval <- qt(1 - alpha/2,df=df) ## approx 95% CI
        upr <- preds$fit + (critval * preds$se.fit)
        lwr <- preds$fit - (critval * preds$se.fit)
        fitt <- preds$fit

        fit2 <- fit$family$linkinv(fitt)
        upr2 <- fit$family$linkinv(upr)
        lwr2 <- fit$family$linkinv(lwr)
        int<-cbind(lwr2,upr2)
        ylim<-c(min(lwr2),max(upr2))
        return(list(fit = fit2,int = int, ylim = ylim) ) 
       }

#### 


leu<-c(13,5,5,3,4,18)
total<-c(391,205,156,50,35,51)
(prop <- leu/total)
dosis<-0:5
A<-data.frame(leu = leu,total = total, prop = prop,
           d1 =dosis,d2=dosis^2,d3 =dosis^3)

 
fit1<-glm(prop ~d1,family = binomial,weights = total,data = A)
AIC(fit1);deviance(fit1);df.residual(fit1)
summary(fit1)
fit2<-glm(prop ~d1 + d2,family = binomial,weights = total,data = A)
AIC(fit2);deviance(fit2);df.residual(fit2)
summary(fit2)
fit3<-glm(prop ~d1 + d2 + d3,family = binomial,weights = total,data = A)
AIC(fit3);deviance(fit3);df.residual(fit3)
summary(fit3)


anova(fit3,test = "Chisq")


fit4<-glm(prop ~d2,family = binomial,weights = total,data = A)
summary(fit4)
AIC(fit4);deviance(fit4);df.residual(fit4)


newdata1<-data.frame(d1= seq(min(dosis), max(dosis), length =100))
newdata2<-data.frame(d1= seq(min(dosis), max(dosis), length =100),
                     d2= seq(min(dosis), max(dosis), length =100)^2)
newdata3<-data.frame(d1= seq(min(dosis), max(dosis), length =100),
                     d2= seq(min(dosis), max(dosis), length =100)^2,
                     d3= seq(min(dosis), max(dosis), length =100)^3)
newdata4<-data.frame(d2= seq(min(dosis), max(dosis), length =100)^2)
ci1 <- cofInt(fit1,newdata1,alpha = 0.05)
ci2 <- cofInt(fit2,newdata2,alpha = 0.05)
ci3 <- cofInt(fit3,newdata3,alpha = 0.05)
ci4 <- cofInt(fit4,newdata4,alpha = 0.05)

names(ci1)


par(mfrow = c(1,3))
plotR(newdata1$d1,ci1$fit,type = "l",ylim = ci1$ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad Mod1", lwd = 2)
matlines(newdata1$d1,ci1$int, col = "#0082FF", lty=2,type = "l")
points(dosis,prop,pch = 20,col="red",cex = 1.5)

plotR(newdata2$d1,ci2$fit,type = "l",ylim = ci2$ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad Mod2", lwd = 2)
matlines(newdata2$d1,ci2$int, col = "#0082FF", lty=2,type = "l")
points(dosis,prop,pch = 20,col="red",cex = 1.5)

plotR(newdata3$d1,ci3$fit,type = "l",ylim = ci3$ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad Mod3", lwd = 2)
matlines(newdata3$d1,ci3$int, col = "#0082FF", lty=2,type = "l")
points(dosis,prop,pch = 20,col="red",cex = 1.5)

par(mfrow = c(1,1))
plotR(newdata3$d1,ci4$fit,type = "l",ylim = ci4$ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad Mod 4", lwd = 2)
matlines(newdata3$d1,ci4$int, col = "#0082FF", lty=2,type = "l")
points(dosis,prop,pch = 20,col="red",cex = 1.5)


fit5<-glm(prop ~d3,family = binomial,weights = total,data = A)
summary(fit5)

newdata5<-data.frame(d3= seq(min(dosis), max(dosis), length =100)^3)
ci5 <- cofInt(fit5,newdata5,alpha = 0.05)

plotR(newdata3$d1,ci5$fit,type = "l",ylim = ci5$ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad Mod 4", lwd = 2)
matlines(newdata3$d1,ci5$int, col = "#0082FF", lty=2,type = "l")
points(dosis,prop,pch = 20,col="red",cex = 1.5)
