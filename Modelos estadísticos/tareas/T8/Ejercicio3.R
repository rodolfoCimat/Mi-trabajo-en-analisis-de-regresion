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


leu<-c(13,5,5,3,4,18)
total<-c(391,205,156,50,35,51)
(prop <- leu/total)
options(contrasts= c("contr.treatment", "contr.treatment"))
dosrad<-c("0","1-9","10-49","50-99","100-199","200+")
dosrad<-ordered(dosrad
                ,levels=c("0","1-9","10-49","50-99","100-199","200+"))

Data <- data.frame(prop = prop, w = total,dosrad = dosrad)
head(Data)

matplot( xtabs( prop ~ dosrad, data=Data), pch=20, lty=5,
        type="b", lwd=2, col="black", axes=FALSE,
        xlab="Dosis de Radiación", ylab="Proporción",
        main = "Proporción de muertes por Leucemia por Grupo")
axis(side=1, at=1:6, labels=levels(Data$dosrad))
axis(side=2, las=1); box()

#matplot(xtabs( log(prop/(1-prop)) ~ dosrad, data=Data), pch=20, lty=5,
#        type="b", lwd=2, col="black", axes=FALSE,
#        xlab="Dosis de Radiación", ylab="Proporción",
#        main = "Logit de la Proporción de muertes por Leucemia por Grupo")
#axis(side=1, at=1:6, labels=levels(Data$dosrad))
#axis(side=2, las=1); box()

#doses<-c(0,mean(c(1,9)), mean(c(10,49)), mean(c(50,99)), mean(c(100,199)),200 )
doses<-c(0,1,10,50,100,200)
Data$doses <-doses

###LOGISTIC
fitLOG<-glm(prop ~ doses, family=binomial,weights=w, data=Data)
summary(fitLOG)
AIC(fitLOG);df.residual(fitLOG);deviance(fitLOG)

###PROBIT
fitPROB<-update(fitLOG, family=binomial(link="probit") )
summary(fitPROB)
AIC(fitPROB);df.residual(fitPROB);deviance(fitPROB)

###CLOGLOG
fitCLL<-update(fitLOG,family=binomial(link="cloglog") )
summary(fitCLL)
AIC(fitCLL);df.residual(fitCLL);deviance(fitCLL)

tr.array <- rbind( coef(fitLOG), coef(fitPROB), coef(fitCLL))
tr.array <- cbind(tr.array, 
                  c(deviance(fitLOG),deviance(fitPROB), deviance(fitCLL)),
                  c(df.residual(fitLOG),df.residual(fitPROB), df.residual(fitCLL)),
                  c(AIC(fitLOG),AIC(fitPROB),AIC(fitCLL)))
rownames(tr.array) <- c("Logit","Probit","Comp log-log")
round(tr.array,4)

neWx<-data.frame(doses=seq(0,max(doses),length.out= 1000))

preds <- predict(fitLOG, newdata = neWx, type = "link", se.fit = TRUE)

critval <- qnorm(0.975) ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- fitLOG$family$linkinv(fit)
upr2 <- fitLOG$family$linkinv(upr)
lwr2 <- fitLOG$family$linkinv(lwr)
int<-cbind(lwr2,upr2)
ylim<-c(min(lwr2),max(upr2))

head(Data)

plotR(neWx$doses,fit2,type = "l",ylim = ylim, 
       xlab = "Dosis de Radiación", ylab = TeX("$\\pi$"),
       main= "Probabilidad", lwd = 2)
matlines(neWx$doses, int, col = "#0082FF", lty=2,type = "l")
points(doses,prop,pch = 20,col="red")


###Modelo Poisson para proporciones 

Data1 <- data.frame(leu = leu, Total = total,dosrad = dosrad)
Data1$doses<-Data$doses
head(Data1)
#matplot(xtabs( log(prop) ~ dosrad, data=Data), pch=20, lty=5,
#        type="b", lwd=2, col="black", axes=FALSE,
#        xlab="Dosis de Radiación", ylab="Proporción",
#        main = "Ln de la Proporción de muertes por Leucemia por Grupo")
#axis(side=1, at=1:6, labels=levels(Data$dosrad))
#axis(side=2, las=1); box()


Data1$doses<-Data$doses
dlc.m1 <- glm( leu ~ offset(log(total)) + doses,
               family=poisson, data=Data1)
summary(dlc.m1)

neWx<-data.frame(doses=seq(0,max(doses),length.out= 1000))

##ConfBands 95%
(preds <- predict(dlc.m1, newdata = neWx, type = "link", se.fit = TRUE))

critval <- qnorm(0.975) ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- dlc.m1$family$linkinv(fit)
upr2 <- dlc.m1$family$linkinv(upr)
lwr2 <- dlc.m1$family$linkinv(lwr)
int<-cbind(lwr2,upr2)
ylim<-c(min(lwr2),max(upr2))

plotR(neWx$doses,fit2,type = "l",ylim = ylim)
matlines(neWx$doses, int, col = "#0082FF", lty=2,type = "l")
points( , pch = 20)




