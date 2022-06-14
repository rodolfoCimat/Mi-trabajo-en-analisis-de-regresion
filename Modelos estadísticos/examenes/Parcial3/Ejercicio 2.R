Y<-c(2,3,6,7,8,9,10,12,15)
X<-c(-1,-1,0,0,0,0,1,1,1)
B<-data.frame(X = X, Y = Y)

fit1 <- glm(Y ~ X,family = "poisson",data = B)
summary(fit1)
fit2 <- glm(Y ~ X,family = poisson("identity"),data = B)
summary(fit2)

newdata1<-data.frame(X= seq(-1, 1, length = 100))

(pred1 <- cofInt(fit1,alpha = 0.05))
(pred2 <- cofInt(fit2,alpha = 0.05))

ci1 <- cofInt(fit1,newdata1,alpha = 0.05)
ci2 <- cofInt(fit2,newdata1,alpha = 0.05)

par(mfrow = c(1,2))
plotR(newdata1$X,ci1$fit,type = "l",ylim = ci1$ylim, 
       xlab = "X", ylab = TeX("$\\mu$"),
       main= "Pred Mod1", lwd = 2)
matlines(newdata1$X,ci1$int, col = "#0082FF", lty=2,type = "l")
points(X,Y,pch = 20,col="red",cex = 1.5)

plotR(newdata1$X,ci2$fit,type = "l",ylim = ci2$ylim, 
       xlab = "X", ylab = TeX("$\\mu$"),
       main= "Pred Mod2", lwd = 2)
matlines(newdata1$X,ci2$int, col = "#0082FF", lty=2,type = "l")
points(X,Y,pch = 20,col="red",cex = 1.5)



par(mfrow = c(1,2))
plotR(X,Y,pch = 20, main = "XvsY", ylab = "Y", xlab = "X")
FIT <-lm( )
plotR(X,log(Y),pch = 20, main = "Xvslog(Y)", ylab = "ln(Y)", xlab = "X")

par(mfrow = c(1,2))
plotQQqres(fit1)
plotQQqres(fit2)

