##Tarea6 Ej 1 y 3
#####################Paquetería
rm(list =ls())
library(ggplot2)   
library(tseries)      
library(GGally)
library(latex2exp)
library(plot3D)
library(dplyr)
library(glmnet)
library(Matrix)
########################Ejercicios###############################
#######################Funciones necesarias######################
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
     colnames(RR) <- c("Valor Singular", paste("var",1:ncol(X)),"Indice de Cond")
     RR<-as_tibble(RR,.name_repair = NULL)
     return(RR)
     }
###########################################################################

###Ridge Regression

ridgeDirecto <- function(X, y, k)solve( t(X) %*% X + k * diag(ncol(X)), t(X) %*% y )
ridgeMC <- function(X,y,k){
           p <- ncol(X)
           n <- ncol(X)
           X <- rbind(X, sqrt(k) * diag(p))
           y <- as.matrix(c(y, rep(0, p)), nr = n)
           betaR <- solve(t(X) %*% X)%*%t(X)%*% y 
           return(betaR)
          }

ridgeDf <- function(X,y,k){
           d<-svd(X)$d
           df<-sapply(k,function(x)sum((d^2)/(x+d^2)))
           beta<-sapply(k,function(x)ridgeDirecto(X,y,x))
           return(list(df=df,beta = beta)) 
           }

Grupos <-function(n.g,n.o){ 
         size<- rep(0,n.g)
         res<- n.g
         sapply(1:n.g, function(i){ 
               size[i]<<-ifelse(n.o/res - floor(n.o/res) <= 0.5,floor(n.o/res),ceiling(n.o/res))
               n.o<<-n.o - size[i]
               res<<-res - 1
             }
           )
         return(size)
         }

TrainTest <- function(X,y,ord,ind){
              n.a<-nrow(X)
              n.g<-length(ind)-1
              Tt<-vector(mode = "list", length = n.g)
              sapply(1:n.g, function(i){
                subse<<-ord[(ind[i]+1):ind[i+1]]  
                XT<<- X[-subse, ]
                mT<<- colMeans(XT)    
                XT<<- apply(XT,2,function(x)x - mean(x))
                yT<<- as.matrix(y[-subse])
                alp<<- mean(yT)
                Xt <<- t(t(X[subse, ]) - mT)
                yt <<- as.matrix(y[subse]) 
                Tt[[i]][[1]] <<- as.matrix(XT)
                Tt[[i]][[2]] <<- as.matrix(yT)
                Tt[[i]][[3]] <<- alp
                Tt[[i]][[4]] <<- matrix(Xt, nc = ncol(X))
                Tt[[i]][[5]] <<- as.matrix(yt)
                }
               ) 
              return(Tt)
             } 

KxVC <- function(X,y,n.g,km,kM,nk = 1e3){
        n<-nrow(X)
        ind<- c(0, cumsum(Grupos(n.g,n)))  
        ord<- sample(1:n)
        VC<-rep(0,nk) 
        eVC <- rep(0, n.g)
        kk<-seq(km,kM, length.out = nk)
        Tt<-TrainTest(X,y,ord,ind) 
        sapply(1:nk,function(i){sapply(1:n.g,function(j){
                bR<<-ridgeDirecto(
                     Tt[[j]][[1]],Tt[[j]][[2]]- Tt[[j]][[3]],kk[i]
                     ) 
                eVC[j]<<- mean(
                    (Tt[[j]][[5]]-Tt[[j]][[3]]-Tt[[j]][[4]] %*% bR)^2
                  )
               }       
             )
             VC[i]<<-mean(eVC)
            }
           )
        return(VC) 
       }

KxVCN<- function(X,y,km = 0, kM = 10,nk = 1e3, nr = 100,n.g = nrow(X)){
         VC<-rep(0,nk)
         sapply(1:nr,function(i) VC<<-VC + KxVC(X,y,n.g,km,kM,nk))
         return(
         data.frame(  
           k = seq(km,kM,length.out = nk),
           VC = VC/nr
          )
         )
        }

###VIFS Calculados de acuerdo a las notas del Dr. Rogelio
ridgeVIF <- function(X, k){
    p <- ncol(X)
    diag(solve(t(X) %*% X + k * diag(p),
         t(X) %*% X) %*% solve(t(X) %*% X + k * diag(p)))
}

############################## Ejercicio 3 ##################################

dat1<-read.table(file.choose(),header = TRUE, sep = ",")

X<-as.matrix(Ces1(dat1[,-6]))
y<-as.matrix(dat1[ ,6] - mean(dat1[,6]))

dat11<-data.frame(X,y)
head(X)

###Hoerl 1975

(bet <- lm(y~.-y-1,data = dat11)$coeff)
(s<-summary(lm(y~.-y-1,data = dat11))$sigma)
(s1<-sum((y - X%*%as.matrix(bet))^2)/(nrow(X) - ncol(X)))

p <- ncol(X)

##cota
(2*s^2)/sum(bet^2)
#\(6.837\cdot 10^{-4}\)

##k
(kHo <- (p*s^2)/sum(bet^2))
# \(1.709\cdot 10^{-3} \)  


km<- 1e-10
kM<- 10
kk<-seq(km,kM,length.out = 1e3)

##Gráfica de la traza de ridge 

A<-ridgeDf(X,y,kk)
(m<-nrow(A$beta))
(Min<-min(apply(A$beta,1,min)))
(Max<-max(apply(A$beta,1,max)))

plot.new()
par(mfrow = c(1,2))
plotR(A$df,A$beta[1,],type = "l",col = 1,lwd = 2,ylim = c(Min,Max)
     ,ylab = TeX("$\\delta_{(k)}$"), xlab = "df") 
sapply(2:m,function(i)lines(A$df,A$beta[i,],col=i,lwd=2))
legend("topleft", 
       legend= c(TeX("$\\hat{\\delta_{1(k)}$"),TeX("$\\delta_{2(k)}$"),TeX("$\\delta_{3(k)}$"), 
                 TeX("$\\delta_{4(k)}$"),TeX("$\\delta_{5(k)}$")),
       col = 1:m, lwd = rep(1,m), box.lwd = 2,cex = 0.90
       ) 
 

plotR(kk,A$beta[1,],type = "l",col = 1,lwd = 2,ylim = c(Min,Max),
      ,ylab = "",xlab = "k")
sapply(2:m,function(i)lines(kk,A$beta[i,],col=i,lwd=2))
legend("topright", 
       legend= c(TeX("$\\delta_{1(k)}$"),TeX("$\\delta_{2(k)}$"),TeX("$\\delta_{3(k)}$"), 
                 TeX("$\\delta_{4(k)}$"),TeX("$\\delta_{5(k)}$")),
       col = 1:m, lwd = rep(1,m), box.lwd = 2,cex = 0.90
       )  
### CV 

kM<-0.5 
kk<-seq(km,kM,length.out = 1e3)

set.seed(123)
system.time(
A <- KxVCN(X,y,km,kM,nr = 10)       
)   

par(mfrow = c(1,1))
plotR(A[,1],A[,2],type = "l", xlab = TeX("k")
      ,ylab = "err", main = "Validación Cruzada", lwd = 2)



min(A[,2])
i<-which.min(A[,2])
(k<- A[i,1])
###0.08458458
abline(v = k, col = "red", lwd = 2)

bR11 <- ridgeDirecto(X,y,kHo) 
bR11
#  [,1]
#  9743.096
#  4860.627
# 11479.145
# -2156.766
# -2540.559
bR12 <- ridgeDirecto(X,y,k) 
bR12
# [,1]
# 6778.872
# 5570.405
# 6968.7470
# 2825.6089
# -139.7006



###VIFS (R^2 en Python). 
round(ridgeVIF(X,kHo),3) 
round(ridgeVIF(X,k),3) 

############################ Ejercicio 4 #################################

dat2<-read.table(file.choose(),header = TRUE, sep = ",")

X<-as.matrix(Ces1(dat2[2:5]))
y<-as.matrix(dat2[ ,6]- mean(dat2[,6]))
km<- 1e-10
kM<- 35
kk<-seq(km,kM,length.out = 1e4)
##Traza de ridge 

A<-ridgeDf(X,y,kk)
(m<-nrow(A$beta))
(Min<-min(apply(A$beta,1,min)))
(Max<-max(apply(A$beta,1,max)))

plot.new()
par(mfrow = c(1,2))
plotR(A$df,A$beta[1,],type = "l",col = 1,lwd = 2,ylim = c(Min,Max)
     ,ylab = TeX("$\\alpha_{(k)}$"), xlab = "df") 
sapply(2:m,function(i)lines(A$df,A$beta[i,],col=i,lwd=2))
legend("topleft", 
       legend= c(TeX("$\\alpha_{1(k)}$"),TeX("$\\alpha_{2(k)}$"),TeX("$\\alpha_{3(k)}$"), 
                 TeX("$\\alpha_{4(k)}$")),
       col = 1:m, lwd = rep(1,m), box.lwd = 2,cex = 0.90
       ) 

plotR(kk,A$beta[1,],type = "l",col = 1,lwd = 2,ylim = c(Min,Max),
      ,ylab = "",xlab = "k")
sapply(2:m,function(i)lines(kk,A$beta[i,],col=i,lwd=2))
legend("topright", 
       legend= c(TeX("$\\alpha_{1(k)}$"),TeX("$\\alpha_{2(k)}$"),TeX("$\\alpha_{3(k)}$"), 
                 TeX("$\\alpha_{4(k)}$")),
       col = 1:m, lwd = rep(1,m), box.lwd = 2,cex = 0.90
       ) 
 
###Hoerl 1975

bet <- lm(y~ X[,1] +  X[,2] + X[,3]+ X[,4]-1)$coeff
s<-summary(lm(y~ X[,1] +  X[,2] + X[,3]+ X[,4]-1))$sigma
p <- ncol(X)


##cota 
(2*s^2)/sum(bet^2)
# 5.811695\cdot 10^{-3}


##
(kHo <- (p*s^2)/sum(bet^2))
# 0.01162339


### CV 

kM<-0.5
kk<-seq(km,kM,length.out = 1e3)

set.seed(123)
system.time(
A <- KxVCN(X,y,km,kM,nr = 10)       
)   


par(mfrow = c(1,1))
plotR(A[,1],A[,2],type = "l", xlab = TeX("k")
      ,ylab = "err", main = "Validación Cruzada", lwd = 2)


min(A[,2])
i<-which.min(A[,2])
(k<- A[i,1])
#0.01001001
abline(v = k, col = "red", lwd = 2)

bR21 <- ridgeMC(X,y,kHo) 
round(bR21,3)
bR22 <- ridgeMC(X,y,k) 
round(bR22,3)

###VIFS (R^2 en Python). 
round(ridgeVIF(X,kHo),3) 
round(ridgeVIF(X,k),3) 