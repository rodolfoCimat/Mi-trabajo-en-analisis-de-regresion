##############Funciones 
rm(list =ls())
library(ggplot2)   
library(tseries)      
library(GGally)
library(latex2exp)
library(plot3D)
library(dplyr)
library(Matrix)
library(car)
library(regclass)


#######Plots Chiditas 
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

####Dumies t test 
ttestR<-function(lm,a,np){
        est<-t(as.matrix(a))%*%as.matrix(lm$coefficients) 
        sd<-sqrt(t(as.matrix(a))%*%vcov(lm)%*%as.matrix(a))
        t<-est/sd
        p<- 1- pt(abs(t),df=np)+pt(-abs(t),df=np)
        return(list(sd = sd,tstat=t,pval = p))
      } 
#######Inflience Measures 
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


#####Multicolinearity 


Ces1<-function(A){
    A<-apply(A,2,function(x)(x - mean(x)))
    apply(A,2,function(x)x/(sum(x^2))^(1/2))
    } 

Es<-function(A)apply(A,2,function(x)x/(sum(x^2))^(1/2))


VIF<-function(X)diag(cor(X))

cond<-function(X){
      XS<-Es(cbind(rep(1,nrow(X)),X))
      XSC<-Ces1(X)
      val1<-eigen(t(XS)%*%XS)$values
      val2<-eigen(t(XSC)%*%XSC)$values
      return(list(S = max(val1)/min(val1),CS = max(val2)/min(val2),
                  VS = sort(round(val1,3)), VSC= sort(round(val2,3))   
              )
             ) 
      }

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

####RidgeReg 

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

#linearHypothesis(model, hypothesis.matrix, rhs=NULL,
#test=c("F", "Chisq"), vcov.=NULL,
#white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
#singular.ok=FALSE, ...)

#EXAMPLES

K <-cbind(c(0,0,1,0,0,0),c(0,0,0,0,1,0))
lm7
linearHypothesis(lm7,t(K))

ttestR(lm=lm7,a=c(0,0,0,0,1,-1),np = 9)

####Probelma 1
data<-read.csv(file.choose())
head(data)

X<-data[ ,7:9]
y<-data[ ,6]

mod<-lm(lcpr ~ lrqr + lydr + ltcr, data = data)

PPdat<- data[ , 7:9]
summary(mod)

ggpairs(PPdat,lower= list(continuous = my_fn))


ggpairs( data= data)
(R<-cor(X))
###VIFS 
round( diag(solve(R)), 3)

##Centrada y escalada 
XCE<-Ces1(X)
YCE<-y - mean(y) 
dataC<-data.frame( cbind(XCE,YCE))

##Datos excalados
XE<-Es(cbind(rep(1,nrow(X)),X))


##Numeros de condición 
cond(X)

##Descomposición de varianzas
VarDes(XE)

##Eliminar posiblemente 
summary(mod)

mod1<-lm(YCE~lrqr+ltcr-1, data = dataC) 
summary(mod1)
mod2<-lm(YCE~lydr+ltcr-1, data = dataC) 
summary(mod2)


names(data)
mod2<-lm(lcpr~lydr + ltcr, data = data) 
summary(mod2) 

names(data)
mod3<-lm(lcpr~lydr + ltcr-1, data = data) 
summary(mod3) 
names(X)
cor(X[2:3]) 
diag(solve(cor(X[2:3])))



###modelo Ridge 
##Estim Hoerl

mod1<-lm(YCE~lrqr+lydr+ltcr-1, data = dataC) 
(s<-summary(mod1)$sigma) 
(b<-mod1$coef)
s^2

p<-3

(k<-(p*s^2)/sum(b^2))

##Parametros ridge
delta<-ridgeDirecto(XCE,YCE,k)
delta

###R^2 y VIFS
SSres <- sum((YCE - XCE%*%delta)^2)
SStot<- sum( (YCE)^2)
(R<-1-SSres/SStot)

ridgeVIF(XCE,k)

###data original modelo completo sin modificación  
head(data)
dataOr<-data[ ,2:5]

summary(lm(CPR~ ., data = dataOr))