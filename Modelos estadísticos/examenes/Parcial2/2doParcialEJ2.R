data<-read.csv(file.choose())
head(data)


X<-data[ ,-11]
head(X)
y<-data[ ,11]

mod<-lm(FAT~., data = data)

PPdat<- data.frame(X)
summary(mod)

ggpairs(PPdat,lower= list(continuous = my_fn))


(R<-cor(X))
###VIFS 
round( diag(solve(R)), 3)

##Centrada y escalada 
XCE<-Ces1(X)
YCE<-y - mean(y) 
dataC<-data.frame(cbind(XCE,YCE))

##Datos excalados
XE<-Es(cbind(rep(1,nrow(X)),X))


##Numeros de condición 
cond(X)

##
VarDes(XE)[ ,2:13]


###  - c(CWT, WTWAT, DEP) 
modC<- lm(YCE~.-WTWAT-CWT-DEP-1 ,data = dataC)
summary(modC)

###Modelo anterior en Datos originales  
modC<- lm(FAT~.-WTWAT-DEP-CWT,data = data)
summary(modC)


###Modelo anterior en Datos originales sin intercepto 
modC<- lm(FAT~.-WTWAT-DEP-CWT-1,data = data)
summary(modC)

###VIFS
modC<- lm(FAT~.-WTWAT-DEP-CWT,data = data)
###
vif(modC) 


###Modelo con mayor poder predictivo 
modC<- lm(FAT~.-WTWAT-DEP-CWT-DPSL-LEA-1,data = data)
summary(modC)

###VIFS
modC<- lm(FAT~.-WTWAT-DEP-CWT-DPSL-LEA,data = data)
###
vif(modC) 


############################Ridge 
###k Hoerl
modC<- lm(YCE~.-1,data = dataC)
(s<-summary(modC)$sigma) 
(b<-modC$coef)
s^2
p<-10
(k<-(p*s^2)/sum(b^2))


XCE
YCE
##Parametros ridge
delta<-ridgeDirecto(XCE,YCE,k)
delta <- as.matrix(delta)

###R^2 y VIFS
SSres <- sum((YCE - XCE%*%delta)^2)
SStot<- sum( (YCE)^2)
(R<-1-SSres/SStot)

ridgeVIF(XCE,k)