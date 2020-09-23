##COINTEGRACION 
#Clase de Reposicion

#Instalar URCA, READXL y DYNLM

#El archivo relaciona 
# 1)ConsumtoT= B1 + B2IngresoT + B3Inversiont + At
#                    (+)           (-)

#Los vectores de cointegracion pueden ser 1 - #de variables.

#2)At= Ct - {B0 + B1Ing + B2Inv}
#3)deltaYt = u1 + Phi111deltaYt_1 + ... + Phi11pdeltaYt_p
#           + Phi121deltaXt_1 + ... + deltaXt_p - alpha1*At_1 + Vt

#4) Vt = si(alpha1*Xt) + ut
# El Beta gorro seria = solve(t(alpha1*Xt)%*%(alpha1*Xt))
##------------------------------------------------------------------------##

data=PIB_cons_inv
y=data$Consumo
inv=data$Inversion
ing=data$PIB
T=length(y)

const=matrix(1,T,1)
x=cbind(const, ing, inv)
x

#viendo Tao 3 y valor critico es para rechazar o no rechazar

plot(y)
summary(ur.df(y, type = c("trend"), lags = 3))
#Raiz Unitaria

plot(ing)
summary(ur.df(ing, type = c("trend"), lags = 7))
#Raiz unitaria 


plot(inv)
summary(ur.df(inv, type = c("trend"), lags =4))
#Raiz Unitaria


#EXISTENCIA DE RAIZ UNITARIA 

#Relacion de consumo y el ingreso:

#Son de largo plazo
betahat=solve(t(x)%*%x)%*%t(x)%*%y
betahat

egorro=y-x%*%betahat

##Prueba de DF para los terminos de perturbacion 

plot(egorro)
summary(ur.df(egorro, type = c("none"), lags = 10))
#Si estan cointegradas. 
#No se pueden evaluar con los VC de DF que estan ahi, sino con los valores criticos
#de Engle - Yoo, hay que buscar la tabla en internet. 
#El valor critico es -3.25 y nuestro estadistico dio -2.3

dy=as.ts(diff(y, lag = 1))
dx1=as.ts(diff(ing, lag = 1))
dx2=as.ts(diff(inv, lag = 1))
e=as.ts(egorro)

#El moñito es de regresion, antes del moñito es la varibale independiente 
#y despues las variables explicativaas
reg1=dynlm(dy ~ L(dy,1)+ L(dy, 2)+ L(dy, 3)+ L(dy, 4)+ L(dx1,1)+ L(dx1, 2)+ L(dx1, 3)+ L(dx1, 4)+ L(dx2,1)+ L(dx2, 2)+ L(dx2, 3)+ L(dx2, 4)+ L(e, 1))
reg1

#encontramos el valor de alpha negativo ( el valor de L(e,1))


#TERCERA ETAPA
residMCE=residuals(reg1)
residMCE

Box.test(residMCE, type = c("Box-Pierce"))
#Tenemos ruido blanco, por lo tanto el MCE esta bien hecho 

x3e=-0.254553*x[6:72,]
#Esta multiplicando por cada una de mis variabels que tengo en X
#empezamos desde 6 por las diferencias

b=solve(t(x3e)%*%x3e)%*%t(x3e)%*%residMCE
b

b3=b+betahat
b3
#Es mas o menos, dependiendo de los signos
#Vector de cointegracion para y 
