rm(list = ls())

setwd("C:/Users/Juan Andres/OneDrive/Fuentes/Codigo R")

load(".simulacionLSTAR")
load(".simulacionESTAR")

lag_2=function(y,n.lag,n.total.lags){
  if(n.lag==0){
    y=y
  }
  else{
    y=y[-(1:n.lag)]
  }
  y=y[-((length(y)-(n.total.lags-n.lag-1)):length(y))]
  return(y)
}


nonlia_test=function(yt,n.expansion,t.variable,n.lags,i.variables,stEzt){
  if(!(n.expansion==3|n.expansion==4)){
    stop("Numero de expansion no permitido")
  }
  if(!(i.variables==F|is.list(i.variables))){
    stop("i.variables debe ser False o una lista de las variables independientes")
  }
  if(!is.numeric(n.lags)){
    stop("n.lags debe ser numerico(lags de la variable dependiente)")
  }
  if(!is.logical(stEzt)){
    stop("stEzt debe ser un argumento logico")
  }
  if(!is.logical(i.variables)){
    warning("Las variables exogenas deben estar en una lista con la longitud adecuada")
  }
  ly=length(yt)
  yt=as.matrix(yt,ly,1)
  st=as.matrix(t.variable)
  ZT=matrix()
  if(i.variables!=F){
    zt00=i.variables
    zt0=list()
    for (i in 1:length(zt00)) {
      zz=zt00[[i]]
      zz=zz[-((ly-n.lags+1):ly)]
      zt0[[i]]=zz
    }
  }
  yt_0=list()
  for (i in 0:(n.lags-1)) {
    yt_0[[i+1]]=lag_2(yt,i,n.lags)
  }
  yt=yt[-(1:n.lags)]
  zt_0=matrix(1,length(yt),1)
  
  YT=matrix()
  i=1
  repeat{
    if(i==1){
      YT=yt_0[[i]]}
    if(i!=1){
      YT=cbind(YT,yt_0[[i]])  
    }
    if(i==n.lags) break
    i=i+1
  }
  
  if(i.variables==F){
    ZT=cbind(zt_0,YT)
    ZT=as.matrix(ZT)
    ZT_ni=ZT[,-1]
  }
  
  if(i.variables!=F){
    i=1
    repeat{
      if(i==1){
        ZT=cbind(zt_0,YT,zt0[[i]])
      }
      if(i!=1){
        ZT=cbind(ZT,zt0[[i]])  
      }
      if(i==length(zt0)) break
      i=i+1
      ZT=as.matrix(ZT)
      ZT_ni=ZT[,-1]
      rm(i)
    }
  }
  if (stEzt==TRUE){
    ZT=as.matrix(ZT_ni)
  }
  if (n.expansion==3){
    ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT));ZT_test3=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=ZT[,i]*st
      ZT_test2[,i]=ZT[,i]*(st^2)
      ZT_test3[,i]=ZT[,i]*(st^3)
    }
    ZT_test=as.matrix(cbind(ZT,ZT_test1,ZT_test2,ZT_test3))
  }
  
  if (n.expansion==4){
    ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT));ZT_test3=matrix(NA,nrow(ZT),ncol(ZT))
    ZT_test4=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=matrix(ZT[,i]*st,length(yt),1)
      ZT_test2[,i]=matrix(ZT[,i]*(st^2),length(yt),1)
      ZT_test3[,i]=matrix(ZT[,i]*(st^3),length(yt),1)
      ZT_test4[,i]=matrix(ZT[,i]*(st^4),length(yt),1)
    }
    ZT_test=as.matrix(cbind(ZT_ni,ZT_test1,ZT_test2,ZT_test3,ZT_test4))
  }
  test1=lm(yt~., data =as.data.frame(ZT))
  test2=lm(yt~., data =as.data.frame(ZT_test))
  
  SSR_0=sum(((summary(test1))[["residuals"]])^2)
  SSR_1=sum(((summary(test2))[["residuals"]])^2)
  
  if(stEzt==TRUE){
    test1=lm(yt~.-1, data =as.data.frame(ZT))
    test2=lm(((summary(test1))[["residuals"]])~.-1, data =as.data.frame(ZT_test))
    
    SSR_0=sum(((summary(test1))[["residuals"]])^2)
    SSR_1=sum(((summary(test2))[["residuals"]])^2)
  }
  m=n.lags
  ly2=nrow(ZT)
  if(i.variables!=F){
    m=n.lags+length(zt00)
  }
  
  lm_1=((SSR_0-SSR_1)/(3*(m+1)))/(SSR_1/(ly2-4*(m+1)))
  F_statistic=qf(0.95,(3*(m+1)),(ly2-4*(m+1)))
  
  if(stEzt==TRUE){
    lm_1= ((SSR_0-SSR_1)/(3*m))/(SSR_1/(ly2-4*m-1)) 
    F_statistic=qf(0.95,3*m,(ly2-4*m-1))
  }
  
  if(lm_1>=F_statistic){
    print("Se rechaza la hipotesis nula: La serie presenta comportamiento no lineal")
    print(paste("Estadistico de prueba",lm_1,sep = " "))
    print(paste("Estadistico de contraste",F_statistic,sep = " "))
    return(Estadistico_de_prueba=lm_1)
    }
  
  if(lm_1<F_statistic){ 
    print("No se rechaza la hipotesis nula: No hay evidencia para rechazar comportamiento lineal en la serie")  
    print(paste("Estadistico de prueba",lm_1,sep = " "))
    print(paste("Estadistico de contraste",F_statistic,sep = " "))
    return(Estadistico_de_prueba=lm_1)
    }
  
  }


nonlia_test(yt=lynx,n.expansion=3,t.variable=lag_2(lynx,2,4),n.lags=4,i.variables=F,stEzt=T)

nonlia_test(yt=LSTAR,n.expansion=3,t.variable=tt1,n.lags=2,stEzt = F,i.variables=F)

ar1=c()
ar1[1]=0.5

for (i in 2:1000) {
  ar1[i]=0.1*ar1[i-1]+rnorm(1)
} 

plot(ar1,type = "l")

jj=rnorm(1000)
nonlia_test(yt=ar1,n.expansion=4,t.variable=jj,n.lags=1,stEzt = F,i.variables=F)

jj=rnorm(114)
#yt=log(lynx)

plot(yt,type="l")

nonlia_test(yt=LSTAR2,n.expansion=3,t.variable=tt1,n.lags=2,stEzt = F,i.variables=F)

nonlia_test(yt=ESTAR,n.expansion=3,t.variable=tt2,n.lags=2,stEzt = F,i.variables=F)


save(list = ls(),file = ".non-linearity-test")
