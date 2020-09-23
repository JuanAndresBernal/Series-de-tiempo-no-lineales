rm(list = ls())

setwd("C:/Users/Juan Andres/OneDrive/Fuentes/Codigo R")

load(".simulacionLSTAR")
load(".simulacionESTAR")


######El valor del delay que minimize el p value

trans_func_test=function(yt,t.variable,n.lags,i.variables,stEzt){
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
  if(stEzt==F){
    st=st[-(1:n.lags)]
  }
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
  ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=ZT[,i]*st
      ZT_test2[,i]=ZT[,i]*(st^2)
    }
    ZT_h0=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2))
    rm(ZT_test1,ZT_test2)
    ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT));ZT_test3=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=ZT[,i]*st
      ZT_test2[,i]=ZT[,i]*(st^2)
      ZT_test3[,i]=ZT[,i]*(st^3)
    }
    ZT_h1=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2,ZT_test3))
    rm(ZT_test1,ZT_test2,ZT_test3)
    
  test1=lm(yt~., data =ZT_h0)
  test2=lm(((summary(test1))[["residuals"]])~., data =ZT_h1)
  
  SSR_0=sum(((summary(test1))[["residuals"]])^2)
  SSR_1=sum(((summary(test2))[["residuals"]])^2)
  
  if(stEzt==TRUE){
    test1=lm(yt~.-1, data =ZT_h0)
    test2=lm(yt~.-1, data =ZT_h1)
    
    SSR_0=sum(((summary(test1))[["residuals"]])^2)
    SSR_1=sum(((summary(test2))[["residuals"]])^2)
    
  }
  m=n.lags
  ly2=nrow(ZT)
  if(i.variables!=F){
    m=n.lags+length(zt00)
  }
  
  lm_1=((SSR_0-SSR_1)/((m+1)))/(SSR_1/(ly-4*(m+1)))
  F_statistic=qf(0.95,(3*(m+1)),(ly-4*(m+1)))
  
  if(stEzt==TRUE){
    lm_1= ((SSR_0-SSR_1)/(m))/(SSR_1/(ly2-4*m)) 
    F_statistic=qf(0.95,m,(ly2-4*m))
  }
  
  if(lm_1>=F_statistic){
    print("Se rechaza la hipotesis nula: Existe evidencia para rechazar a la familia ESTAR")
    print(lm_1)
    print(F_statistic)
    print(pf(lm_1,m,ly2-4*m))
    temp=list(Estadistico_de_prueba=lm_1,Estadistico_de_contraste=F_statistic,Suma_residuales_H0=SSR_0,Suma_residuales_H1=SSR_1)
    return(temp)
  }
  if(lm_1<F_statistic){ 
    print("No se rechaza la hipotesis nula: No hay evidencia para rechazar a la familiar ESTAR")  
    print(lm_1)
    print(F_statistic)
    print(pf(lm_1,m,ly2-4*m))
  }
  
    ZT_test1=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=ZT[,i]*st
    }
    ZT_h0=as.data.frame(cbind(ZT_ni,ZT_test1))
    rm(ZT_test1)
    ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT))
    for (i in 1:ncol(ZT)) {
      ZT_test1[,i]=ZT[,i]*st
      ZT_test2[,i]=ZT[,i]*(st^2)
    }
    ZT_h1=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2))
    rm(ZT_test1,ZT_test2)
    
    test1=lm(yt~., data =ZT_h0)
    test2=lm(yt~., data =ZT_h1)
    
    SSR_0=sum(((summary(test1))[["residuals"]])^2)
    SSR_1=sum(((summary(test2))[["residuals"]])^2)
    
    if(stEzt==TRUE){
      test1=lm(yt~.-1, data =ZT_h0)
      test2=lm(yt~.-1, data =ZT_h1)
      
      SSR_0=sum(((summary(test1))[["residuals"]])^2)
      SSR_1=sum(((summary(test2))[["residuals"]])^2)
    }
    m=n.lags
    ly2=nrow(ZT)
    if(i.variables!=F){
      m=n.lags+length(zt00)
    }
    
    lm_2=((SSR_0-SSR_1)/((m+1)))/(SSR_1/(ly-3*(m+1)))
    F_statistic=qf(0.95,(3*(m+1)),(ly-3*(m+1)))
    
    if(stEzt==TRUE){
      lm_2=((SSR_0-SSR_1)/(m))/(SSR_1/(ly2-3*m)) 
      F_statistic=qf(0.95,m,(ly2-3*m))
    }
    
    if(lm_2>=F_statistic){
      print("Se rechaza la hipotesis nula: Existe evidencia para rechazar a la familia ESTAR")
      print(lm_2)
      print(F_statistic)
    }
    
    if(lm_2<F_statistic){ 
      print("No se rechaza la hipotesis nula: Se considera evidencia a favor del LSTAR")  
      print(lm_2)
      print(F_statistic)
    }
    
      ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT))
      for (i in 1:ncol(ZT)) {
        ZT_test1[,i]=ZT[,i]*st
        ZT_test2[,i]=ZT[,i]*(st^2)
      }
      ZT_h1=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2))
      rm(ZT_test1,ZT_test2)
      
      test1=lm(yt~., data =as.data.frame(ZT_ni))
      test2=lm(yt~., data =ZT_h1)
      
      SSR_0=sum(((summary(test1))[["residuals"]])^2)
      SSR_1=sum(((summary(test2))[["residuals"]])^2)
      
      if(stEzt==TRUE){
        test1=lm(yt~.-1, data =ZT_h0)
        test2=lm(yt~.-1, data =ZT_h1)
        
        SSR_0=sum(((summary(test1))[["residuals"]])^2)
        SSR_1=sum(((summary(test2))[["residuals"]])^2)
      }
      m=n.lags
      ly2=nrow(ZT)
      if(i.variables!=F){
        m=n.lags+length(zt00)
      }
      
      lm_3=((SSR_0-SSR_1)/((m+1)))/(SSR_1/(ly-2*(m+1)))
      F_statistic=qf(0.95,(3*(m+1)),(ly-2*(m+1)))
      
      if(stEzt==TRUE){
        lm_3=((SSR_0-SSR_1)/(m))/(SSR_1/(ly-2*m)) 
        F_statistic=qf(0.95,m,(ly-2*m))
      }
      
      if(lm_3>=F_statistic){
        print("Se rechaza la hipotesis nula: Soporta la seleccion del modelo LSTAR")
        print(lm_3)
        print(F_statistic)
      }
      if(lm_3<F_statistic){ 
        print("No se rechaza la hipotesis nula: El modelo podria ser un ESTAR")  
        print(lm_3)
        print(F_statistic)
      }
      return(lm_3)
    }



trans_func_test(yt=LSTAR,t.variable=tt,n.lags=2,i.variables=F,stEzt=F)

trans_func_test(yt=ESTAR,t.variable=tt,n.lags=2,i.variables=F,stEzt=F)

trans_func_test(yt=lynx,t.variable=lag_2(lynx,2,4),n.lags=4,i.variables=F,stEzt=T)




