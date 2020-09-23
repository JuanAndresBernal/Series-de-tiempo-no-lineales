rm(list = ls())

setwd("C:/Users/Juan Andres/OneDrive/Fuentes/Codigo R")

load(".simulacionLSTAR")
load(".simulacionESTAR")
load(".non-linearity-test")

jorda_escribano=function(yt,t.variable,n.lags,i.variables,stEzt,c,alpha){
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
    ZT_test2[,i]=ZT[,i]*(st^3)
  }
  ZT_h0=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2))
  rm(ZT_test1,ZT_test2)
  ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT));ZT_test3=matrix(NA,nrow(ZT),ncol(ZT))
  ZT_test4=matrix(NA,nrow(ZT),ncol(ZT))
  for (i in 1:ncol(ZT)) {
    ZT_test1[,i]=ZT[,i]*st
    ZT_test2[,i]=ZT[,i]*(st^2)
    ZT_test3[,i]=ZT[,i]*(st^3)
    ZT_test4[,i]=ZT[,i]*(st^4)
  }
  ZT_h1=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2,ZT_test3,ZT_test4))
  rm(ZT_test1,ZT_test2,ZT_test3,ZT_test4)
  
  test01=lm(yt~., data =ZT_h0)
  test02=lm(yt~., data =ZT_h1)
  
  SSR_00=sum(((summary(test01))[["residuals"]])^2)
  SSR_01=sum(((summary(test02))[["residuals"]])^2)
  
  if(stEzt==TRUE){
    test01=lm(yt~.-1, data =ZT_h0)
    test02=lm(yt~.-1, data =ZT_h1)
    
    SSR_00=sum(((summary(test01))[["residuals"]])^2)
    SSR_01=sum(((summary(test02))[["residuals"]])^2)
  }
  m=n.lags
  ly2=nrow(ZT)
  if(i.variables!=F){
    m=n.lags+length(zt00)
  }
  
  lm_1=((SSR_00-SSR_01)/(2*(m+1)))/(SSR_01/(ly-5*(m+1)))
  p_value1=pf(lm_1,2*(m+1),(ly-5*(m+1)),lower.tail = F)
  
  if(stEzt==TRUE){
    lm_1= ((SSR_00-SSR_01)/(2*m))/(SSR_01/(ly-5*m)) 
    p_value1=pf(lm_1,2*m,(ly-5*m) ,lower.tail = F)
  }
  
  ZT_test1=matrix(NA,nrow(ZT),ncol(ZT));ZT_test2=matrix(NA,nrow(ZT),ncol(ZT))
  for (i in 1:ncol(ZT)) {
    ZT_test1[,i]=ZT[,i]*(st^2)
    ZT_test2[,i]=ZT[,i]*(st^4)
  }
  ZT_h0=as.data.frame(cbind(ZT_ni,ZT_test1,ZT_test2))
  rm(ZT_test1,ZT_test2)

  test1=lm(yt~., data =ZT_h0)
  
  SSR_0=sum(((summary(test1))[["residuals"]])^2)
  SSR_1=sum(((summary(test02))[["residuals"]])^2)
  
  if(stEzt==TRUE){
    test01=lm(yt~.-1, data =ZT_h0)
    test02=lm(yt~.-1, data =ZT_h1)
    
    SSR_00=sum(((summary(test01))[["residuals"]])^2)
    SSR_01=sum(((summary(test02))[["residuals"]])^2)
  }
  m=n.lags
  ly2=nrow(ZT)
  if(i.variables!=F){
    m=n.lags+length(zt00)
  }
  
  lm_2=((SSR_0-SSR_1)/(2*(m+1)))/(SSR_1/(ly-5*(m+1)))
  p_value2=pf(lm_2,2*(m+1),(ly-5*(m+1)),lower.tail = F)
  
  if(stEzt==TRUE){
    lm_2= ((SSR_0-SSR_1)/(2*m))/(SSR_1/(ly-5*m)) 
    p_value2=pf(lm_2,2*m,(ly-5*m) ,lower.tail = F)
  }
  
  if(c==0){
    if(p_value1<p_value2){
      print("Existe evidencia para seleccionar un ESTAR")
     }
    if(p_value1>p_value2){
      print("Existe evidencia para seleccionar un LSTAR")
     }
  }
  
  if(c!=0){
    if(stEzt==T){ 
        if(lm_2>qf(1-alpha,2*m,(ly-5*m))&&lm_1<qf(1-alpha,2*m,(ly-5*m))){
          print("Existe evidencia para seleccionar un LSTAR con c=0")
        }
        if(lm_2<qf(1-alpha,2*m,(ly-5*m))&&lm_1>qf(1-alpha,2*m,(ly-5*m))){
          print("Existe evidencia para seleccionar un LSTAR con c=0")
        }
    }
    
    if(stEzt==F){
        if(lm_2>qf(1-alpha,2*(m+1),(ly-5*(m+1)))&&lm_1<qf(1-alpha,2*(m+1),(ly-5*(m+1)))){
          print("Existe evidencia para seleccionar un LSTAR con c=0")
        }
        if(lm_2<qf(1-alpha,2*(m+1),(ly-5*(m+1)))&&lm_1>qf(1-alpha,2*(m+1),(ly-5*(m+1)))){
          print("Existe evidencia para seleccionar un LSTAR con c=0")
        }
    }
    
  }  
    
  return(list(P_VALUE1=p_value1,P_VALUE2=p_value2,ss1=SSR_00,ss2=SSR_01,ss3=SSR_0,ss4=SSR_1,((summary(test02))[["coefficients"]])))
}



jorda_escribano(yt=LSTAR,t.variable=tt1,n.lags=2,i.variables=F,stEzt=F,c=0,alpha = 0.05)

jorda_escribano(yt=LSTAR2,t.variable=tt1,n.lags=2,i.variables=F,stEzt=F,c=0,alpha = 0.05)

jorda_escribano(yt=ESTAR,t.variable=tt2,n.lags=2,i.variables=F,stEzt=F,c=0,alpha = 0.05)

jorda_escribano(yt=lynx,t.variable=lag_2(lynx,2,4),n.lags=4,i.variables=F,stEzt=T,c=0,alpha = 0.05)



plot(lynx)



lm(rnorm(10)~rnorm(10))


View(lm(rnorm(10)~rnorm(10)))
