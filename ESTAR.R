setwd("C:/Users/Juan Andres/OneDrive/Fuentes/Codigo R")

#####Transicion suave#####
#Funcion de transicion LSTAR Y ESTAR
L=function(st,gam,c){
  1/(1+exp(-gam*(st-c)))
}
E=function(st,gam,c){
  1-exp(-gam*((st-c)^2))
}

n=1000
ESTAR=matrix(NA,n,1)
ESTAR[1:3]=c(0.2,0.122,0.446)
tt=c()
tt[1:2]=c(0.5,0.2)
for (i in 3:1000) {
  tt[i]=0.5*tt[i-1]+0.3*tt[i-2]+rnorm(1)
}

#tt=seq(0,2,0.002)
#tt=tt[-1]

c=0
for (i in 4:n) {
  ESTAR[i]=((0.5+0.8897*ESTAR[i-1]-0.4858*ESTAR[i-2])*(1-exp(-10*((tt[i-1]-c)^2)))+(0.2+0.7*ESTAR[i-1]-0.6*ESTAR[i-2])*(1-(1-exp(-10*((tt[i-1]-c)^2))))+rnorm(1))
}

plot(as.ts(ESTAR))

assign("tt2",tt)

rm(n,tt,L,i,E,c)

save(list = ls(all.names = TRUE), file = ".simulacionESTAR")
