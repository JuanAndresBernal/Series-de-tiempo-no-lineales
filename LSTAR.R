rm(list = ls())
setwd("C:/Users/Juan Andres/OneDrive/Fuentes/Codigo R")

#####Transicion suave#####
#Funcion de transicion logistica

G=function(st){
  1/(1+exp(-10*(st-0.92)))
}


n=1000
LSTAR=matrix(NA,n,1)
LSTAR[1:3]=c(0.7,0.122,0.446)
tt=c()
tt[1]=0.5
for (i in 2:1000) {
  tt[i]=0.5*tt[i-1]+rnorm(1)
}

for (i in 4:n) {
    LSTAR[i]=((0.8897*LSTAR[i-1]-0.4858*LSTAR[i-2])*(1/(1+exp(-10*(tt[i-1]))))
    +(6+0.3*LSTAR[i-1]-0.6*LSTAR[i-2]+0.21*LSTAR[i-3])*(1-(1/(1+exp(-10*(tt[i-1])))))+rnorm(1))
}
plot(as.ts(LSTAR),type = "l")


LSTAR2=matrix(NA,n,1)
LSTAR2[1:3]=c(0.7,0.122,0.446)

for (i in 4:n) {
  LSTAR2[i]=((1+0.8897*LSTAR2[i-1]-0.4858*LSTAR2[i-2])*(1/(1+exp(-10*(tt[i-1]))))
            +(1+0.3*LSTAR2[i-1]-0.6*LSTAR2[i-2]+0.21*LSTAR2[i-3])*(1-(1/(1+exp(-10*(tt[i-1])))))+rnorm(1))
}
plot(LSTAR2,type = "l")


assign("tt1",tt)

rm(G,i,n,tt)

save(list = ls(all.names = TRUE), file = ".simulacionLSTAR")

