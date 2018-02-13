a=1
b=2
c=3
rm(list=ls())
library(np); library(stats); library(splines); library(gam)
library(SemiPar); library(KernSmooth); library(locfit); library(lokern); library(lpridge)
library(pspline); library(sfsmisc)

library(dlnm); library(locpol)
data(chicagoNMMAPS)
meantemp=c(); rhum=c(); death=c(); o3=c()

for (i in 1:730){
  begin=1+7*(i-1); end=7+7*(i-1)
  meantemp[i]=mean(chicagoNMMAPS[begin:end,10]) 
  rhum[i]=mean(chicagoNMMAPS[begin:end,12])
  death[i]=sum(chicagoNMMAPS[begin:end,7])
  o3[i]=mean(chicagoNMMAPS[begin:end,14])
}
Newdata=data.frame(meantemp,rhum,death, o3)
Newdata=Newdata[1:574,]
Newdata=Newdata[-446,]
meantemp=Newdata[,1]; rhum=Newdata[,2]; death=Newdata[,3]; o3=Newdata[,4]
plot(meantemp, death)
plot(rhum, death)
plot(o3, death)
play=56
for (i in 1:9){print(c)}
wow=9
