rm(list=ls())
library(np); library(stats); library(splines); library(gam)
library(SemiPar); library(KernSmooth); library(locfit); library(lokern); library(lpridge)
library(pspline); library(sfsmisc)
library(dlnm); library(locpol)
data(chicagoNMMAPS)

date=chicagoNMMAPS[,1]; time=chicagoNMMAPS[,2]; year=chicagoNMMAPS[,3]
month=chicagoNMMAPS[,4]; day=chicagoNMMAPS[,5]; dow=chicagoNMMAPS[,6]
death=chicagoNMMAPS[,7]; cvd=chicagoNMMAPS[,8]; respd=chicagoNMMAPS[,9]
meantemp=chicagoNMMAPS[,10]; dptp=chicagoNMMAPS[,11]; rhum=chicagoNMMAPS[,12]
pm10=chicagoNMMAPS[,13]; o3=chicagoNMMAPS[,14]

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

# Program real data  "1" vs 3 change points

theta1=seq(17,28, by=1); n1=length(theta1)
theta2=seq(65,82, by=1); n2=length(theta2)
theta3=seq(15,35, by=1); n3=length(theta3)

mat=matrix(,(n2),2)
res1=c()
sum=0
for (i in 1:n2){
  z1=rhum-theta2[i]; z1=z1*(z1>0)
  X=data.frame(meantemp,rhum, death, o3, z1)
  glmT <- glm(death ~ meantemp + rhum + o3 + z1 , data=X, family = poisson())
  res1=sum((resid(glmT, type='pearson'))^2)/568
  xx=cbind(theta2[i], res1)
  sum=sum+1
  mat[sum,]=xx
}
jpSIM=mat[which.min(mat[,2]),1]; print(jpSIM)
z1=rhum-jpSIM[1]; z1=z1*(z1>0)
data=data.frame(meantemp,rhum,death,o3,z1)
glm0 <- glm(death ~ meantemp + rhum + o3 + z1, data=data, family = poisson())


mat=matrix(,(n1*n2*n3),4)
res1=c()
sum=0
for (i in 1:n1){
  for (j in 1:n2){
    for (k in 1:n3){
      z1=meantemp-theta1[i]; z1=z1*(z1>0)
      z2=rhum-theta2[j]; z2=z2*(z2>0)
      z3=o3-theta3[k]; z3=z3*(z3>0)
      
      X=data.frame(meantemp,rhum,death,o3, z1, z2, z3)
      glmT <- glm(death ~ meantemp + rhum + o3 + z1 + z2 +z3 , data=X, family = poisson())
      res1=sum((resid(glmT, type='pearson'))^2)/566
      xx=cbind(theta1[i],theta2[j],theta3[k], res1)
      sum=sum+1
      mat[sum,]=xx
    }}}
jpSIM=mat[which.min(mat[,4]),1:3]; print(jpSIM)

z1=meantemp-jpSIM[1]; z1=z1*(z1>0)
z2=rhum-jpSIM[2]; z2=z2*(z2>0)
z3=o3-jpSIM[3]; z3=z3*(z3>0)
data=data.frame(meantemp,rhum,death, o3,z1,z2,z3)

glm1 <- glm(death ~ meantemp + rhum + o3 + z1 + z2 + z3, data=data, family = poisson())
T0=(sum(resid(glm0, type='pearson')^2)/568)/(sum(resid(glm1, type='pearson')^2)/566); T0
summary(glm0)

# ++++ Gtting T's for permuted data
ResPerm0=resid(glm0, type='pearson')
fitted0=fitted(glm0)
simulation=200

Tper=c()
for (l in 1:simulation){
  death =fitted0+sample(ResPerm0)
  death =cbind(death)
  
  mat=matrix(,n2,2)
  res1=c()
  sum=0
  for (i in 1:n2){
    z1=rhum-theta2[i]; z1=z1*(z1>0)
    X=data.frame(meantemp, rhum, death, o3, z1)
    glmT <- glm(death ~ meantemp + rhum + o3 + z1 , data=X, family = poisson())
    res1=sum(resid(glmT, type='pearson')^2)/568
    xx=cbind(theta2[i], res1)
    sum=sum+1
    mat[sum,]=xx
  }
  jpSIM=mat[which.min(mat[,2]),1]; print(jpSIM)
  z1=rhum-jpSIM[1]; z1=z1*(z1>0)
  data=data.frame(meantemp,rhum,death,o3,z1)
  glm0 <- glm(death ~ meantemp + rhum + o3 + z1, data=data, family = poisson())
  
  
  mat=matrix(,(n1*n2*n3),4)
  res1=c()
  
  sum=0
  for (i in 1:n1){
    for (j in 1:n2){
      for (k in 1:n3){
        z1=meantemp-theta1[i]; z1=z1*(z1>0)
        z2=rhum-theta2[j]; z2=z2*(z2>0)
        z3=o3-theta3[k]; z3=z3*(z3>0)
        
        X=data.frame(death, meantemp, rhum, o3 , z1, z2, z3)
        glmT <- glm(death ~ meantemp + rhum + o3 + z1 + z2 +z3 , data=X, family = poisson())
        res1=sum(resid(glmT, type='pearson')^2)/566
        xx=cbind(theta1[i],theta2[j],theta3[k], res1)
        sum=sum+1
        mat[sum,]=xx
      }}}
  jpSIM=mat[which.min(mat[,4]),1:3]; print(jpSIM)
  
  z1=meantemp-jpSIM[1]; z1=z1*(z1>0)
  z2=rhum-jpSIM[2]; z2=z2*(z2>0)
  z3=o3-jpSIM[3]; z3=z3*(z3>0)
  data=data.frame(meantemp,rhum,death, o3,z1,z2,z3)
  
  glm1 <- glm(death ~ meantemp + rhum + o3 + z1 + z2 + z3, data=data, family = poisson())
  Tper[l]=(sum(resid(glm0, type='pearson')^2)/568)/(sum(resid(glm1, type='pearson')^2)/566)
  print(Tper)
}
p.value13=sum(Tper >= T0)/(simulation); print(p.value13)
# +++++++++++++++++++++++++