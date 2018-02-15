rm(list=ls())
library(np); library(stats); library(splines); library(gam)
n <- 200
no.simulation=500
#Poisson regression coefficients
beta0 <- 3
beta1 <- 0.5
beta2 <- 0.5
#generate covariate values
x1 <- runif(n=n, min=0, max=1.5)
x2 <- runif(n=n, min=0, max=1.5)
z=x1-0.75; z=z*(z>0)

#compute mu's
mu <- exp(beta0 + beta1*x1 + beta2*x2+3*z)
#generate Y-values
y <- rpois(n=n, lambda=mu)
#data set
data <- data.frame(y=y, x1=x1, x2=x2, z=z)
#plot(data)
plot(x1,y)
# -----------------------------------------------------------------------------
theta <- seq(from=0, to=1.5, by=0.05); n=length(theta)
# -----------------------------------------------------------------------------

# Start the clock!
print("0 vs 1")
fit0<- glm(y ~ x1 + x2, family="poisson")
res0=sum(residuals(fit0)^2) 
ResPerm0=residuals(fit0)
fitted0=fitted(fit0)
res=c()
for (i in 1:n){
  z=x1-theta[i]; z=z*(z>0)
  daty=data.frame(x1, x2, z)
  fit1 <- glm(y ~ x1  + x2 + z, data=daty, family="poisson")
  res[i]=sum(resid(fit1)^2)
  
}
xx=cbind(theta, res); print(" theta res ")
jpSIM=xx[which.min(xx[,2]),1] 
T0=res0/min(res)
print(T0); print(jpSIM)

z=x1-jpSIM; z=z*(z>0)
daty=data.frame(x1, x2, z)
fitF <- glm(y ~ x1  + x2 + z, data=daty, family="poisson")
summary(fitF)
# ----------- Estimation of permuted data
TT=c()
for (j in 1:no.simulation){
  Ny =fitted0+sample(ResPerm0)
  Ny=Ny*(Ny>0); Ny =cbind(Ny)
  fit00<- glm(Ny ~ x1 + x2, family="poisson")
  res00=sum(resid(fit00)^2)
  
  res11=c()
  for (i in 1:n){
    z=x1-theta[i]; z=z*(z>0)
    fit11 <- glm(Ny ~ x1 + z+ x2, , family="poisson")
    res11[i]=sum(resid(fit11)^2)
  }
  TT[j]=res00/min(res11)
  print(j); print(TT[j])
}
TT=c(TT,T0)
p.value111=sum(TT >= T0)/(no.simulation+1)
print(p.value111)

# Stop the clock!
#proc.time() - ptm


x1=rnorm(100)
x2=rnorm(100)
x3=rnorm(100)
x4=rnorm(100)

y=rpois(100, 1)

my.data=data.frame(y,x1,x2,x3,x4)
predictors=colnames(my.data[-1])
target=c("y")

indicators=rbinom(length(predictors), 1, 0.5)

# the more common case
result = formula(paste(target, " ~ ", paste(predictors[indicators == 1], collapse = " + ")))

summ=summary(glm(result,data=my.data,family=poisson()))




























































#data(chicagoNMMAPS)

#date=chicagoNMMAPS[,1]; time=chicagoNMMAPS[,2]; year=chicagoNMMAPS[,3]
#month=chicagoNMMAPS[,4]; day=chicagoNMMAPS[,5]; dow=chicagoNMMAPS[,6]
#death=chicagoNMMAPS[,7]; cvd=chicagoNMMAPS[,8]; respd=chicagoNMMAPS[,9]
#meantemp=chicagoNMMAPS[,10]; dptp=chicagoNMMAPS[,11]; rhum=chicagoNMMAPS[,12]
#pm10=chicagoNMMAPS[,13]; o3=chicagoNMMAPS[,14]

#meantemp=c(); rhum=c(); death=c(); o3=c()
#for (i in 1:730){
#  begin=1+7*(i-1); end=7+7*(i-1)
#  meantemp[i]=mean(chicagoNMMAPS[begin:end,10]) 
#  rhum[i]=mean(chicagoNMMAPS[begin:end,12])
#  death[i]=sum(chicagoNMMAPS[begin:end,7])
#  o3[i]=mean(chicagoNMMAPS[begin:end,14])
#}
#Newdata=data.frame(meantemp,rhum,death, o3)
#Newdata=Newdata[1:574,]
#Newdata=Newdata[-446,]
#meantemp=Newdata[,1]; rhum=Newdata[,2]; death=Newdata[,3]; o3=Newdata[,4]
#plot(meantemp, death)
#plot(rhum, death)
#plot(o3, death)


o3=runif(200)
rhum=runif(200)
z1=rhum-0.5; z1=z1*(z1>0)
meantemp=runif(200)
meanp=exp(1+0.08*o3+0.08*rhum+0.08*meantemp+0.9*z1)#+rnorm(200,0,0.01)
plot(rhum, meanp)

death=rpois(200,meanp)
barplot(table(death))
hist(death)
mean(death); var(death)
plot(meanp,death)
#Three variables are created
# Program real data  "1" vs 3 change points

theta1=seq(0,1, by=0.05); n1=length(theta1)
theta2=seq(0,1, by=0.05); n2=length(theta2)
theta3=seq(0,1, by=0.05); n3=length(theta3)

mat=matrix(,n2,2)
res1=c()
sum=0
for (i in 1:n2){
  z1=rhum-theta2[i]; z1=z1*(z1>0)
  X=data.frame(meantemp,rhum, death, o3, z1)
  glmT <- glm(death ~ meantemp + rhum + o3 + z1 , data=X, family = poisson())
  summary(glmT)
  
  res1=sum((resid(glmT, type='pearson'))^2)/568
  xx=cbind(theta2[i], res1)
  sum=sum+1
  mat[sum,]=xx
}
jpSIM=mat[which.min(mat[,2]),1]; print(jpSIM)
z1=rhum-jpSIM; z1=z1*(z1>0)
data=data.frame(meantemp,rhum,death,o3,z1)
glm0 <- glm(death ~ meantemp + rhum + o3 + z1, data=data, family = poisson())
summary(glm0)


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
# +++++++++++++++++++++++++ ==