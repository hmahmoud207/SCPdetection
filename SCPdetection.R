# Libraries needed
rm(list=ls())
library(np); library(stats); library(splines); library(gam)
# Simulated data
nn <- 100 # sample zise
no.simulation=10000 # number of permuted data sets for the permutation test
#Poisson regression coefficients
beta0 <- 3
beta1 <- 0.5
beta2 <- 0.5
beta3 <- 0.5
#generate covariate values
x1 <- runif(n=nn, min=0, max=1.5)
x2 <- runif(n=nn, min=0, max=1.5)
x3 <- runif(n=nn, min=0, max=1.5)
z1=x1-0.75; z1=z1*(z1>0)
z2=x2-0.75; z2=z2*(z2>0)

#compute mu's
mu <- exp(beta0 + beta1*x1 + beta2*x2 + beta3*x3 + 3*z1 + 2*z2)
#generate Y-values
y <- rpois(n=nn, lambda=mu)

# ----------------------------------------------------
theta <- seq(from=0, to=1.5, by=0.05); n=length(theta) # Change points grid search
# ----------------------------------------------------

indicatorsN=c(1,1,1,0,0,0)
for (s in 1:3){# Loop for testing the sifnificance of simultaneous change points

for (f in 1:2){# Loop for calculating test statistic value for the original data, T0
  if (f==1){indicators=indicatorsN}
  if (f==2){indicators=c(1,1,1,1,1,1)}
  
  res=c(); it=0
  RES=matrix(NA,length(theta)^3,nn+4); FITT=matrix(NA,length(theta)^3,nn+1)
  for (i in 1:n){
    z1=x1-theta[i]; z1=z1*(z1>0)
    
    for (j in 1:n){
      z2=x2-theta[j]; z2=z2*(z2>0)
      
      for (k in 1:n){
        z3=x3-theta[k]; z3=z3*(z3>0)
        
        it=it+1
        
        my.data=data.frame(y,x1,x2,x3,z1,z2,z3)
        predictors=colnames(my.data[-1])
        target=c("y")
        result = formula(paste(target, " ~ ", paste(predictors[indicators == 1], collapse = " + ")))
        fit1=glm(result,data=my.data,family="poisson")
        
        res[it]=sum(resid(fit1)^2) 
        RES[it,]=c(residuals(fit1), theta[i], theta[j], theta[k], res[it])
        FITT[it,]=c(fitted(fit1),res[it]) 
      }}}
  
  if (f==1) {
    jpSIM10=RES[which.min(RES[,(nn+4)]),(nn+1):(nn+3)]
    ResPerm0=RES[which.min(RES[,(nn+4)]),1:nn]
    fitted0=FITT[which.min(FITT[,(nn+1)]),1:nn]
    res0=min(res)} 
  if (f==2) {
    jpSIM20=RES[which.min(RES[,(nn+4)]),(nn+1):(nn+3)]
    res1=min(res)} 
}
T0=res0/res1

# ----------- PERMUTATION TESTING 
TT=c()
for (jj in 1:no.simulation){# Loop for creating N permuted data sets  
  
  Ny =fitted0+sample(ResPerm0)
  Ny=Ny*(Ny>0); Ny =cbind(Ny)
  
  
  for (f in 1:2){# Loop for calculating test statistic values for permuted data sets
    if (f==1){indicators=indicatorsN}
    if (f==2){indicators=c(1,1,1,1,1,1)}
    
    res=c(); it=0; xx=matrix(NA,length(theta)^3,3)
    RES=matrix(NA,length(theta)^3,(nn+4)); FITT=matrix(NA,length(theta)^3,(nn+1))
    for (i in 1:n){
      z1=x1-theta[i]; z1=z1*(z1>0)
      
      for (j in 1:n){
        z2=x2-theta[j]; z2=z2*(z2>0)
        
        for (k in 1:n){
          z3=x3-theta[j]; z3=z3*(z3>0)
          
          it=it+1
          
          my.data=data.frame(Ny,x1,x2,x3,z1,z2,z3)
          predictors=colnames(my.data[-1])
          target=c("Ny")
          result = formula(paste(target, " ~ ", paste(predictors[indicators == 1], collapse = " + ")))
          fit1=glm(result,data=my.data,family="poisson")
          
          res[it]=sum(resid(fit1)^2) 
          RES[it,]=c(residuals(fit1), theta[i], theta[j], theta[k], res[it])
          FITT[it,]=c(fitted(fit1),res[it]) 
        }}}
    
    if (f==1) {
      jpSIM1=RES[which.min(RES[,(nn+4)]),(nn+1):(nn+3)]
      ResPerm0=RES[which.min(RES[,(nn+4)]),1:nn]
      fitted0=FITT[which.min(FITT[,(nn+1)]),1:nn]
      res00=min(res)} 
    if (f==2) {
      jpSIM2=RES[which.min(RES[,(nn+4)]),(nn+1):(nn+3)]
      res11=min(res)} 
  }
  TT[jj]=res00/res11
}
TT=c(TT,T0); p.valueF=sum(TT >= T0)/(no.simulation + 1); print(p.valueF)

if (p.valueF > 0.05) break
if (sum(indicatorsN)==5) break
if (p.valueF < 0.05) {indicatorsN[s+3]=1}
}

# Declare the Simultaneous change points detected and fit the final model  
if (p.valueF >0.05){
    Change.point=jpSIM10; Variable=c("X1","X2","X3")
    print("Simultaneous change points detected are:")
    print(data.frame(Variable,Change.point))
    z1=x1-jpSIM10[1]; z1=z1*(z1>0)
    z2=x2-jpSIM10[2]; z2=z2*(z2>0)
    z3=x3-jpSIM10[3]; z3=z3*(z3>0)
    my.data=data.frame(y,x1,x2,x3,z1,z2,z3)
    predictors=colnames(my.data[-1])
    target=c("y")
    result = formula(paste(target, " ~ ", paste(predictors[indicatorsN == 1], collapse = " + ")))
    fitF=glm(result,data=my.data,family="poisson")
    summary(fitF)} else {
      
      print("Simultaneous change points detected are:")
      Change.point=jpSIM20; Variable=c("X1","X2","X3")
      print(data.frame(Variable,Change.point))
      z1=x1-jpSIM20[1]; z1=z1*(z1>0)
      z2=x2-jpSIM20[2]; z2=z2*(z2>0)
      z3=x3-jpSIM20[3]; z3=z3*(z3>0)
      my.data=data.frame(y,x1,x2,x3,z1,z2,z3)
      predictors=colnames(my.data[-1])
      target=c("y")
      result = formula(paste(target, " ~ ", paste(predictors[indicators == 1], collapse = " + ")))
      fitF=glm(result,data=my.data,family="poisson")
      summary(fitF)}