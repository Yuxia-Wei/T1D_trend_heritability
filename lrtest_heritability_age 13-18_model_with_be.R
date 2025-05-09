#Supplementary Table 3: model fitting for likelihood ratio tests: T1D heritability at age 13-18_model 3 (with be parameter) based on all sibling pairs
#OpenMx
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')

#install.packages("OpenMx")
#install.packages("RcppParallel")
#install.packages("polycor")
#install.packages("drgee")
library(OpenMx)
library(polycor)
setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418") 
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d1318_any_allpair.csv', sep="," , header=TRUE)

# Definition variables to be used - Assumed named "sex_1" and "sex_2"
datWide$by_1=(datWide$birth_year_index-1996)/28
datWide$by_2=(datWide$birth_year-1996)/28

names(datWide)[names(datWide) == "men_index"] <- "sex_1"
names(datWide)[names(datWide) == "men"] <- "sex_2"

defNames <- c('by_1','by_2','sex_1','sex_2')
varNames <- c("t1d1318_index","t1d1318")
datfulsib <- datWide[ datWide$Relation=="0010 - Helsyskon" , c(varNames,defNames) ]
# OpenMx wants the variables in a particular format:
datfulsib[,varNames] <- mxFactor( datfulsib[,varNames] , levels=c(0,1) )
### Observed information
### Prevalence
prev_1 <- mean( unlist(datWide[,varNames]) , na.rm=T )
prev_1 
# Threshold
thresh_1 <- qnorm( 1-prev_1 )
thresh_1 

### Modelling heritability 
### Specify AE model
AEmod <- mxModel('AE' , 
# The means: fixed at 0???
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=FALSE , values = 0 , labels=c('mean','mean') , name='expMean') ,
# Thresholds: free to  be estimated 
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=TRUE , values = thresh_1 , labels=c('thrfulsib','thrfulsib') , name='expThrfulsib') ,
# Regression coefficients
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=TRUE , values = .01 , labels='betSex', name='BetSex') ,
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=TRUE , values = .01 , labels='betBy', name='BetBy') ,

# Parameters to be estimated
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.1 ,  name='a') ,
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.01 ,  name='be') ,

# Variance fixed==1, the e is function of a and c
  mxAlgebra( sqrt( 1 - a^2 ) , name='e' ),

# Some parameters for outputs
  mxAlgebra( a^2+e^2 , name='Var'),
  mxAlgebra( a^2/Var , name ='h2'),
  mxAlgebra( e^2/Var , name ='e2'),  
# The expected covariance matrices
  mxData( datfulsib[ , c(varNames,defNames) ] , type = 'raw' ),
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_1') , name='By1') ,
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_2') , name='By2') ,

  mxAlgebra( rbind( cbind(a^2+(e+be*By1)^2 , .5*a^2  ),

                    cbind( .5*a^2 , a^2+(e+be*By2)^2 ) ) , name='expCovfulsib' ) ,
  
### RALF : You don't need submodels when you only have one dataset to read in
# Definition variables - Assumed named "sex_1" and "sex_2"
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=F , labels=c('data.sex_1','data.sex_2') , name='Sex') ,
# Regression model for mean
  mxAlgebra( expMean + BetSex %x% Sex+ BetBy %x% cbind(By1,By2) , name='eMean' ),

# Expected means and covariance
  mxExpectationNormal( covariance = 'expCovfulsib' , means = 'eMean' , thresholds='expThrfulsib' , dimnames = varNames , threshnames=varNames ),
# Fit function
  mxFitFunctionML() 
)
# Fit the model
mxOption( NULL , 'Default optimizer' , 'CSOLNP' ) # May need to try another optimizer
AE_t1d1318_any_onlybe <- mxRun( AEmod , intervals=FALSE )

summary(AE_t1d1318_any_onlybe)
save(AE_t1d1318_any_onlybe, file="AE_t1d1318_any_onlybe.RData") #names to be changed
###########################################comparing different models using likelihood ratio tests###############################################
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d1318_any.RData")
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d1318_any_onlyba.RData")
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d1318_any_nomod.RData")

com1=mxCompare(AE_t1d1318_any, AE_t1d1318_any_onlyba)
com2=mxCompare(AE_t1d1318_any, AE_t1d1318_any_onlybe)
com3=mxCompare(AE_t1d1318_any, AE_t1d1318_any_nomod)
com4=mxCompare(AE_t1d1318_any_onlyba, AE_t1d1318_any_onlybe)
com1
com2
com3
com4
model="Model 1"
a=summary(AE_t1d1318_any)
beta1=a[["parameters"]][["Estimate"]][5]
beta2=a[["parameters"]][["Estimate"]][6]
minus2LL=com1@results[["minus2LL"]][1]
AIC=com1@results[["AIC"]][1]
p_value=1
model1=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model="Model 2"
b=summary(AE_t1d1318_any_onlyba)
beta1=b[["parameters"]][["Estimate"]][5]
beta2=0
minus2LL=com1@results[["minus2LL"]][2]
p_value=com1@results[["p"]][2]
AIC=com1@results[["AIC"]][2]
model2=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model="Model 3"
c=summary(AE_t1d1318_any_onlybe)
beta1=0
beta2=c[["parameters"]][["Estimate"]][5]
minus2LL=com2@results[["minus2LL"]][2]
p_value=com2@results[["p"]][2]
AIC=com2@results[["AIC"]][2]
model3=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model="Model 4"
beta1=0
beta2=0
minus2LL=com3@results[["minus2LL"]][2]
p_value=com3@results[["p"]][2]
AIC=com3@results[["AIC"]][2]
model4=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model_fitting_t1d1318_any=data.frame(rbind(model1,model2,model3,model4))
library(openxlsx)
write.xlsx(model_fitting_t1d1318_any,file="model_fitting_t1d1318_any.xlsx",overwrite=TRUE)
getwd()


