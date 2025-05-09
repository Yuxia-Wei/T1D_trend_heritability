#Supplementary Table 3: model fitting for likelihood ratio tests: T1D heritability at age 0-18_model 3 (with be parameter) based on all sibling pairs

#Updated by using type=symm
#OpenMx
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')

#install.packages("OpenMx")
#install.packages("RcppParallel")
#install.packages("polycor")
#install.packages("drgee")
library(OpenMx)
#solution to "there is no package called OpenMx": restart Rstudio after installing OpenMx
library(polycor)
setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418") 
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d18_any_allpair.csv', sep=";" , header=TRUE)

# Definition variables to be used - Assumed named "sex_1" and "sex_2"
datWide$by_1=(datWide$birth_year_index-1996)/28
datWide$by_2=(datWide$birth_year-1996)/28

names(datWide)[names(datWide) == "men_index"] <- "sex_1"
names(datWide)[names(datWide) == "men"] <- "sex_2"

defNames <- c('by_1','by_2','sex_1','sex_2')
varNames <- c("t1d18_index","t1d18")
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
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.3 ,  name='a') ,
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.01 ,  name='be') ,

# Variance fixed==1, the e is function of a and c
  mxAlgebra( sqrt( 1 - a^2 ) , name='e' ),

# Some parameters for outputs
  mxAlgebra( a^2+(e+(1982-1996)/28*be)^2 , name='Var_1982'),
  mxAlgebra( a^2/Var_1982 , name ='h2_1982'),
  mxAlgebra( (e+(1982-1996)/28*be)^2/Var_1982 , name ='e2_1982'),

  mxAlgebra( a^2+(e+(1983-1996)/28*be)^2 , name='Var_1983'),
  mxAlgebra( a^2/Var_1983 , name ='h2_1983'),
  mxAlgebra( (e+(1983-1996)/28*be)^2/Var_1983 , name ='e2_1983'),

  mxAlgebra( a^2+(e+(1984-1996)/28*be)^2 , name='Var_1984'),
  mxAlgebra( a^2/Var_1984 , name ='h2_1984'),
  mxAlgebra( (e+(1984-1996)/28*be)^2/Var_1984 , name ='e2_1984'),

  mxAlgebra( a^2+(e+(1985-1996)/28*be)^2 , name='Var_1985'),
  mxAlgebra( a^2/Var_1985 , name ='h2_1985'),
  mxAlgebra( (e+(1985-1996)/28*be)^2/Var_1985 , name ='e2_1985'),

  mxAlgebra( a^2+(e+(1986-1996)/28*be)^2 , name='Var_1986'),
  mxAlgebra( a^2/Var_1986 , name ='h2_1986'),
  mxAlgebra( (e+(1986-1996)/28*be)^2/Var_1986 , name ='e2_1986'),

  mxAlgebra( a^2+(e+(1987-1996)/28*be)^2 , name='Var_1987'),
  mxAlgebra( a^2/Var_1987 , name ='h2_1987'),
  mxAlgebra( (e+(1987-1996)/28*be)^2/Var_1987 , name ='e2_1987'),

  mxAlgebra( a^2+(e+(1988-1996)/28*be)^2 , name='Var_1988'),
  mxAlgebra( a^2/Var_1988 , name ='h2_1988'),
  mxAlgebra( (e+(1988-1996)/28*be)^2/Var_1988 , name ='e2_1988'),

  mxAlgebra( a^2+(e+(1989-1996)/28*be)^2 , name='Var_1989'),
  mxAlgebra( a^2/Var_1989 , name ='h2_1989'),
  mxAlgebra( (e+(1989-1996)/28*be)^2/Var_1989 , name ='e2_1989'),

  mxAlgebra( a^2+(e+(1990-1996)/28*be)^2 , name='Var_1990'),
  mxAlgebra( a^2/Var_1990 , name ='h2_1990'),
  mxAlgebra( (e+(1990-1996)/28*be)^2/Var_1990 , name ='e2_1990'),

  mxAlgebra( a^2+(e+(1991-1996)/28*be)^2 , name='Var_1991'),
  mxAlgebra( a^2/Var_1991 , name ='h2_1991'),
  mxAlgebra( (e+(1991-1996)/28*be)^2/Var_1991 , name ='e2_1991'),

  mxAlgebra( a^2+(e+(1992-1996)/28*be)^2 , name='Var_1992'),
  mxAlgebra( a^2/Var_1992 , name ='h2_1992'),
  mxAlgebra( (e+(1992-1996)/28*be)^2/Var_1992 , name ='e2_1992'),

  mxAlgebra( a^2+(e+(1993-1996)/28*be)^2 , name='Var_1993'),
  mxAlgebra( a^2/Var_1993 , name ='h2_1993'),
  mxAlgebra( (e+(1993-1996)/28*be)^2/Var_1993 , name ='e2_1993'),

  mxAlgebra( a^2+(e+(1994-1996)/28*be)^2 , name='Var_1994'),
  mxAlgebra( a^2/Var_1994 , name ='h2_1994'),
  mxAlgebra( (e+(1994-1996)/28*be)^2/Var_1994 , name ='e2_1994'),

  mxAlgebra( a^2+(e+(1995-1996)/28*be)^2 , name='Var_1995'),
  mxAlgebra( a^2/Var_1995 , name ='h2_1995'),
  mxAlgebra( (e+(1995-1996)/28*be)^2/Var_1995 , name ='e2_1995'),

  mxAlgebra( a^2+e^2 , name='Var'),
  mxAlgebra( a^2/Var , name ='h2'),
  mxAlgebra( e^2/Var , name ='e2'),  

  mxAlgebra( a^2+(e+(1997-1996)/28*be)^2 , name='Var_1997'),
  mxAlgebra( a^2/Var_1997 , name ='h2_1997'),
  mxAlgebra( (e+(1997-1996)/28*be)^2/Var_1997 , name ='e2_1997'),

  mxAlgebra( a^2+(e+(1998-1996)/28*be)^2 , name='Var_1998'),
  mxAlgebra( a^2/Var_1998 , name ='h2_1998'),
  mxAlgebra( (e+(1998-1996)/28*be)^2/Var_1998 , name ='e2_1998'),

  mxAlgebra( a^2+(e+(1999-1996)/28*be)^2 , name='Var_1999'),
  mxAlgebra( a^2/Var_1999 , name ='h2_1999'),
  mxAlgebra( (e+(1999-1996)/28*be)^2/Var_1999 , name ='e2_1999'),
  
  mxAlgebra( a^2+(e+(2000-1996)/28*be)^2 , name='Var_2000'),
  mxAlgebra( a^2/Var_2000 , name ='h2_2000'),
  mxAlgebra( (e+(2000-1996)/28*be)^2/Var_2000 , name ='e2_2000'),

  mxAlgebra( a^2+(e+(2001-1996)/28*be)^2 , name='Var_2001'),
  mxAlgebra( a^2/Var_2001 , name ='h2_2001'),
  mxAlgebra( (e+(2001-1996)/28*be)^2/Var_2001 , name ='e2_2001'),

  mxAlgebra( a^2+(e+(2002-1996)/28*be)^2 , name='Var_2002'),
  mxAlgebra( a^2/Var_2002 , name ='h2_2002'),
  mxAlgebra( (e+(2002-1996)/28*be)^2/Var_2002 , name ='e2_2002'),

  mxAlgebra( a^2+(e+(2003-1996)/28*be)^2 , name='Var_2003'),
  mxAlgebra( a^2/Var_2003 , name ='h2_2003'),
  mxAlgebra( (e+(2003-1996)/28*be)^2/Var_2003 , name ='e2_2003'),

  mxAlgebra( a^2+(e+(2004-1996)/28*be)^2 , name='Var_2004'),
  mxAlgebra( a^2/Var_2004 , name ='h2_2004'),
  mxAlgebra( (e+(2004-1996)/28*be)^2/Var_2004 , name ='e2_2004'),

  mxAlgebra( a^2+(e+(2005-1996)/28*be)^2 , name='Var_2005'),
  mxAlgebra( a^2/Var_2005 , name ='h2_2005'),
  mxAlgebra( (e+(2005-1996)/28*be)^2/Var_2005 , name ='e2_2005'),

  mxAlgebra( a^2+(e+(2006-1996)/28*be)^2 , name='Var_2006'),
  mxAlgebra( a^2/Var_2006 , name ='h2_2006'),
  mxAlgebra( (e+(2006-1996)/28*be)^2/Var_2006 , name ='e2_2006'),

  mxAlgebra( h2_2006-h2_1985,name='h2dif_0685' ),
  mxAlgebra( h2_2006-h2_1990,name='h2dif_0690' ),
  mxAlgebra( h2_2006-h2_1995,name='h2dif_0695' ),

  mxAlgebra( h2_2000-h2_1985,name='h2dif_0085' ),
  mxAlgebra( h2_2000-h2_1990,name='h2dif_0090' ),
  mxAlgebra( h2_2000-h2_1995,name='h2dif_0095' ),
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
AE_t1d18_any_onlybe <- mxRun( AEmod , intervals=TRUE )
#All fit attempts resulted in errors - check starting values or model specification
summary(AE_t1d18_any_onlybe)
save(AE_t1d18_any_onlybe, file="AE_t1d18_any_onlybe.RData") #names to be changed

##################comparing different models and export parameters#########
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d18_any.RData")
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d18_any_onlyba.RData")
load("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/AE_t1d18_any_nomod.RData")

com1=mxCompare(AE_t1d18_any, AE_t1d18_any_onlyba)
com2=mxCompare(AE_t1d18_any, AE_t1d18_any_onlybe)
com3=mxCompare(AE_t1d18_any, AE_t1d18_any_nomod)
com4=mxCompare(AE_t1d18_any_onlyba, AE_t1d18_any_onlybe)
com1
com2
com3
com4
model="Model 1"
a=summary(AE_t1d18_any)
beta1=a[["parameters"]][["Estimate"]][5]
beta2=a[["parameters"]][["Estimate"]][6]
minus2LL=com1@results[["minus2LL"]][1]
AIC=com1@results[["AIC"]][1]
p_value=1
model1=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model="Model 2"
b=summary(AE_t1d18_any_onlyba)
beta1=b[["parameters"]][["Estimate"]][5]
beta2=0
minus2LL=com1@results[["minus2LL"]][2]
p_value=com1@results[["p"]][2]
AIC=com1@results[["AIC"]][2]
model2=cbind(model,beta1,beta2,minus2LL,p_value,AIC)

model="Model 3"
c=summary(AE_t1d18_any_onlybe)
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

model_fitting_t1d18_any=rbind(model1,model2,model3,model4)
model_fitting_t1d18_any=data.frame(model_fitting_t1d18_any)
library(openxlsx)
write.xlsx(model_fitting_t1d18_any,file="model_fitting_t1d18_any.xlsx",overwrite=TRUE)
getwd()