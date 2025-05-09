#Supplementary Table 6: point esitmate for T1D heritability in the simulation analysis
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
setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240515") 
#change any to other numbers to run analyses for other birth cohorts
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d18_any_simulated.csv' , sep=";" , header=TRUE) 
#1549357 pairs
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

  mxAlgebra( rbind( cbind(a^2+e^2 , .5*a^2  ),

                    cbind( .5*a^2 , a^2+e^2 ) ) , name='expCovfulsib' ) ,
  
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=F , labels=c('data.sex_1','data.sex_2') , name='Sex') ,
# Regression model for mean
  mxAlgebra( expMean + BetSex %x% Sex+ BetBy %x% cbind(By1,By2) , name='eMean' ),

# Expected means and covariance
  mxExpectationNormal( covariance = 'expCovfulsib' , means = 'eMean' , thresholds='expThrfulsib' , dimnames = varNames , threshnames=varNames ),
# Fit function
  mxFitFunctionML() ,
# Profile likelhood confidence bounds
  mxCI(c('h2','e2'))
)
# Fit the model
mxOption( NULL , 'Default optimizer' , 'CSOLNP' ) 
AE_t1d18_any_nomod_simulated <- mxRun( AEmod , intervals=TRUE ) #h2: 0.59 (95% CI: 0.56, 0.61) if the increasing incidence is 0.004 in the generation code; 0.56[95% CI: 0.53, 0.58] if increasing incidence is 0.005 in the code

summary(AE_t1d18_any_nomod_simulated)
save(AE_t1d18_any_nomod_simulated, file="AE_t1d18_any_nomod_simulated.RData") #names to be changed



