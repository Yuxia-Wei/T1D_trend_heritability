#Supplementary Table 4: model fitting for likelihood ratio tests: T1D heritability at age 7-12_model 2(with ba parameter) based on one sibling pair from each family
#OpenMx
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')

#install.packages("OpenMx")
#install.packages("RcppParallel")
#install.packages("polycor")
#install.packages("drgee")
library(OpenMx)
#solution to "there is no package called OpenMx": restart Rstudio after installing OpenMx
setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418") 
#change any to other numbers to run analyses for other birth cohorts
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d712_any.csv', sep="," , header=TRUE)

# Definition variables to be used - Assumed named "sex_1" and "sex_2"
datWide$by_1=(datWide$birth_year_index-1996)/28
datWide$by_2=(datWide$birth_year-1996)/28

names(datWide)[names(datWide) == "men_index"] <- "sex_1"
names(datWide)[names(datWide) == "men"] <- "sex_2"

defNames <- c('by_1','by_2','sex_1','sex_2')
varNames <- c("t1d712_index","t1d712")
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
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.01 ,  name='ba') ,

# Variance fixed==1, the e is function of a and c
  mxAlgebra( sqrt( 1 - a^2 ) , name='e' ),

# Some parameters for outputs
  mxAlgebra( (a+(1982-1996)/28*ba)^2+e^2 , name='Var_1982'),
  mxAlgebra( (a+(1982-1996)/28*ba)^2/Var_1982 , name ='h2_1982'),
  mxAlgebra( e^2/Var_1982 , name ='e2_1982'),

  mxAlgebra( (a+(1983-1996)/28*ba)^2+e^2 , name='Var_1983'),
  mxAlgebra( (a+(1983-1996)/28*ba)^2/Var_1983 , name ='h2_1983'),
  mxAlgebra( e^2/Var_1983 , name ='e2_1983'),

  mxAlgebra( (a+(1984-1996)/28*ba)^2+e^2 , name='Var_1984'),
  mxAlgebra( (a+(1984-1996)/28*ba)^2/Var_1984 , name ='h2_1984'),
  mxAlgebra( e^2/Var_1984 , name ='e2_1984'),

  mxAlgebra( (a+(1985-1996)/28*ba)^2+e^2 , name='Var_1985'),
  mxAlgebra( (a+(1985-1996)/28*ba)^2/Var_1985 , name ='h2_1985'),
  mxAlgebra( e^2/Var_1985 , name ='e2_1985'),

  mxAlgebra( (a+(1986-1996)/28*ba)^2+e^2 , name='Var_1986'),
  mxAlgebra( (a+(1986-1996)/28*ba)^2/Var_1986 , name ='h2_1986'),
  mxAlgebra( e^2/Var_1986 , name ='e2_1986'),

  mxAlgebra( (a+(1987-1996)/28*ba)^2+e^2 , name='Var_1987'),
  mxAlgebra( (a+(1987-1996)/28*ba)^2/Var_1987 , name ='h2_1987'),
  mxAlgebra( e^2/Var_1987 , name ='e2_1987'),

  mxAlgebra( (a+(1988-1996)/28*ba)^2+e^2 , name='Var_1988'),
  mxAlgebra( (a+(1988-1996)/28*ba)^2/Var_1988 , name ='h2_1988'),
  mxAlgebra( e^2/Var_1988 , name ='e2_1988'),

  mxAlgebra( (a+(1989-1996)/28*ba)^2+e^2 , name='Var_1989'),
  mxAlgebra( (a+(1989-1996)/28*ba)^2/Var_1989 , name ='h2_1989'),
  mxAlgebra( e^2/Var_1989 , name ='e2_1989'),

  mxAlgebra( (a+(1990-1996)/28*ba)^2+e^2 , name='Var_1990'),
  mxAlgebra( (a+(1990-1996)/28*ba)^2/Var_1990 , name ='h2_1990'),
  mxAlgebra( e^2/Var_1990 , name ='e2_1990'),

  mxAlgebra( (a+(1991-1996)/28*ba)^2+e^2 , name='Var_1991'),
  mxAlgebra( (a+(1991-1996)/28*ba)^2/Var_1991 , name ='h2_1991'),
  mxAlgebra( e^2/Var_1991 , name ='e2_1991'),

  mxAlgebra( (a+(1992-1996)/28*ba)^2+e^2 , name='Var_1992'),
  mxAlgebra( (a+(1992-1996)/28*ba)^2/Var_1992 , name ='h2_1992'),
  mxAlgebra( e^2/Var_1992 , name ='e2_1992'),

  mxAlgebra( (a+(1993-1996)/28*ba)^2+e^2 , name='Var_1993'),
  mxAlgebra( (a+(1993-1996)/28*ba)^2/Var_1993 , name ='h2_1993'),
  mxAlgebra( e^2/Var_1993 , name ='e2_1993'),

  mxAlgebra( (a+(1994-1996)/28*ba)^2+e^2 , name='Var_1994'),
  mxAlgebra( (a+(1994-1996)/28*ba)^2/Var_1994 , name ='h2_1994'),
  mxAlgebra( e^2/Var_1994 , name ='e2_1994'),

  mxAlgebra( (a+(1995-1996)/28*ba)^2+e^2 , name='Var_1995'),
  mxAlgebra( (a+(1995-1996)/28*ba)^2/Var_1995 , name ='h2_1995'),
  mxAlgebra( e^2/Var_1995 , name ='e2_1995'),

  mxAlgebra( a^2+e^2 , name='Var'),
  mxAlgebra( a^2/Var , name ='h2'),
  mxAlgebra( e^2/Var , name ='e2'),  

  mxAlgebra( (a+(1997-1996)/28*ba)^2+e^2 , name='Var_1997'),
  mxAlgebra( (a+(1997-1996)/28*ba)^2/Var_1997 , name ='h2_1997'),
  mxAlgebra( e^2/Var_1997 , name ='e2_1997'),

  mxAlgebra( (a+(1998-1996)/28*ba)^2+e^2 , name='Var_1998'),
  mxAlgebra( (a+(1998-1996)/28*ba)^2/Var_1998 , name ='h2_1998'),
  mxAlgebra( e^2/Var_1998 , name ='e2_1998'),

  mxAlgebra( (a+(1999-1996)/28*ba)^2+e^2 , name='Var_1999'),
  mxAlgebra( (a+(1999-1996)/28*ba)^2/Var_1999 , name ='h2_1999'),
  mxAlgebra( e^2/Var_1999 , name ='e2_1999'),
  
  mxAlgebra( (a+(2000-1996)/28*ba)^2+e^2 , name='Var_2000'),
  mxAlgebra( (a+(2000-1996)/28*ba)^2/Var_2000 , name ='h2_2000'),
  mxAlgebra( e^2/Var_2000 , name ='e2_2000'),

  mxAlgebra( (a+(2001-1996)/28*ba)^2+e^2 , name='Var_2001'),
  mxAlgebra( (a+(2001-1996)/28*ba)^2/Var_2001 , name ='h2_2001'),
  mxAlgebra( e^2/Var_2001 , name ='e2_2001'),

  mxAlgebra( (a+(2002-1996)/28*ba)^2+e^2 , name='Var_2002'),
  mxAlgebra( (a+(2002-1996)/28*ba)^2/Var_2002 , name ='h2_2002'),
  mxAlgebra( e^2/Var_2002 , name ='e2_2002'),

  mxAlgebra( (a+(2003-1996)/28*ba)^2+e^2 , name='Var_2003'),
  mxAlgebra( (a+(2003-1996)/28*ba)^2/Var_2003 , name ='h2_2003'),
  mxAlgebra( e^2/Var_2003 , name ='e2_2003'),

  mxAlgebra( (a+(2004-1996)/28*ba)^2+e^2 , name='Var_2004'),
  mxAlgebra( (a+(2004-1996)/28*ba)^2/Var_2004 , name ='h2_2004'),
  mxAlgebra( e^2/Var_2004 , name ='e2_2004'),

  mxAlgebra( (a+(2005-1996)/28*ba)^2+e^2 , name='Var_2005'),
  mxAlgebra( (a+(2005-1996)/28*ba)^2/Var_2005 , name ='h2_2005'),
  mxAlgebra( e^2/Var_2005 , name ='e2_2005'),

  mxAlgebra( (a+(2006-1996)/28*ba)^2+e^2 , name='Var_2006'),
  mxAlgebra( (a+(2006-1996)/28*ba)^2/Var_2006 , name ='h2_2006'),
  mxAlgebra( e^2/Var_2006 , name ='e2_2006'),

  mxAlgebra( (a+(2007-1996)/28*ba)^2+e^2 , name='Var_2007'),
  mxAlgebra( (a+(2007-1996)/28*ba)^2/Var_2007 , name ='h2_2007'),
  mxAlgebra( e^2/Var_2007 , name ='e2_2007'),

  mxAlgebra( (a+(2008-1996)/28*ba)^2+e^2 , name='Var_2008'),
  mxAlgebra( (a+(2008-1996)/28*ba)^2/Var_2008 , name ='h2_2008'),
  mxAlgebra( e^2/Var_2008 , name ='e2_2008'),

  mxAlgebra( (a+(2009-1996)/28*ba)^2+e^2 , name='Var_2009'),
  mxAlgebra( (a+(2009-1996)/28*ba)^2/Var_2009 , name ='h2_2009'),
  mxAlgebra( e^2/Var_2009 , name ='e2_2009'),

  mxAlgebra( (a+0.5*ba)^2+e^2 , name='Var_2010'),
  mxAlgebra( (a+0.5*ba)^2/Var_2010 , name ='h2_2010'),
  mxAlgebra( e^2/Var_2010 , name ='e2_2010'),
# The expected covariance matrices
  mxData( datfulsib[ , c(varNames,defNames) ] , type = 'raw' ),
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_1') , name='By1') ,
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_2') , name='By2') ,

  mxAlgebra( rbind( cbind((a+ba*By1)^2+e^2 , .5*(a+ba*By1)*(a+ba*By2)  ),

                    cbind( .5*(a+ba*By1)*(a+ba*By2) , (a+ba*By2)^2+e^2 ) ) , name='expCovfulsib' ) ,
  
### RALF : You don't need submodels when you only have one dataset to read in
# Definition variables - Assumed named "sex_1" and "sex_2"
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=F , labels=c('data.sex_1','data.sex_2') , name='Sex') ,
# Regression model for mean
  mxAlgebra( expMean + BetSex %x% Sex+ BetBy %x% cbind(By1,By2) , name='eMean' ),

# Expected means and covariance
  mxExpectationNormal( covariance = 'expCovfulsib' , means = 'eMean' , thresholds='expThrfulsib' , dimnames = varNames , threshnames=varNames ),
# Fit function
  mxFitFunctionML() 
# RALF : Note that these CI's will be for birthyear==0. Because of the moderation due to ba and be the estimates will differ
#        for different values of birthyear.  
)
# Fit the model
mxOption( NULL , 'Default optimizer' , 'CSOLNP' ) # May need to try another optimizer
AE_t1d712_any_onepair_onlyba <- mxRun( AEmod , intervals=FALSE )
summary(AE_t1d712_any_onepair_onlyba)
save(AE_t1d712_any_onepair_onlyba, file="AE_t1d712_any_onepair_onlyba.RData") #"any" means any childhood-onset T1D regardless of conflicting types of diagnosis

#To export results to an excel file
Birth_year=1982
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1982"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1982"]]@result
year1982=data.frame(cbind(Birth_year,h2,e2))
names(year1982)=c("Birth_year","Heritability","Environment")

Birth_year=1983
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1983"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1983"]]@result
year1983=data.frame(cbind(Birth_year,h2,e2))
names(year1983)=c("Birth_year","Heritability","Environment")

Birth_year=1984
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1984"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1984"]]@result
year1984=data.frame(cbind(Birth_year,h2,e2))
names(year1984)=c("Birth_year","Heritability","Environment")

Birth_year=1985
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1985"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1985"]]@result
year1985=data.frame(cbind(Birth_year,h2,e2))
names(year1985)=c("Birth_year","Heritability","Environment")

Birth_year=1986
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1986"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1986"]]@result
year1986=data.frame(cbind(Birth_year,h2,e2))
names(year1986)=c("Birth_year","Heritability","Environment")

Birth_year=1987
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1987"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1987"]]@result
year1987=data.frame(cbind(Birth_year,h2,e2))
names(year1987)=c("Birth_year","Heritability","Environment")

Birth_year=1988
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1988"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1988"]]@result
year1988=data.frame(cbind(Birth_year,h2,e2))
names(year1988)=c("Birth_year","Heritability","Environment")

Birth_year=1989
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1989"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1989"]]@result
year1989=data.frame(cbind(Birth_year,h2,e2))
names(year1989)=c("Birth_year","Heritability","Environment")

Birth_year=1990
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1990"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1990"]]@result
year1990=data.frame(cbind(Birth_year,h2,e2))
names(year1990)=c("Birth_year","Heritability","Environment")

Birth_year=1991
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1991"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1991"]]@result
year1991=data.frame(cbind(Birth_year,h2,e2))
names(year1991)=c("Birth_year","Heritability","Environment")

Birth_year=1992
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1992"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1992"]]@result
year1992=data.frame(cbind(Birth_year,h2,e2))
names(year1992)=c("Birth_year","Heritability","Environment")

Birth_year=1993
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1993"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1993"]]@result
year1993=data.frame(cbind(Birth_year,h2,e2))
names(year1993)=c("Birth_year","Heritability","Environment")

Birth_year=1994
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1994"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1994"]]@result
year1994=data.frame(cbind(Birth_year,h2,e2))
names(year1994)=c("Birth_year","Heritability","Environment")

Birth_year=1995
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1995"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1995"]]@result
year1995=data.frame(cbind(Birth_year,h2,e2))
names(year1995)=c("Birth_year","Heritability","Environment")

Birth_year=1996
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2"]]@result
year1996=data.frame(cbind(Birth_year,h2,e2))
names(year1996)=c("Birth_year","Heritability","Environment")

Birth_year=1997
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1997"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1997"]]@result
year1997=data.frame(cbind(Birth_year,h2,e2))
names(year1997)=c("Birth_year","Heritability","Environment")

Birth_year=1998
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1998"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1998"]]@result
year1998=data.frame(cbind(Birth_year,h2,e2))
names(year1998)=c("Birth_year","Heritability","Environment")

Birth_year=1999
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_1999"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_1999"]]@result
year1999=data.frame(cbind(Birth_year,h2,e2))
names(year1999)=c("Birth_year","Heritability","Environment")

Birth_year=2000
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2000"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2000"]]@result
year2000=data.frame(cbind(Birth_year,h2,e2))
names(year2000)=c("Birth_year","Heritability","Environment")

Birth_year=2001
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2001"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2001"]]@result
year2001=data.frame(cbind(Birth_year,h2,e2))
names(year2001)=c("Birth_year","Heritability","Environment")

Birth_year=2002
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2002"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2002"]]@result
year2002=data.frame(cbind(Birth_year,h2,e2))
names(year2002)=c("Birth_year","Heritability","Environment")

Birth_year=2003
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2003"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2003"]]@result
year2003=data.frame(cbind(Birth_year,h2,e2))
names(year2003)=c("Birth_year","Heritability","Environment")

Birth_year=2004
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2004"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2004"]]@result
year2004=data.frame(cbind(Birth_year,h2,e2))
names(year2004)=c("Birth_year","Heritability","Environment")

Birth_year=2005
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2005"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2005"]]@result
year2005=data.frame(cbind(Birth_year,h2,e2))
names(year2005)=c("Birth_year","Heritability","Environment")

Birth_year=2006
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2006"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2006"]]@result
year2006=data.frame(cbind(Birth_year,h2,e2))
names(year2006)=c("Birth_year","Heritability","Environment")

Birth_year=2007
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2007"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2007"]]@result
year2007=data.frame(cbind(Birth_year,h2,e2))
names(year2007)=c("Birth_year","Heritability","Environment")

Birth_year=2008
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2008"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2008"]]@result
year2008=data.frame(cbind(Birth_year,h2,e2))
names(year2008)=c("Birth_year","Heritability","Environment")

Birth_year=2009
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2009"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2009"]]@result
year2009=data.frame(cbind(Birth_year,h2,e2))
names(year2009)=c("Birth_year","Heritability","Environment")

Birth_year=2010
h2=AE_t1d712_any_onepair_onlyba@algebras[["h2_2010"]]@result
e2=AE_t1d712_any_onepair_onlyba@algebras[["e2_2010"]]@result
year2010=data.frame(cbind(Birth_year,h2,e2))
names(year2010)=c("Birth_year","Heritability","Environment")

heritability_point_t1d712_any_onepair_onlyba=rbind(year1982,year1983,year1984,year1985, year1986,year1987,year1988,year1989,year1990,year1991,year1992,year1993,year1994,year1995,year1996,year1997,year1998,year1999,year2000,year2001,year2002,year2003,year2004,year2005,year2006,year2007,year2008,year2009,year2010)

library(openxlsx)
write.xlsx(heritability_point_t1d712_any_onepair_onlyba,file="heritability_point_t1d712_any_onepair_onlyba.xlsx",overwrite=TRUE)
