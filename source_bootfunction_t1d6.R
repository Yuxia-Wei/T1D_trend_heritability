#Figure 3(B): source function for bootstrapping to calculate 95%CIs for heritability (age 0-6) in each birth year based on model 1 (with ba and be parameter)
aeFun <- function( datWideIn ){
datWideIn$by_1=(datWideIn$birth_year_index-1996)/28
datWideIn$by_2=(datWideIn$birth_year-1996)/28

names(datWideIn)[names(datWideIn) == "men_index"] <- "sex_1"
names(datWideIn)[names(datWideIn) == "men"] <- "sex_2"

defNames <- c('by_1','by_2','sex_1','sex_2')
varNames <- c("t1d6_index","t1d6")
datfulsib <- datWideIn[ datWideIn$Relation=="0010 - Helsyskon" , c(varNames,defNames) ]
# OpenMx wants the variables in a particular format:
datfulsib[,varNames] <- mxFactor( datfulsib[,varNames] , levels=c(0,1) )
### Observed information
### Prevalence
prev_1 <- mean( unlist(datWideIn[,varNames]) , na.rm=T )
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
  mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.01 ,  name='be') ,

# Variance fixed==1, the e is function of a and c
  mxAlgebra( sqrt( 1 - a^2 ) , name='e' ),

# Some parameters for outputs
  mxAlgebra( (a+(1982-1996)/28*ba)^2+(e+(1982-1996)/28*be)^2 , name='Var_1982'),
  mxAlgebra( (a+(1982-1996)/28*ba)^2/Var_1982 , name ='h2_1982'),
  mxAlgebra( (e+(1982-1996)/28*be)^2/Var_1982 , name ='e2_1982'),

  mxAlgebra( (a+(1983-1996)/28*ba)^2+(e+(1983-1996)/28*be)^2 , name='Var_1983'),
  mxAlgebra( (a+(1983-1996)/28*ba)^2/Var_1983 , name ='h2_1983'),
  mxAlgebra( (e+(1983-1996)/28*be)^2/Var_1983 , name ='e2_1983'),

  mxAlgebra( (a+(1984-1996)/28*ba)^2+(e+(1984-1996)/28*be)^2 , name='Var_1984'),
  mxAlgebra( (a+(1984-1996)/28*ba)^2/Var_1984 , name ='h2_1984'),
  mxAlgebra( (e+(1984-1996)/28*be)^2/Var_1984 , name ='e2_1984'),

  mxAlgebra( (a+(1985-1996)/28*ba)^2+(e+(1985-1996)/28*be)^2 , name='Var_1985'),
  mxAlgebra( (a+(1985-1996)/28*ba)^2/Var_1985 , name ='h2_1985'),
  mxAlgebra( (e+(1985-1996)/28*be)^2/Var_1985 , name ='e2_1985'),

  mxAlgebra( (a+(1986-1996)/28*ba)^2+(e+(1986-1996)/28*be)^2 , name='Var_1986'),
  mxAlgebra( (a+(1986-1996)/28*ba)^2/Var_1986 , name ='h2_1986'),
  mxAlgebra( (e+(1986-1996)/28*be)^2/Var_1986 , name ='e2_1986'),

  mxAlgebra( (a+(1987-1996)/28*ba)^2+(e+(1987-1996)/28*be)^2 , name='Var_1987'),
  mxAlgebra( (a+(1987-1996)/28*ba)^2/Var_1987 , name ='h2_1987'),
  mxAlgebra( (e+(1987-1996)/28*be)^2/Var_1987 , name ='e2_1987'),

  mxAlgebra( (a+(1988-1996)/28*ba)^2+(e+(1988-1996)/28*be)^2 , name='Var_1988'),
  mxAlgebra( (a+(1988-1996)/28*ba)^2/Var_1988 , name ='h2_1988'),
  mxAlgebra( (e+(1988-1996)/28*be)^2/Var_1988 , name ='e2_1988'),

  mxAlgebra( (a+(1989-1996)/28*ba)^2+(e+(1989-1996)/28*be)^2 , name='Var_1989'),
  mxAlgebra( (a+(1989-1996)/28*ba)^2/Var_1989 , name ='h2_1989'),
  mxAlgebra( (e+(1989-1996)/28*be)^2/Var_1989 , name ='e2_1989'),

  mxAlgebra( (a+(1990-1996)/28*ba)^2+(e+(1990-1996)/28*be)^2 , name='Var_1990'),
  mxAlgebra( (a+(1990-1996)/28*ba)^2/Var_1990 , name ='h2_1990'),
  mxAlgebra( (e+(1990-1996)/28*be)^2/Var_1990 , name ='e2_1990'),

  mxAlgebra( (a+(1991-1996)/28*ba)^2+(e+(1991-1996)/28*be)^2 , name='Var_1991'),
  mxAlgebra( (a+(1991-1996)/28*ba)^2/Var_1991 , name ='h2_1991'),
  mxAlgebra( (e+(1991-1996)/28*be)^2/Var_1991 , name ='e2_1991'),

  mxAlgebra( (a+(1992-1996)/28*ba)^2+(e+(1992-1996)/28*be)^2 , name='Var_1992'),
  mxAlgebra( (a+(1992-1996)/28*ba)^2/Var_1992 , name ='h2_1992'),
  mxAlgebra( (e+(1992-1996)/28*be)^2/Var_1992 , name ='e2_1992'),

  mxAlgebra( (a+(1993-1996)/28*ba)^2+(e+(1993-1996)/28*be)^2 , name='Var_1993'),
  mxAlgebra( (a+(1993-1996)/28*ba)^2/Var_1993 , name ='h2_1993'),
  mxAlgebra( (e+(1993-1996)/28*be)^2/Var_1993 , name ='e2_1993'),

  mxAlgebra( (a+(1994-1996)/28*ba)^2+(e+(1994-1996)/28*be)^2 , name='Var_1994'),
  mxAlgebra( (a+(1994-1996)/28*ba)^2/Var_1994 , name ='h2_1994'),
  mxAlgebra( (e+(1994-1996)/28*be)^2/Var_1994 , name ='e2_1994'),

  mxAlgebra( (a+(1995-1996)/28*ba)^2+(e+(1995-1996)/28*be)^2 , name='Var_1995'),
  mxAlgebra( (a+(1995-1996)/28*ba)^2/Var_1995 , name ='h2_1995'),
  mxAlgebra( (e+(1995-1996)/28*be)^2/Var_1995 , name ='e2_1995'),

  mxAlgebra( a^2+e^2 , name='Var'),
  mxAlgebra( a^2/Var , name ='h2'),
  mxAlgebra( e^2/Var , name ='e2'),  

  mxAlgebra( (a+(1997-1996)/28*ba)^2+(e+(1997-1996)/28*be)^2 , name='Var_1997'),
  mxAlgebra( (a+(1997-1996)/28*ba)^2/Var_1997 , name ='h2_1997'),
  mxAlgebra( (e+(1997-1996)/28*be)^2/Var_1997 , name ='e2_1997'),

  mxAlgebra( (a+(1998-1996)/28*ba)^2+(e+(1998-1996)/28*be)^2 , name='Var_1998'),
  mxAlgebra( (a+(1998-1996)/28*ba)^2/Var_1998 , name ='h2_1998'),
  mxAlgebra( (e+(1998-1996)/28*be)^2/Var_1998 , name ='e2_1998'),

  mxAlgebra( (a+(1999-1996)/28*ba)^2+(e+(1999-1996)/28*be)^2 , name='Var_1999'),
  mxAlgebra( (a+(1999-1996)/28*ba)^2/Var_1999 , name ='h2_1999'),
  mxAlgebra( (e+(1999-1996)/28*be)^2/Var_1999 , name ='e2_1999'),
  
  mxAlgebra( (a+(2000-1996)/28*ba)^2+(e+(2000-1996)/28*be)^2 , name='Var_2000'),
  mxAlgebra( (a+(2000-1996)/28*ba)^2/Var_2000 , name ='h2_2000'),
  mxAlgebra( (e+(2000-1996)/28*be)^2/Var_2000 , name ='e2_2000'),

  mxAlgebra( (a+(2001-1996)/28*ba)^2+(e+(2001-1996)/28*be)^2 , name='Var_2001'),
  mxAlgebra( (a+(2001-1996)/28*ba)^2/Var_2001 , name ='h2_2001'),
  mxAlgebra( (e+(2001-1996)/28*be)^2/Var_2001 , name ='e2_2001'),

  mxAlgebra( (a+(2002-1996)/28*ba)^2+(e+(2002-1996)/28*be)^2 , name='Var_2002'),
  mxAlgebra( (a+(2002-1996)/28*ba)^2/Var_2002 , name ='h2_2002'),
  mxAlgebra( (e+(2002-1996)/28*be)^2/Var_2002 , name ='e2_2002'),

  mxAlgebra( (a+(2003-1996)/28*ba)^2+(e+(2003-1996)/28*be)^2 , name='Var_2003'),
  mxAlgebra( (a+(2003-1996)/28*ba)^2/Var_2003 , name ='h2_2003'),
  mxAlgebra( (e+(2003-1996)/28*be)^2/Var_2003 , name ='e2_2003'),

  mxAlgebra( (a+(2004-1996)/28*ba)^2+(e+(2004-1996)/28*be)^2 , name='Var_2004'),
  mxAlgebra( (a+(2004-1996)/28*ba)^2/Var_2004 , name ='h2_2004'),
  mxAlgebra( (e+(2004-1996)/28*be)^2/Var_2004 , name ='e2_2004'),

  mxAlgebra( (a+(2005-1996)/28*ba)^2+(e+(2005-1996)/28*be)^2 , name='Var_2005'),
  mxAlgebra( (a+(2005-1996)/28*ba)^2/Var_2005 , name ='h2_2005'),
  mxAlgebra( (e+(2005-1996)/28*be)^2/Var_2005 , name ='e2_2005'),

  mxAlgebra( (a+(2006-1996)/28*ba)^2+(e+(2006-1996)/28*be)^2 , name='Var_2006'),
  mxAlgebra( (a+(2006-1996)/28*ba)^2/Var_2006 , name ='h2_2006'),
  mxAlgebra( (e+(2006-1996)/28*be)^2/Var_2006 , name ='e2_2006'),

  mxAlgebra( (a+(2007-1996)/28*ba)^2+(e+(2007-1996)/28*be)^2 , name='Var_2007'),
  mxAlgebra( (a+(2007-1996)/28*ba)^2/Var_2007 , name ='h2_2007'),
  mxAlgebra( (e+(2007-1996)/28*be)^2/Var_2007 , name ='e2_2007'),

  mxAlgebra( (a+(2008-1996)/28*ba)^2+(e+(2008-1996)/28*be)^2 , name='Var_2008'),
  mxAlgebra( (a+(2008-1996)/28*ba)^2/Var_2008 , name ='h2_2008'),
  mxAlgebra( (e+(2008-1996)/28*be)^2/Var_2008 , name ='e2_2008'),

  mxAlgebra( (a+(2009-1996)/28*ba)^2+(e+(2009-1996)/28*be)^2 , name='Var_2009'),
  mxAlgebra( (a+(2009-1996)/28*ba)^2/Var_2009 , name ='h2_2009'),
  mxAlgebra( (e+(2009-1996)/28*be)^2/Var_2009 , name ='e2_2009'),

  mxAlgebra( (a+0.5*ba)^2+(e+0.5*be)^2 , name='Var_2010'),
  mxAlgebra( (a+0.5*ba)^2/Var_2010 , name ='h2_2010'),
  mxAlgebra( (e+0.5*be)^2/Var_2010 , name ='e2_2010'),

  mxAlgebra( h2_2000-h2_1985,name='h2dif_0085' ),
  mxAlgebra( h2_2000-h2_1990,name='h2dif_0090' ),
  mxAlgebra( h2_2000-h2_1995,name='h2dif_0095' ),

  mxAlgebra( h2_2010-h2_1985,name='h2dif_1085' ),
  mxAlgebra( h2_2010-h2_1990,name='h2dif_1090' ),
  mxAlgebra( h2_2010-h2_1995,name='h2dif_1095' ),

# The expected covariance matrices
  mxData( datfulsib[ , c(varNames,defNames) ] , type = 'raw' ),
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_1') , name='By1') ,
  mxMatrix(type='Full' , nrow=1 , ncol=1 , free=F , labels=c('data.by_2') , name='By2') ,

  mxAlgebra( rbind( cbind((a+ba*By1)^2+(e+be*By1)^2 , .5*(a+ba*By1)*(a+ba*By2)  ),

                    cbind( .5*(a+ba*By1)*(a+ba*By2) , (a+ba*By2)^2+(e+be*By2)^2 ) ) , name='expCovfulsib' ) ,
  
# Definition variables - Assumed named "sex_1" and "sex_2"
  mxMatrix(type='Full' , nrow=1 , ncol=2 , free=F , labels=c('data.sex_1','data.sex_2') , name='Sex') ,
# Regression model for mean
  mxAlgebra(expMean + BetSex %x% Sex+ BetBy %x% cbind(By1,By2) , name='eMean' ),

# Expected means and covariance
  mxExpectationNormal(covariance = 'expCovfulsib' , means = 'eMean' , thresholds='expThrfulsib' , dimnames = varNames , threshnames=varNames ),
# Fit function
  mxFitFunctionML() 
)
# Fit the model
mxOption( NULL , 'Default optimizer' , 'CSOLNP' ) # May need to try another optimizer
#mxOption( NULL , 'Default optimizer' , 'SLSQP' ) # try another optimizer

AE_t1d6 <- mxRun( AEmod , intervals=TRUE )

return( c(mxEval(h2_1982,AE_t1d6),mxEval(e2_1982,AE_t1d6),mxEval(h2_1983,AE_t1d6),mxEval(e2_1983,AE_t1d6),mxEval(h2_1984,AE_t1d6),mxEval(e2_1984,AE_t1d6),mxEval(h2_1985,AE_t1d6),mxEval(e2_1985,AE_t1d6),mxEval(h2_1986,AE_t1d6),mxEval(e2_1986,AE_t1d6),mxEval(h2_1987,AE_t1d6),mxEval(e2_1987,AE_t1d6),mxEval(h2_1988,AE_t1d6),mxEval(e2_1988,AE_t1d6),mxEval(h2_1989,AE_t1d6),mxEval(e2_1989,AE_t1d6),mxEval(h2_1990,AE_t1d6),mxEval(e2_1990,AE_t1d6),mxEval(h2_1991,AE_t1d6),mxEval(e2_1991,AE_t1d6),mxEval(h2_1992,AE_t1d6),mxEval(e2_1992,AE_t1d6),mxEval(h2_1993,AE_t1d6),mxEval(e2_1993,AE_t1d6),mxEval(h2_1994,AE_t1d6),mxEval(e2_1994,AE_t1d6),mxEval(h2_1995,AE_t1d6),mxEval(e2_1995,AE_t1d6),mxEval(h2,AE_t1d6),mxEval(e2,AE_t1d6),mxEval(h2_1997,AE_t1d6),mxEval(e2_1997,AE_t1d6),mxEval(h2_1998,AE_t1d6),mxEval(e2_1998,AE_t1d6),mxEval(h2_1999,AE_t1d6),mxEval(e2_1999,AE_t1d6),mxEval(h2_2000,AE_t1d6),mxEval(e2_2000,AE_t1d6),mxEval(h2_2001,AE_t1d6),mxEval(e2_2001,AE_t1d6),mxEval(h2_2002,AE_t1d6),mxEval(e2_2002,AE_t1d6),mxEval(h2_2003,AE_t1d6),mxEval(e2_2003,AE_t1d6),mxEval(h2_2004,AE_t1d6),mxEval(e2_2004,AE_t1d6),mxEval(h2_2005,AE_t1d6),mxEval(e2_2005,AE_t1d6),mxEval(h2_2006,AE_t1d6),mxEval(e2_2006,AE_t1d6),mxEval(h2_2007,AE_t1d6),mxEval(e2_2007,AE_t1d6),mxEval(h2_2008,AE_t1d6),mxEval(e2_2008,AE_t1d6),mxEval(h2_2009,AE_t1d6),mxEval(e2_2009,AE_t1d6),mxEval(h2_2010,AE_t1d6),mxEval(e2_2010,AE_t1d6),mxEval(h2dif_0085,AE_t1d6),mxEval(h2dif_0090,AE_t1d6),mxEval(h2dif_0095,AE_t1d6),mxEval(h2dif_1085,AE_t1d6),mxEval(h2dif_1090,AE_t1d6),mxEval(h2dif_1095,AE_t1d6)) )
}

### Function to select random families
drawfamiliesFun <- function( dat ){
  dattemp <-  getdata( dat , 
    cluster( dat , clustername = 'familyid' , size = length(unique(dat$familyid)) ,method = 'srswr' ) 
    )
  dattemp <- dattemp[ rep(row.names(dattemp),dattemp$Replicates) , ]
  return(dattemp)
}


