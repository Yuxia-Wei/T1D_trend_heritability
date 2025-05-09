#Figure 3 (A), Figure 3 (B): source function for bootstrapping to calculate 95%CIs for heritability (age 13-18) in each birth year based on model 4 (without ba or be parameter)
# Definition variables to be used - Assumed named "sex_1" and "sex_2"
aeFun <- function( datWideIn ){
  datWideIn$by_1=(datWideIn$birth_year_index-1996)/28
  datWideIn$by_2=(datWideIn$birth_year-1996)/28

  names(datWideIn)[names(datWideIn) == "men_index"] <- "sex_1"
  names(datWideIn)[names(datWideIn) == "men"] <- "sex_2"

  defNames <- c('by_1','by_2','sex_1','sex_2')
  varNames <- c("t1d1318_index","t1d1318")
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
    mxMatrix(type='Symm' , nrow=1 , ncol=1 , free=TRUE , values = 0.1 ,  name='a') ,
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
  AE_t1d1318_nomod <- mxRun( AEmod , intervals=TRUE )
  return( c(mxEval(h2,AE_t1d1318_nomod),mxEval(e2,AE_t1d1318_nomod)) )
}

### Function to select random families
drawfamiliesFun <- function( dat ){
  dattemp <-  getdata( dat , 
    cluster( dat , clustername = 'familyid' , size = length(unique(dat$familyid)) ,method = 'srswr' ) 
    )
  dattemp <- dattemp[ rep(row.names(dattemp),dattemp$Replicates) , ]
  return(dattemp)
}


