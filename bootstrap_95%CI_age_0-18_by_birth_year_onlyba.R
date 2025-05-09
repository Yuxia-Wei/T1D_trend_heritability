#supplementary Figure 3(A): bootstrapping process to calculate 95%CIs for heritability (age 0-18) in each birth year based on model 2 (with only ba parameter)
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
#options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')
library(sampling)
library(OpenMx)
library(parallel)
#solution to "there is no package called OpenMx": restart Rstudio after installing OpenMx

################################T1D at age 0-18 years#######################
### Source functions
source( 'W:/C6_Carlsson/Yuxia Wei/do/familial/250114_trend/source_bootfunction_t1d18_onlyba.R' )
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d18_any_allpair.csv' , sep=";" , header=TRUE) 
# Function to use to get estimates from a bootstrap replicate
bootfun <- function( doit=1 ){
  #aeFun( drawfamiliesFun( dat=datWide2 ) )
  aeFun( drawfamiliesFun( dat=datWide ) )
}
# Test that it works
#bootfun( 1 )
#savetest <- lapply( rep(1,3) , bootfun)

# Set up parallel computing clusters (using socketing)
cl <- makeCluster( 3 )
# Step one, set up for analysis
clusterEvalQ( cl , {
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
library(sampling)
library(OpenMx)
source( 'W:/C6_Carlsson/Yuxia Wei/do/familial/250114_trend/source_bootfunction_t1d18_onlyba.R' )

### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d18_any_allpair.csv' , sep=";" , header=TRUE) 
}
)
# Run bootstrap
N <- 1000
t0 <- Sys.time()
aeboot <- parLapply( cl , rep(1,N) , fun=bootfun )
t1 <- Sys.time()
t1-t0
stopCluster( cl )

# Data into useful format
aebootResults <- matrix( NA , length(aeboot) , 64 )
colnames(aebootResults) <- c('h2_1982','e2_1982','h2_1983','e2_1983','h2_1984','e2_1984','h2_1985','e2_1985','h2_1986','e2_1986','h2_1987','e2_1987','h2_1988','e2_1988','h2_1989','e2_1989','h2_1990','e2_1990','h2_1991','e2_1991','h2_1992','e2_1992','h2_1993','e2_1993','h2_1994','e2_1994','h2_1995','e2_1995','h2','e2','h2_1997','e2_1997','h2_1998','e2_1998','h2_1999','e2_1999','h2_2000','e2_2000','h2_2001','e2_2001','h2_2002','e2_2002','h2_2003','e2_2003','h2_2004','e2_2004','h2_2005','e2_2005','h2_2006','e2_2006','h2_2007','e2_2007','h2_2008','e2_2008','h2_2009','e2_2009','h2_2010','e2_2010','h2dif_0085','h2dif_0090','h2dif_0095','h2dif_1085','h2dif_1090','h2dif_1095')
for( i in 1:length(aeboot) ){
  aebootResults[i,] <- unlist( aeboot[[i]] )
}
aebootResults <- as.data.frame( aebootResults )

saveRDS( aebootResults , 'W:/C6_Carlsson/Yuxia Wei/results/familial/trend_250114/bootrun_t1d18_onlyba.Rds' )
aebootResults <- readRDS('W:/C6_Carlsson/Yuxia Wei/results/familial/trend_250114/bootrun_t1d18_onlyba.Rds' )

# Check distribution
windows()
par(mfrow=c(2,3))
hist( aebootResults$h2_1982 )
hist( aebootResults$h2_1983 )
hist( aebootResults$h2_1984 )
hist( aebootResults$h2_1985 )
hist( aebootResults$h2_1986 )
hist( aebootResults$h2_1987 )

windows()
par(mfrow=c(2,3))
hist( aebootResults$h2_1988 )
hist( aebootResults$h2_1989 )
hist( aebootResults$h2_1990 )
hist( aebootResults$h2_1991 )
hist( aebootResults$h2_1992 )
hist( aebootResults$h2_1993 )

windows()
par(mfrow=c(2,3))
hist( aebootResults$h2_1994 )
hist( aebootResults$h2_1995 )
hist( aebootResults$h2 )
hist( aebootResults$h2_1997 )
hist( aebootResults$h2_1998 )
hist( aebootResults$h2_1999 )

windows()
par(mfrow=c(2,3))
hist( aebootResults$h2_2000 )
hist( aebootResults$h2_2001 )
hist( aebootResults$h2_2002 )
hist( aebootResults$h2_2003 )
hist( aebootResults$h2_2004 )
hist( aebootResults$h2_2005 )

windows()
par(mfrow=c(2,3))
hist( aebootResults$h2_2006 )
hist( aebootResults$h2_2007 )
hist( aebootResults$h2_2008 )
hist( aebootResults$h2_2009 )
hist( aebootResults$h2_2010 )

windows()
par(mfrow=c(2,3))
hist( aebootResults$h2dif_0085 )
hist( aebootResults$h2dif_0090 )
hist( aebootResults$h2dif_0095 )
hist( aebootResults$h2dif_1085 )
hist( aebootResults$h2dif_1090 )
hist( aebootResults$h2dif_1095 )

#To get 95% bootstrap intervals
#If using SLSQP as the optimizer:
replace h2_CI_t1d18_onlyba_trend with h2dif_CI_t1d18_onlyba_trend_SLSQP
and h2dif_CI_t1d18_onlyba_trend with h2dif_CI_t1d18_onlyba_trend_SLSQP

#trend of heritability
h2_CI_t1d18_onlyba_trend=data.frame(rbind(quantile( aebootResults$h2_1982 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1983 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1984 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1985 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1986 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1987 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1988 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1989 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1990 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1991 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1992 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1993 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1994 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1995 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1997 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1998 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_1999 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2000 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2001 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2002 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2003 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2004 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2005 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2006 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2007 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2008 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2009 , c(.025,0.5,.975) ),
		quantile( aebootResults$h2_2010 , c(.025,0.5,.975) )))
h2_CI_t1d18_onlyba_trend$Birth_year=seq(1982,2010,1)
names(h2_CI_t1d18_onlyba_trend)=c("h2_lci","h2_median","h2_uci","Birth_year")
#Re-arrange column orders
h2_CI_t1d18_onlyba_trend=subset(h2_CI_t1d18_onlyba_trend,select=c(Birth_year,h2_median,h2_lci,h2_uci))

# trend of heritability difference
h2dif_CI_t1d18_onlyba_trend=data.frame(rbind(quantile( aebootResults$h2dif_0085 , c(.025,0.5,.975) ),
				quantile( aebootResults$h2dif_0090 , c(.025,0.5,.975) ),
				quantile( aebootResults$h2dif_0095 , c(.025,0.5,.975) ),
				quantile( aebootResults$h2dif_1085 , c(.025,0.5,.975) ),
				quantile( aebootResults$h2dif_1090 , c(.025,0.5,.975) ),
				quantile( aebootResults$h2dif_1095 , c(.025,0.5,.975) ) ))
h2dif_CI_t1d18_onlyba_trend$Birth_year=c("2000 vs 1985", "2000 vs 1990", 
								"2000 vs 1995", "2010 vs 1985",
								"2010 vs 1990","2010 vs 1995")
names(h2dif_CI_t1d18_onlyba_trend)=c("h2dif_lci","h2dif_median","h2dif_uci","Birth_year")
#Re-arrange column orders
h2dif_CI_t1d18_onlyba_trend=subset(h2dif_CI_t1d18_onlyba_trend,select=c(Birth_year,h2dif_median,h2dif_lci,h2dif_uci))

h2dif_CI_t1d18_onlyba_trend$CI=paste("(",round(h2dif_CI_t1d18_onlyba_trend$h2dif_lci,digit=3),", ",round(h2dif_CI_t1d18_onlyba_trend$h2dif_uci,digit=3),")",sep="")
#trend of e2
e2_CI_t1d18_onlyba_trend=data.frame(rbind(quantile( aebootResults$e2_1982 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1983 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1984 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1985 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1986 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1987 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1988 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1989 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1990 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1991 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1992 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1993 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1994 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1995 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1997 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1998 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_1999 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2000 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2001 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2002 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2003 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2004 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2005 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2006 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2007 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2008 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2009 , c(.025,0.5,.975) ),
		quantile( aebootResults$e2_2010 , c(.025,0.5,.975) )))
e2_CI_t1d18_onlyba_trend$Birth_year=seq(1982,2010,1)
names(e2_CI_t1d18_onlyba_trend)=c("e2_lci","e2_median","e2_uci","Birth_year")
#Re-arrange column orders
e2_CI_t1d18_onlyba_trend=subset(e2_CI_t1d18_onlyba_trend,select=c(Birth_year,e2_median,e2_lci,e2_uci))

library(openxlsx)
write.xlsx(h2_CI_t1d18_onlyba_trend,file="h2_CI_t1d18_onlyba_trend.xlsx",overwrite=TRUE)
write.xlsx(h2dif_CI_t1d18_onlyba_trend,file="h2dif_CI_t1d18_onlyba_trend.xlsx",overwrite=TRUE)
write.xlsx(e2_CI_t1d18_onlyba_trend,file="e2_CI_t1d18_onlyba_trend.xlsx",overwrite=TRUE)
getwd()

