#Figure 3(A): bootstrapping process to calculate 95%CIs for heritability (age 7-12) in each birth year based on model 4 (without ba or be parameter)
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
#options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')
library(sampling)
library(OpenMx)
library(parallel)
#solution to "there is no package called OpenMx": restart Rstudio after installing OpenMx

################################T1D at age 7-12 years#######################
### Source functions
source( 'W:/C6_Carlsson/Yuxia Wei/do/familial/240418_trend/source_bootfunction_t1d712_nomoderation.R' )
### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d712_any_allpair.csv' , sep="," , header=TRUE) 
# Function to use to get estimates from a bootstrap replicate
bootfun <- function( doit=1 ){
  #aeFun( drawfamiliesFun( dat=datWide2 ) )
  aeFun( drawfamiliesFun( dat=datWide ) )
}
# Test that it works
#bootfun( 1 )
#savetest <- lapply( rep(1,3) , bootfun )

# Set up parallel computing clusters (using socketing)
cl <- makeCluster( 3 )
# Step one, set up for analysis
clusterEvalQ( cl , {
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
library(sampling)
library(OpenMx)
source( 'W:/C6_Carlsson/Yuxia Wei/do/familial/240418_trend/source_bootfunction_t1d712_nomoderation.R' )

### Load data
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d712_any_allpair.csv' , sep="," , header=TRUE) 
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
aebootResults <- matrix( NA , length(aeboot) , 2 )
colnames(aebootResults) <- c('h2','e2')
for( i in 1:length(aeboot) ){
  aebootResults[i,] <- unlist( aeboot[[i]] )
}
aebootResults <- as.data.frame( aebootResults )
# Save or load
saveRDS( aebootResults , 'W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/bootrun_t1d712_no_moderation.Rds' )
aebootResults <- readRDS('W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240418/bootrun_t1d712_no_moderation.Rds' )
# Check distribution
windows()
par(mfrow=c(1,2))
hist( aebootResults$h2)
hist( aebootResults$e2)

# To get 95% bootstrap intervals
quantile( aebootResults$h2 , c(.025,0.5,.975) ) #median 0.66 (95% CI:0.61,0.72)
quantile( aebootResults$e2 , c(.025,0.5,.975) ) #median 0.34(95% CI:0.28, 0.39)
# To get standard errors for parameters
sd( aebootResults$h2)
sd( aebootResults$e2)


