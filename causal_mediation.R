#Table 1: 
#use causal mediation analysis to quantify the proportion of increasing T1D incidence explined by changing prevalence of environmental factors
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')

library(survival)
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240515") 
library(openxlsx)
#######################maternal smoking during pregnancy#########################
#To check cumulative inciddence in those without missing data on maternal smoking during pregnancy
data_with_nomis_smok <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data_with_nomis_smok$year2000_smok=factor(data_with_nomis_smok$year2000_smok)
data_with_nomis_smok=subset(data_with_nomis_smok,smoking1!=9)

KM <- survfit(Surv(age_stop, t1d18_any) ~ year2000_smok, data = data_with_nomis_smok)
cum_year2000_smok_18yr=data.frame(summary(KM, times = 18.99)$surv)
names(cum_year2000_smok_18yr)="survival"
cum_year2000_smok_18yr$cumulative_incidence=(1-cum_year2000_smok_18yr$survival)
cum_year2000_smok_18yr[1,"year2000_smok"]=1983 #0.005837704 overall and 0.005865067 in those without missing maternal smoking
cum_year2000_smok_18yr[2,"year2000_smok"]=2000 #0.009299934 overall and 0.009448788 in those without missing maternal smoking
write.xlsx(cum_year2000_smok_18yr,file="cum_year2000_smok_18yr_nomis_smok.xlsx",overwrite=TRUE)
getwd()  

#data analysis
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,smoking1!=9)
fitM <- glm(formula=smoking1~year2000_smok+men, data=data , family='binomial' )
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000_smok+men+smoking1, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000_smok"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type="response") #for a binomial model, type="response" gives predicted probability of M given different values of L in different individuals
predT0 <- predict(object=fitT, newdata=data0, type="risk")  #Choices are the linear predictor ("lp"), the risk score exp(lp) ("risk") Ask Tomas: What is this?
data1 <- data
data1[, "year2000_smok"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="response")
predT1 <- predict(object=fitT, newdata=data1, type="risk")  
data00 <- data0
data00[, M] <- 0
predT00 <- predict(object=fitT, newdata=data00, type="risk")
data01 <- data0
data01[, M] <- 1
predT01 <- predict(object=fitT, newdata=data01, type="risk")
data10 <- data1
data10[, M] <- 0
predT10 <- predict(object=fitT, newdata=data10, type="risk")
data11 <- data1
data11[, M] <- 1
predT11 <- predict(object=fitT, newdata=data11, type="risk")

#S(t)=exp(-H(t))   H(t) here is H0(tj)*predT
for(j in 1:nt){ 
tj <- Time[j]  
S0 <- mean(exp(-H0(tj)*predT0))  #survival function when X=0 and M equals to the value it has when X=0 (in the original dataset, n=10000)?  mean: average over the population values for covariates
S00 <- exp(-H0(tj)*predT00) #survival function when X=0 and M=0 (n=1)
S01 <- exp(-H0(tj)*predT01) #survival function when X=0 and M=1 (n=1)
S0M1 <- mean(S00*(1-predM1)+S01*predM1) #survival function when X=0 and M equals to the value it has when X=1, PredM1 is the probability of M=1 when X=1
S10 <- exp(-H0(tj)*predT10) #survival function when X=1 and M=0 
S11 <- exp(-H0(tj)*predT11) #survival function when X=1 and M=1
S1M0 <- mean(S10*(1-predM0)+S11*predM0) #survival function when X=1 and M equals to the value it has when X=0, PredM0 is the probability of M=1 when X=0
S1 <- mean(exp(-H0(tj)*predT1)) #survival function when X=1 and M equals to the value it has when X=1 (in the original dataset)?
est[j, "S0"] <- S0 
est[j, "S0M1"] <- S0M1
est[j, "S1M0"] <- S1M0
est[j, "S1"] <- S1     
}
mediation_smoking=as.data.frame(est)
mediation_smoking$TOTm <- (1- mediation_smoking$S1) - (1- mediation_smoking$S0)
mediation_smoking$NDEm <- (1- mediation_smoking$S1) - (1- mediation_smoking$S0M1)
mediation_smoking$NIEm <- (1- mediation_smoking$S0M1) - (1- mediation_smoking$S0)
mediation_smoking$proportion=mediation_smoking$NIEm/mediation_smoking$TOTm
mediation_smoking[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_smoking,file="mediation_smoking.xlsx",overwrite=TRUE)
getwd()  
#mediation proportion after additionally adjusted for maternal age at delivery and maternal educational level: 1.3%
########Plot####################
mediation_smoking1=subset(mediation_smoking,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_smoking1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1983)", 
  main="(C) Maternal smoking during pregnancy", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))
#Present exact proportions in the footnote
#Mediation proportion in boys: 0.02252570; in girls: 0.05399526
################################

#####################mediation by maternal age at delivery#######################################################
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
fitM <- glm(formula=age_mor~year2000+men, data=data , family='gaussian')
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+age_mor, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type="response") #for a binomial model, type="response" gives predicted probability of M given different values of L in different individuals
predT0 <- predict(object=fitT, newdata=data0, type="risk") 

data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="response")
predT1 <- predict(object=fitT, newdata=data1, type="risk") 

data01 <- data0
data01[, M] <- predM1
predT01 <- predict(object=fitT, newdata=data01, type="risk")

data10 <- data1
data10[, M] <- predM0
predT10 <- predict(object=fitT, newdata=data10, type="risk")

#S(t)=exp(-H(t))   H(t) here is H0(tj)*predT
for(j in 1:nt){ 
	tj <- Time[j]  
	S0 <- mean(exp(-H0(tj)*predT0))  #survival function when X=0 and M equals to the value it has when X=0 (in the original dataset, n=10000)?  mean: average over the population values for covariates
	S0M1 <- mean(exp(-H0(tj)*predT01)) #survival function when X=0 and M equals to the value it has when X=1
	S1M0 <- mean(exp(-H0(tj)*predT10)) #survival function when X=1 and M equals to the value it has when X=0
	S1 <- mean(exp(-H0(tj)*predT1)) #survival function when X=1 and M equals to the value it has when X=1 (in the original dataset)?

	est[j, "S0"] <- S0 
	est[j, "S0M1"] <- S0M1
	est[j, "S1M0"] <- S1M0
	est[j, "S1"] <- S1     
}

mediation_age_mor=as.data.frame(est)
mediation_age_mor$TOTm <- (1- mediation_age_mor$S1) - (1- mediation_age_mor$S0)
mediation_age_mor$NDEm <- (1- mediation_age_mor$S1) - (1- mediation_age_mor$S0M1)
mediation_age_mor$NIEm <- (1- mediation_age_mor$S0M1) - (1- mediation_age_mor$S0)
mediation_age_mor$proportion=mediation_age_mor$NIEm/mediation_age_mor$TOTm
mediation_age_mor[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_age_mor,file="mediation_age_mor.xlsx",overwrite=TRUE)
getwd()

mediation_age_mor1=subset(mediation_age_mor,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_age_mor1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="(A) Maternal age at delivery", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))
#maternal education
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
data=subset(data,edu_3g_mor!=9)
data$edu_3g_mor=as.factor(data$edu_3g_mor)
data$edu_3g_mor2=relevel(data$edu_3g_mor,ref="0")
fitM <- multinom(edu_3g_mor2 ~ year2000 + men, data = data)
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+edu_3g_mor2, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type = "probs") 
predT0 <- predict(object=fitT, newdata=data0, type="risk") 

data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="probs")
predT1 <- predict(object=fitT, newdata=data1, type="risk") 

data00 <- data0
data00[, "edu_3g_mor2"] <- 0
data00$edu_3g_mor2=as.factor(data00$edu_3g_mor2)
predT00 <- predict(object=fitT, newdata=data00, type="risk")

data01 <- data0
data01[, "edu_3g_mor2"] <- 1
data01$edu_3g_mor2=as.factor(data01$edu_3g_mor2)
predT01 <- predict(object=fitT, newdata=data01, type="risk")

data02 <- data0
data02[, "edu_3g_mor2"] <- 2
data02$edu_3g_mor2=as.factor(data02$edu_3g_mor2)
predT02 <- predict(object=fitT, newdata=data02, type="risk")


data10 <- data1
data10[, "edu_3g_mor2"] <- 0
data10$edu_3g_mor2=as.factor(data10$edu_3g_mor2)
predT10 <- predict(object=fitT, newdata=data10, type="risk")

data11 <- data1
data11[, "edu_3g_mor2"] <- 1
data11$edu_3g_mor2=as.factor(data11$edu_3g_mor2)
predT11 <- predict(object=fitT, newdata=data11, type="risk")

data12 <- data1
data12[, "edu_3g_mor2"] <- 2
data12$edu_3g_mor2=as.factor(data12$edu_3g_mor2)
predT12 <- predict(object=fitT, newdata=data12, type="risk")

#S(t)=exp(-H(t))   H(t) here is H0(tj)*predT
predM0_f=data.frame(predM0)
predM1_f=data.frame(predM1)

for(j in 1:nt){ 
tj <- Time[j]  
S0 <- mean(exp(-H0(tj)*predT0))  
S00 <- exp(-H0(tj)*predT00) 
S01 <- exp(-H0(tj)*predT01) 
S02 <- exp(-H0(tj)*predT02) 
S0M1 <- mean(S00*predM1_f[,"X0"]+S01*predM1_f[,"X1"]+S02*predM1_f[,"X2"])  
S10 <- exp(-H0(tj)*predT10) 
S11 <- exp(-H0(tj)*predT11) 
S12 <- exp(-H0(tj)*predT12) 
S1M0 <- mean(S10*predM0_f[,"X0"]+S11*predM0_f[,"X1"]+S12*predM0_f[,"X2"])  
S1 <- mean(exp(-H0(tj)*predT1)) 
est[j, "S0"] <- S0 
est[j, "S0M1"] <- S0M1
est[j, "S1M0"] <- S1M0
est[j, "S1"] <- S1     
}
mediation_edu_mor=as.data.frame(est)
mediation_edu_mor$TOTm <- (1- mediation_edu_mor$S1) - (1- mediation_edu_mor$S0)
mediation_edu_mor$NDEm <- (1- mediation_edu_mor$S1) - (1- mediation_edu_mor$S0M1)
mediation_edu_mor$NIEm <- (1- mediation_edu_mor$S0M1) - (1- mediation_edu_mor$S0)
mediation_edu_mor$proportion=mediation_edu_mor$NIEm/mediation_edu_mor$TOTm
mediation_edu_mor[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_edu_mor,file="mediation_edu_mor.xlsx",overwrite=TRUE)
getwd()
#Plot
mediation_edu_mor1=subset(mediation_edu_mor,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_edu_mor1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="(B) Maternal education", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))

##########################maternal Bbacterial infection#########################
#binary
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
data=subset(data,infection_type!=9)

data[data$infection_type==2,"bacteria_infec_mor"]=1
data[data$infection_type!=2,"bacteria_infec_mor"]=0

fitM <- glm(formula=bacteria_infec_mor~year2000+men, data=data , family='binomial' )
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+bacteria_infec_mor, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type="response") #for a binomial model, type="response" gives predicted probability of M given different values of L in different individuals
predT0 <- predict(object=fitT, newdata=data0, type="risk")  #Choices are the linear predictor ("lp"), the risk score exp(lp) ("risk") Ask Tomas: What is this?
data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="response")
predT1 <- predict(object=fitT, newdata=data1, type="risk")  
data00 <- data0
data00[, M] <- 0
predT00 <- predict(object=fitT, newdata=data00, type="risk")
data01 <- data0
data01[, M] <- 1
predT01 <- predict(object=fitT, newdata=data01, type="risk")
data10 <- data1
data10[, M] <- 0
predT10 <- predict(object=fitT, newdata=data10, type="risk")
data11 <- data1
data11[, M] <- 1
predT11 <- predict(object=fitT, newdata=data11, type="risk")

#S(t)=exp(-H(t))   H(t) here is H0(tj)*predT
for(j in 1:nt){ 
tj <- Time[j]  
S0 <- mean(exp(-H0(tj)*predT0))  #survival function when X=0 and M equals to the value it has when X=0 (in the original dataset, n=10000)?  mean: average over the population values for covariates
S00 <- exp(-H0(tj)*predT00) #survival function when X=0 and M=0 (n=1)
S01 <- exp(-H0(tj)*predT01) #survival function when X=0 and M=1 (n=1)
S0M1 <- mean(S00*(1-predM1)+S01*predM1) #survival function when X=0 and M equals to the value it has when X=1, PredM1 is the probability of M=1 when X=1
S10 <- exp(-H0(tj)*predT10) #survival function when X=1 and M=0 
S11 <- exp(-H0(tj)*predT11) #survival function when X=1 and M=1
S1M0 <- mean(S10*(1-predM0)+S11*predM0) #survival function when X=1 and M equals to the value it has when X=0, PredM0 is the probability of M=1 when X=0
S1 <- mean(exp(-H0(tj)*predT1)) #survival function when X=1 and M equals to the value it has when X=1 (in the original dataset)?
est[j, "S0"] <- S0 
est[j, "S0M1"] <- S0M1
est[j, "S1M0"] <- S1M0
est[j, "S1"] <- S1     
}
mediation_bacteria_infec_mor=as.data.frame(est)
mediation_bacteria_infec_mor$TOTm <- (1- mediation_bacteria_infec_mor$S1) - (1- mediation_bacteria_infec_mor$S0)
mediation_bacteria_infec_mor$NDEm <- (1- mediation_bacteria_infec_mor$S1) - (1- mediation_bacteria_infec_mor$S0M1)
mediation_bacteria_infec_mor$NIEm <- (1- mediation_bacteria_infec_mor$S0M1) - (1- mediation_bacteria_infec_mor$S0)
mediation_bacteria_infec_mor$proportion=mediation_bacteria_infec_mor$NIEm/mediation_bacteria_infec_mor$TOTm
mediation_bacteria_infec_mor[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_bacteria_infec_mor,file="mediation_bacteria_infec_mor.xlsx",overwrite=TRUE)
getwd()  
#Plot
mediation_bacteria_infec_mor1=subset(mediation_bacteria_infec_mor,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_bacteria_infec_mor1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="(D) Maternal bacterial infection during pregnancy", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))
#Present exact proportions in the footnote

########################################gestational age########################################
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
data=subset(data,ga_5g!=9)
data$ga_5g=as.factor(data$ga_5g)
data$ga_5g2=relevel(data$ga_5g,ref="4")
fitM <- multinom(ga_5g2 ~ year2000 + men, data = data)
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+ga_5g2, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type = "probs") 
predT0 <- predict(object=fitT, newdata=data0, type="risk") 

data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="probs")
predT1 <- predict(object=fitT, newdata=data1, type="risk") 

data01 <- data0
data01[, "ga_5g2"] <- 1
data01$ga_5g2=as.factor(data01$ga_5g2)
predT01 <- predict(object=fitT, newdata=data01, type="risk")

data02 <- data0
data02[, "ga_5g2"] <- 2
data02$ga_5g2=as.factor(data02$ga_5g2)
predT02 <- predict(object=fitT, newdata=data02, type="risk")

data03 <- data0
data03[, "ga_5g2"] <- 3
data03$ga_5g2=as.factor(data03$ga_5g2)
predT03 <- predict(object=fitT, newdata=data03, type="risk")

data04 <- data0
data04[, "ga_5g2"] <- 4
data04$ga_5g2=as.factor(data04$ga_5g2)
predT04 <- predict(object=fitT, newdata=data04, type="risk")

data05 <- data0
data05[, "ga_5g2"] <- 5
data05$ga_5g2=as.factor(data05$ga_5g2)
predT05 <- predict(object=fitT, newdata=data05, type="risk")


data11 <- data1
data11[, "ga_5g2"] <- 1
data11$ga_5g2=as.factor(data11$ga_5g2)
predT11 <- predict(object=fitT, newdata=data11, type="risk")

data12 <- data1
data12[, "ga_5g2"] <- 2
data12$ga_5g2=as.factor(data12$ga_5g2)
predT12 <- predict(object=fitT, newdata=data12, type="risk")

data13 <- data1
data13[, "ga_5g2"] <- 3
data13$ga_5g2=as.factor(data13$ga_5g2)
predT13 <- predict(object=fitT, newdata=data13, type="risk")

data14 <- data1
data14[, "ga_5g2"] <- 4
data14$ga_5g2=as.factor(data14$ga_5g2)
predT14 <- predict(object=fitT, newdata=data14, type="risk")

data15 <- data1
data15[, "ga_5g2"] <- 5
data15$ga_5g2=as.factor(data15$ga_5g2)
predT15 <- predict(object=fitT, newdata=data15, type="risk")

predM0_f=data.frame(predM0)
predM1_f=data.frame(predM1)

for(j in 1:nt){ 
  tj <- Time[j]  
  S0 <- mean(exp(-H0(tj)*predT0))  
  S01 <- exp(-H0(tj)*predT01) 
  S02 <- exp(-H0(tj)*predT02) 
  S03 <- exp(-H0(tj)*predT03) 
  S04 <- exp(-H0(tj)*predT04) 
  S05 <- exp(-H0(tj)*predT05) 
  S0M1 <- mean(S01*predM1_f[,"X1"]+S02*predM1_f[,"X2"]+S03*predM1_f[,"X3"]+S04*predM1_f[,"X4"]+S05*predM1_f[,"X5"])  
  S11 <- exp(-H0(tj)*predT11) 
  S12 <- exp(-H0(tj)*predT12) 
  S13 <- exp(-H0(tj)*predT13) 
  S14 <- exp(-H0(tj)*predT14) 
  S15 <- exp(-H0(tj)*predT15) 
  S1M0 <- mean(S11*predM0_f[,"X1"]+S12*predM0_f[,"X2"]+S13*predM0_f[,"X3"]+S14*predM0_f[,"X4"]+S15*predM0_f[,"X5"])  
  S1 <- mean(exp(-H0(tj)*predT1)) 
  est[j, "S0"] <- S0 
  est[j, "S0M1"] <- S0M1
  est[j, "S1M0"] <- S1M0
  est[j, "S1"] <- S1     
}

mediation_ga_5g=as.data.frame(est)
mediation_ga_5g$TOTm <- (1- mediation_ga_5g$S1) - (1- mediation_ga_5g$S0)
mediation_ga_5g$NDEm <- (1- mediation_ga_5g$S1) - (1- mediation_ga_5g$S0M1)
mediation_ga_5g$NIEm <- (1- mediation_ga_5g$S0M1) - (1- mediation_ga_5g$S0)
mediation_ga_5g$proportion=mediation_ga_5g$NIEm/mediation_ga_5g$TOTm
mediation_ga_5g[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_ga_5g,file="mediation_ga_5g.xlsx",overwrite=TRUE)
getwd()
#Plot
mediation_ga_5g1=subset(mediation_ga_5g,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_ga_5g1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="(E) Gestational age", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))

#birth weight
data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
data=subset(data,bw_5g!=9)
data$bw_5g=as.factor(data$bw_5g)
data$bw_5g2=relevel(data$bw_5g,ref="4")
fitM <- multinom(bw_5g2 ~ year2000 + men, data = data)
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+bw_5g2, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type = "probs") 
predT0 <- predict(object=fitT, newdata=data0, type="risk") 

data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="probs")
predT1 <- predict(object=fitT, newdata=data1, type="risk") 

data01 <- data0
data01[, "bw_5g2"] <- 1
data01$bw_5g2=as.factor(data01$bw_5g2)
predT01 <- predict(object=fitT, newdata=data01, type="risk")

data02 <- data0
data02[, "bw_5g2"] <- 2
data02$bw_5g2=as.factor(data02$bw_5g2)
predT02 <- predict(object=fitT, newdata=data02, type="risk")

data03 <- data0
data03[, "bw_5g2"] <- 3
data03$bw_5g2=as.factor(data03$bw_5g2)
predT03 <- predict(object=fitT, newdata=data03, type="risk")

data04 <- data0
data04[, "bw_5g2"] <- 4
data04$bw_5g2=as.factor(data04$bw_5g2)
predT04 <- predict(object=fitT, newdata=data04, type="risk")

data05 <- data0
data05[, "bw_5g2"] <- 5
data05$bw_5g2=as.factor(data05$bw_5g2)
predT05 <- predict(object=fitT, newdata=data05, type="risk")


data11 <- data1
data11[, "bw_5g2"] <- 1
data11$bw_5g2=as.factor(data11$bw_5g2)
predT11 <- predict(object=fitT, newdata=data11, type="risk")

data12 <- data1
data12[, "bw_5g2"] <- 2
data12$bw_5g2=as.factor(data12$bw_5g2)
predT12 <- predict(object=fitT, newdata=data12, type="risk")

data13 <- data1
data13[, "bw_5g2"] <- 3
data13$bw_5g2=as.factor(data13$bw_5g2)
predT13 <- predict(object=fitT, newdata=data13, type="risk")

data14 <- data1
data14[, "bw_5g2"] <- 4
data14$bw_5g2=as.factor(data14$bw_5g2)
predT14 <- predict(object=fitT, newdata=data14, type="risk")

data15 <- data1
data15[, "bw_5g2"] <- 5
data15$bw_5g2=as.factor(data15$bw_5g2)
predT15 <- predict(object=fitT, newdata=data15, type="risk")

predM0_f=data.frame(predM0)
predM1_f=data.frame(predM1)

for(j in 1:nt){ 
  tj <- Time[j]  
  S0 <- mean(exp(-H0(tj)*predT0))  
  S01 <- exp(-H0(tj)*predT01) 
  S02 <- exp(-H0(tj)*predT02) 
  S03 <- exp(-H0(tj)*predT03) 
  S04 <- exp(-H0(tj)*predT04) 
  S05 <- exp(-H0(tj)*predT05) 
  S0M1 <- mean(S01*predM1_f[,"X1"]+S02*predM1_f[,"X2"]+S03*predM1_f[,"X3"]+S04*predM1_f[,"X4"]+S05*predM1_f[,"X5"])  
  S11 <- exp(-H0(tj)*predT11) 
  S12 <- exp(-H0(tj)*predT12) 
  S13 <- exp(-H0(tj)*predT13) 
  S14 <- exp(-H0(tj)*predT14) 
  S15 <- exp(-H0(tj)*predT15) 
  S1M0 <- mean(S11*predM0_f[,"X1"]+S12*predM0_f[,"X2"]+S13*predM0_f[,"X3"]+S14*predM0_f[,"X4"]+S15*predM0_f[,"X5"])  
  S1 <- mean(exp(-H0(tj)*predT1)) 
  est[j, "S0"] <- S0 
  est[j, "S0M1"] <- S0M1
  est[j, "S1M0"] <- S1M0
  est[j, "S1"] <- S1     
}

mediation_bw_5g=as.data.frame(est)
mediation_bw_5g$TOTm <- (1- mediation_bw_5g$S1) - (1- mediation_bw_5g$S0)
mediation_bw_5g$NDEm <- (1- mediation_bw_5g$S1) - (1- mediation_bw_5g$S0M1)
mediation_bw_5g$NIEm <- (1- mediation_bw_5g$S0M1) - (1- mediation_bw_5g$S0)
mediation_bw_5g$proportion=mediation_bw_5g$NIEm/mediation_bw_5g$TOTm
mediation_bw_5g[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_bw_5g,file="mediation_bw_5g.xlsx",overwrite=TRUE)
getwd()
#Plot
mediation_bw_5g1=subset(mediation_bw_5g,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_bw_5g1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="(F) Birth weight", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))

##########################mode of delivery####################################
setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_250507") 
library(survival)
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

data <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_mediation.csv' , sep="," , header=TRUE)
data=subset(data,year2000==0 | year2000==1)
data=subset(data,caesarean1!=9)

fitM <- glm(formula=caesarean1~year2000+men, data=data , family='binomial' )
fitT <- coxph(formula=Surv(age_stop, t1d18_any)~year2000+men+caesarean1, data=data)

Time=seq(0,18.5,0.5)

n <- nrow(data)
M <- as.character(fitM$call$formula[2])
nt <- length(Time)
fit.detail <- coxph.detail(object=fitT)
dH0 <- fit.detail$hazard
H0 <- stepfun(fit.detail$time, c(0, cumsum(dH0)))
est <- matrix(nrow=nt, ncol=4)
rownames(est) <- Time
colnames(est) <- c("S0", "S0M1", "S1M0", "S1")

data0 <- data
data0[, "year2000"] <- 0
predM0 <- predict(object=fitM, newdata=data0, type="response") #for a binomial model, type="response" gives predicted probability of M given different values of L in different individuals
predT0 <- predict(object=fitT, newdata=data0, type="risk")  #Choices are the linear predictor ("lp"), the risk score exp(lp) ("risk") Ask Tomas: What is this?
data1 <- data
data1[, "year2000"] <- 1
predM1 <- predict(object=fitM, newdata=data1, type="response")
predT1 <- predict(object=fitT, newdata=data1, type="risk")  
data00 <- data0
data00[, M] <- 0
predT00 <- predict(object=fitT, newdata=data00, type="risk")
data01 <- data0
data01[, M] <- 1
predT01 <- predict(object=fitT, newdata=data01, type="risk")
data10 <- data1
data10[, M] <- 0
predT10 <- predict(object=fitT, newdata=data10, type="risk")
data11 <- data1
data11[, M] <- 1
predT11 <- predict(object=fitT, newdata=data11, type="risk")

#S(t)=exp(-H(t))   H(t) here is H0(tj)*predT
for(j in 1:nt){ 
tj <- Time[j]  
S0 <- mean(exp(-H0(tj)*predT0))  #survival function when X=0 and M equals to the value it has when X=0 (in the original dataset, n=10000)?  mean: average over the population values for covariates
S00 <- exp(-H0(tj)*predT00) #survival function when X=0 and M=0 (n=1)
S01 <- exp(-H0(tj)*predT01) #survival function when X=0 and M=1 (n=1)
S0M1 <- mean(S00*(1-predM1)+S01*predM1) #survival function when X=0 and M equals to the value it has when X=1, PredM1 is the probability of M=1 when X=1
S10 <- exp(-H0(tj)*predT10) #survival function when X=1 and M=0 
S11 <- exp(-H0(tj)*predT11) #survival function when X=1 and M=1
S1M0 <- mean(S10*(1-predM0)+S11*predM0) #survival function when X=1 and M equals to the value it has when X=0, PredM0 is the probability of M=1 when X=0
S1 <- mean(exp(-H0(tj)*predT1)) #survival function when X=1 and M equals to the value it has when X=1 (in the original dataset)?
est[j, "S0"] <- S0 
est[j, "S0M1"] <- S0M1
est[j, "S1M0"] <- S1M0
est[j, "S1"] <- S1     
}
mediation_caesarean=as.data.frame(est)
mediation_caesarean$TOTm <- (1- mediation_caesarean$S1) - (1- mediation_caesarean$S0)
mediation_caesarean$NDEm <- (1- mediation_caesarean$S1) - (1- mediation_caesarean$S0M1)
mediation_caesarean$NIEm <- (1- mediation_caesarean$S0M1) - (1- mediation_caesarean$S0)
mediation_caesarean$proportion=mediation_caesarean$NIEm/mediation_caesarean$TOTm
mediation_caesarean[,"age"]=seq(0,18.5,0.5)
write.xlsx(mediation_caesarean,file="mediation_caesarean.xlsx",overwrite=TRUE)
getwd()  
#Plot
mediation_caesarean1=subset(mediation_caesarean,select=c(TOTm,NDEm,NIEm))
windows()
matplot(Time, mediation_caesarean1, ylim=c(0, 0.01), col=1:3, xlab="Age in years", ylab="Excess cumulative incidence of T1D (birth year 2000 vs 1982)", 
  main="Caesarean", type="l", lty=1)
legend(x=0, y=0.01, lty=1, col=1:3, legend=c("Total excess cumulative incidence","Natural direct effect","Natural indirect effect"))
#Present exact proportions in the footnote
