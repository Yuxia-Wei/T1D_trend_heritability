#Figure 1B: cumulative incidence of T1D in birth year 1982, 1990, 2000, and 2010
.libPaths("W:/C6_Carlsson/Yuxia Wei/R packages")
options(repos = 'http://nexus.ki.se/repository/cran.r-project.org')
library(openxlsx)
library(ggsurvfit)
library(survival)

setwd("W:/C6_Carlsson/Yuxia Wei/results/familial/trend_240515") 

data1 <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/mbr_cumulative_inc.csv' , sep="," , header=TRUE)
data1$birth_year=factor(data1$birth_year)
#Figure 1 (B)
  data3=subset(data1,birth_year==1982 | birth_year==1990 |birth_year==2000 |birth_year==2010)
  tiff(filename='birth_year_4yr.tiff',  units="in", width=9, height=4, res=300)
  birth_year=survfit(Surv(age_stop, t1d18_any) ~ birth_year, data =data3) %>% 
    ggsurvfit(linewidth = 0.8,type="risk") +
    ggtitle("(B)")+
    add_confidence_interval() +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(hjust = 0, face = "bold"))+
      scale_color_manual(values = c('#F87067', '#00B832','#5C99FF',"orange"),
                       labels = c('Birth year 1982', 'Birth year 1990','Birth year 2000','Birth year 2010')) +
    scale_fill_manual(values = c('#F87067', '#00B832','#5C99FF',"orange"),
                      labels = c('Birth year 1982', 'Birth year 1990','Birth year 2000','Birth year 2010'))+
  labs(y="Cumulative incidence",x="Age in years") 
  birth_year
  dev.off() 

#numbers in Figure 1(B)
#To obtain cumulative incidence at age 18 years
data4=subset(data3,birth_year!=2010)
KM <- survfit(Surv(age_stop, t1d18_any) ~ birth_year, data = data4)
cum_birth_year_18yr=data.frame(summary(KM, times = 18.99)$surv)
names(cum_birth_year_18yr)="survival"
cum_birth_year_18yr$cumulative_incidence=(1-cum_birth_year_18yr$survival)
cum_birth_year_18yr[1,"birth_year"]=1982 
cum_birth_year_18yr[2,"birth_year"]=1990
cum_birth_year_18yr[3,"birth_year"]=2000
write.xlsx(cum_birth_year_18yr,file="cum_birth_year_18yr.xlsx",overwrite=TRUE)
getwd()
