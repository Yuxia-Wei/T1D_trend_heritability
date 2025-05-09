/*Figure 4: Associations of environmental factors (maternal age, education and smoking) with T1D in the traditional cohort and sibling comparison analyses*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\edu_3g_par.csv",clear /*including those meeting exclusion criteria*/
sort lopnr_barn
tempfile edu_3g_par
save `edu_3g_par'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
sort lopnr_barn
keep if birth_year<2011

merge 1:1 lopnr_barn using `edu_3g_par',keepusing(edu_3g_mor edu_3g_far) keep(1 3) nogen
/* 2926490 matched; 0 not matched*/
tab edu_3g_mor,missing
tab edu_3g_far,missing

capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_age_edu_bmi_smok19822010.log",replace

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   

stcox c.age_mor_sd b2.edu_3g_mor b2.edu_3g_far c.bmi_mor_sd i.imputed_bmi i.smoking1 i.men i.birth_year i.country_par i.t1d_par   , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_age_edu_bmi_smok19822010",replace 
local Exposures="Maternal age at delivery"
local HR_CI=string(a[1,1],"%9.2f")+" ("+string(a[5,1],"%9.2f")+", "+string(a[6,1],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Upper-secondary or high school"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Pre-secondary or lower"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Maternal BMI during pregnancy"
local HR_CI=string(a[1,10],"%9.2f")+" ("+string(a[5,10],"%9.2f")+", "+string(a[6,10],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Maternal smoking during pregnancy"
local HR_CI=string(a[1,14],"%9.2f")+" ("+string(a[5,14],"%9.2f")+", "+string(a[6,14],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_age_edu_bmi_smok19822010.dta",clear
export excel using mor_age_edu_bmi_smok19822010,firstrow(variables) replace

/*by sex (for maternal smoking during pregnancy when calulating mediation proportion by childhood adiposity)*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\edu_3g_par.csv",clear
sort lopnr_barn
tempfile edu_3g_par
save `edu_3g_par'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
sort lopnr_barn
keep if birth_year<2011
merge 1:1 lopnr_barn using `edu_3g_par',keepusing(edu_3g_mor edu_3g_far) keep(1 3) nogen

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   

stcox c.age_mor_sd b2.edu_3g_mor b2.edu_3g_far c.bmi_mor_sd i.imputed_bmi b1.smoking1 i.birth_year i.country_par i.t1d_par  if men==1
stcox c.age_mor_sd b2.edu_3g_mor b2.edu_3g_far c.bmi_mor_sd i.imputed_bmi b1.smoking1 i.birth_year i.country_par i.t1d_par  if men==0

/*HR after adjusting for maternal education*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\edu_3g_par.csv",clear
sort lopnr_barn
tempfile edu_3g_par
save `edu_3g_par'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_trend.dta",clear
sort lopnr_barn
merge 1:1 lopnr_barn using `edu_3g_par',keepusing(edu_3g_mor) keep(1 3) nogen
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) enter(date_birth_barn) scale(365.25) failure(t1d18_any==1)   

capture postclose stats 
postfile stats str50 birth_year hr lci uci str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_230815\hr_per_year_adjedu",replace 
/*Adjusting for maternal education*/
stcox c.birth_year b2.edu_3g_mor if birth_year<2001,strata(men) vce(cluster lopnr_mor)
matrix b=r(table)
	local birth_year="Maternal eudcational level"
	local hr=b[1,1]
	local lci=b[5,1]
	local uci=b[6,1]
	local HR_CI=string(b[1,1],"%9.3f")+" ("+string(b[5,1],"%9.3f")+", "+string(b[6,1],"%9.3f")+")"
	post stats ("`birth_year'") (`hr') (`lci') (`uci') ("`HR_CI'") 

	postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_230815"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_230815\hr_per_year_adjedu.dta",clear
export excel using hr_per_year_adjedu,firstrow(variables) replace
