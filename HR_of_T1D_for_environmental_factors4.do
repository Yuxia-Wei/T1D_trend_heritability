/*Figure 4 & Supplementary Figure 5: Associations of environmental factors (serious life events) with T1D in the traditional cohort and sibling comparison analyses*/
/*Since age 0*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\timvar_dor_0yr_all.csv",clear
sort lopnr_barn
local outcome "start_timvar_dor_0yr stop_timvar_dor_0yr" 
foreach x of local outcome {
	gen `x'_new=date(`x',"YMD")
	format `x'_new %td
	drop `x'
	rename `x'_new `x'
}
codebook(timvar_dor_0yr t1d18_any_split start_timvar_dor_0yr stop_timvar_dor_0yr)

tempfile timvar_dor_0yr_all
save `timvar_dor_0yr_all'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\t1d_fulsib_covariate.csv",clear
sort lopnr_barn
tempfile t1d_fulsib_covariate
save `t1d_fulsib_covariate'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
sort lopnr_barn 
merge 1:1 lopnr_barn using `t1d_fulsib_covariate',keepusing(t1d_main_fulsib) keep(1 3) nogen
merge 1:m lopnr_barn using `timvar_dor_0yr_all',keepusing(timvar_dor_0yr t1d18_any_split start_timvar_dor_0yr stop_timvar_dor_0yr) keep(1 3) nogen
/*2,860,108  matched; 1,953,430 from master*/
tab t1d_main_fulsib 
tab timvar_dor_0yr
tab t1d18_any_split,missing

codebook( start_timvar_dor_0yr stop_timvar_dor_0yr)
replace t1d_main_fulsib=0 if t1d_main_fulsib==.
replace timvar_dor_0yr=0 if timvar_dor_0yr==.
replace start_timvar_dor_0yr=date_birth_barn if start_timvar_dor_0yr==.
replace stop_timvar_dor_0yr=min(date_t1d18_any,mdy(4,30,2020)) if stop_timvar_dor_0yr==.
/*it seems that "mdy(4,30,2020)" is not needed here*/

tab t1d18_any_split t1d18_any,missing
replace t1d18_any_split=t1d18_any if t1d18_any_split==.
tab timvar_dor_0yr t1d18_any_split

keep if birth_year<2011
stset stop_timvar_dor_0yr, id(lopnr_barn) origin(time date_birth_barn) enter(time start_timvar_dor_0yr) scale(365.25) failure(t1d18_any_split==1)   
/*The above codes should be rerun before every Cox model*/

capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr_cohort19822010.log",replace
stcox b0.timvar_dor_0yr i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr19822010.dta",clear
export excel using timvar_dor_0yr19822010,firstrow(variables) replace
/*adjusting t1d_main_fulsib*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr_adj_t1dfulsib_cohort19822010.log",replace
stcox b0.timvar_dor_0yr i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 i.t1d_main_fulsib, vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr_adj_t1dfulsib19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_0yr_adj_t1dfulsib19822010.dta",clear
export excel using timvar_dor_0yr_adj_t1dfulsib19822010,firstrow(variables) replace

/*Since age 5*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\timvar_dor_5yr_all.csv",clear
sort lopnr_barn
local outcome "start_timvar_dor_5yr stop_timvar_dor_5yr" 
foreach x of local outcome {
	gen `x'_new=date(`x',"YMD")
	format `x'_new %td
	drop `x'
	rename `x'_new `x'
}
codebook(timvar_dor_5yr t1d18_any_split start_timvar_dor_5yr stop_timvar_dor_5yr)

tempfile timvar_dor_5yr_all
save `timvar_dor_5yr_all'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\t1d_fulsib_covariate.csv",clear
sort lopnr_barn
tempfile t1d_fulsib_covariate
save `t1d_fulsib_covariate'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
sort lopnr_barn 
merge 1:1 lopnr_barn using `t1d_fulsib_covariate',keepusing(t1d_main_fulsib) keep(1 3) nogen
/*30207 matched*/
merge 1:m lopnr_barn using `timvar_dor_5yr_all',keepusing(timvar_dor_5yr t1d18_any_split start_timvar_dor_5yr stop_timvar_dor_5yr) keep(1 3) nogen
/*2,309,424  matched; 2,224,229 from master*/
gen date_5=mdy(month(date_birth_barn),15,birth_year+5) 
format date_5 %td

codebook(timvar_dor_5yr t1d18_any_split start_timvar_dor_5yr stop_timvar_dor_5yr)
replace t1d_main_fulsib=0 if t1d_main_fulsib==.
replace timvar_dor_5yr=0 if timvar_dor_5yr==.
replace start_timvar_dor_5yr=date_5 if start_timvar_dor_5yr==.
replace stop_timvar_dor_5yr=min(date_t1d18_any,mdy(4,30,2020)) if stop_timvar_dor_5yr==.

tab t1d18_any_split t1d18_any,missing
replace t1d18_any_split=t1d18_any if t1d18_any_split==.
tab timvar_dor_5yr t1d18_any_split
count if start_timvar_dor_5yr>=stop_timvar_dor_5yr
/*46,929*/
drop if start_timvar_dor_5yr>=stop_timvar_dor_5yr
keep if birth_year<2011

stset stop_timvar_dor_5yr, id(lopnr_barn) origin(time date_birth_barn) enter(time start_timvar_dor_5yr) scale(365.25) failure(t1d18_any_split==1)   

capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr_cohort19822010.log",replace
stcox b0.timvar_dor_5yr i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr19822010.dta",clear
export excel using timvar_dor_5yr19822010,firstrow(variables) replace

/*adjusting t1d_main_fulsib*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr_adj_t1dfulsib_cohort19822010.log",replace
stcox b0.timvar_dor_5yr i.t1d_main_fulsib i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr_adj_t1dfulsib19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_5yr_adj_t1dfulsib19822010.dta",clear
export excel using timvar_dor_5yr_adj_t1dfulsib19822010,firstrow(variables) replace

/*Since age 10*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\timvar_dor_10yr_all.csv",clear
sort lopnr_barn
local outcome "start_timvar_dor_10yr stop_timvar_dor_10yr" 
foreach x of local outcome {
	gen `x'_new=date(`x',"YMD")
	format `x'_new %td
	drop `x'
	rename `x'_new `x'
}
codebook(timvar_dor_10yr t1d18_any_split start_timvar_dor_10yr stop_timvar_dor_10yr)

tempfile timvar_dor_10yr_all
save `timvar_dor_10yr_all'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\t1d_fulsib_covariate.csv",clear
sort lopnr_barn
tempfile t1d_fulsib_covariate
save `t1d_fulsib_covariate'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
sort lopnr_barn 
merge 1:1 lopnr_barn using `t1d_fulsib_covariate',keepusing(t1d_main_fulsib) keep(1 3) nogen
merge 1:m lopnr_barn using `timvar_dor_10yr_all',keepusing(timvar_dor_10yr t1d18_any_split start_timvar_dor_10yr stop_timvar_dor_10yr) keep(1 3) nogen
/*1,600,497  matched; 2,573,699 from master*/
gen date_10=mdy(month(date_birth_barn),15,birth_year+10) 
format date_10 %td

tab t1d_main_fulsib 
codebook(timvar_dor_10yr t1d18_any_split start_timvar_dor_10yr stop_timvar_dor_10yr)
replace t1d_main_fulsib=0 if t1d_main_fulsib==.
replace timvar_dor_10yr=0 if timvar_dor_10yr==.
replace start_timvar_dor_10yr=date_10 if start_timvar_dor_10yr==.
replace stop_timvar_dor_10yr=min(date_t1d18_any,mdy(4,30,2020)) if stop_timvar_dor_10yr==.

tab t1d18_any_split t1d18,missing
replace t1d18_any_split=t1d18 if t1d18_any_split==.
tab timvar_dor_10yr t1d18_any_split
count if start_timvar_dor_10yr>=stop_timvar_dor_10yr
/*508674*/
drop if start_timvar_dor_10yr>=stop_timvar_dor_10yr
keep if birth_year<2011

stset stop_timvar_dor_10yr, id(lopnr_barn) origin(time date_birth_barn) enter(time start_timvar_dor_10yr) scale(365.25) failure(t1d18_any_split==1)   
/*The above codes should be run before each Cox model*/

capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr_cohort19822010.log",replace
stcox b0.timvar_dor_10yr i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr19822010.dta",clear
export excel using timvar_dor_10yr19822010,firstrow(variables) replace

/*adjusting t1d_main_fulsib*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr_adj_t1dfulsib_cohort19822010.log",replace
stcox b0.timvar_dor_10yr i.t1d_main_fulsib i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr_adj_t1dfulsib19822010",replace 
local Exposures="Death in second-degree relatives"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Death in first-degree relatives"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\timvar_dor_10yr_adj_t1dfulsib19822010.dta",clear
export excel using timvar_dor_10yr_adj_t1dfulsib19822010,firstrow(variables) replace
