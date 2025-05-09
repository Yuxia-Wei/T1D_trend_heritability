/*Figure 4 & Supplementary Figure 4: Associations of environmental factors (infection during the first year of life) with T1D in the traditional cohort and sibling comparison analyses*/
/*infection_1st*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_1st_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
gen date_15month=mdy(month(date_birth_barn)+3,15,birth_year+1) if month(date_birth_barn)<=9
	replace date_15month=mdy(month(date_birth_barn)-9,15,birth_year+2) if month(date_birth_barn)>9
format date_15month %td

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) enter(time date_15month) scale(365.25) failure(t1d18_any==1)   

stcox b0.infection_1st i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_1st_19822010",replace 
local Exposures="Any infection"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_1st_19822010.dta",clear
export excel using infection_1st_19822010,firstrow(variables) replace

/*infection_type_1st*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_1st_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
gen date_15month=mdy(month(date_birth_barn)+3,15,birth_year+1) if month(date_birth_barn)<=9
	replace date_15month=mdy(month(date_birth_barn)-9,15,birth_year+2) if month(date_birth_barn)>9
format date_15month %td

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) enter(time date_15month) scale(365.25) failure(t1d18_any==1)   
stcox b0.infection_type_1st i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_1st_19822010",replace 
local Exposures="Virus infection"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Bacteria infection"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Other types of infection"
local HR_CI=string(a[1,4],"%9.2f")+" ("+string(a[5,4],"%9.2f")+", "+string(a[6,4],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_1st_19822010.dta",clear
export excel using infection_type_1st_19822010,firstrow(variables) replace

/*infection_site_1st*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_1st_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
gen date_15month=mdy(month(date_birth_barn)+3,15,birth_year+1) if month(date_birth_barn)<=9
	replace date_15month=mdy(month(date_birth_barn)-9,15,birth_year+2) if month(date_birth_barn)>9
format date_15month %td

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) enter(time date_15month) scale(365.25) failure(t1d18_any==1)   
stcox b0.infection_site_1st i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_1st_19822010",replace 
local Exposures="Respiratory infection"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Gastrointestinal infection"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Genitourinary infection"
local HR_CI=string(a[1,4],"%9.2f")+" ("+string(a[5,4],"%9.2f")+", "+string(a[6,4],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Other sites of infection"
local HR_CI=string(a[1,5],"%9.2f")+" ("+string(a[5,5],"%9.2f")+", "+string(a[6,5],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_1st_19822010.dta",clear
export excel using infection_site_1st_19822010,firstrow(variables) replace
