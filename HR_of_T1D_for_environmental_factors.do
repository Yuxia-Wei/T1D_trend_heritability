/*Figure 4& Supplementary Figure 4: Associations of environmental factors with T1D in the traditional cohort and sibling comparison analyses*/
/****************************data analysis*********************************/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_single_birth_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox i.mor_single_birth i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
/*no association*/
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_single_birth_19822010",replace 
local Exposures="Single mother at childbirth"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mor_single_birth_19822010.dta",clear
export excel using mor_single_birth_19822010,firstrow(variables) replace

/*order_3g*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\order_3g_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox i.order_3g i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\order_3g_19822010",replace 
local Exposures="Second-born"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Third or higher"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\order_3g_19822010.dta",clear
export excel using order_3g_19822010,firstrow(variables) replace

/*infection_preg*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_preg_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   

stcox b0.infection_preg i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_preg_19822010",replace 
local Exposures="Any infection"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_preg_19822010.dta",clear
export excel using infection_preg_19822010,firstrow(variables) replace

/*infection_type*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox b0.infection_type i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_19822010",replace 
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
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_type_19822010.dta",clear
export excel using infection_type_19822010,firstrow(variables) replace

/*infection_site*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox b0.infection_site i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_19822010",replace 
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
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\infection_site_19822010.dta",clear
export excel using infection_site_19822010,firstrow(variables) replace

/*mood_preg*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mood_preg_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   

stcox b0.mood_preg i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mood_preg_19822010",replace 
local Exposures="Maternal mood disorders during pregnancy"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\mood_preg_19822010.dta",clear
export excel using mood_preg_19822010,firstrow(variables) replace

/*mode of delivery*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240123\caesarean1_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   

stcox b0.caesarean1 i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240123\caesarean1_19822010",replace 
local Exposures="Maternal mood disorders during pregnancy"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240123"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240123\caesarean1_19822010.dta",clear
export excel using caesarean1_19822010,firstrow(variables) replace

/*gestational age*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\ga_5g_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox b4.ga_5g i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\ga_5g_19822010",replace 
local Exposures="Very preterm"
local HR_CI=string(a[1,1],"%9.2f")+" ("+string(a[5,1],"%9.2f")+", "+string(a[6,1],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Preterm"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Early term"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Postterm"
local HR_CI=string(a[1,5],"%9.2f")+" ("+string(a[5,5],"%9.2f")+", "+string(a[6,5],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\ga_5g_19822010.dta",clear
export excel using ga_5g_19822010,firstrow(variables) replace

/*birth weight*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\bw_5g_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox b4.bw_5g i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\bw_5g_19822010",replace 
local Exposures="Birth weight<1500"
local HR_CI=string(a[1,1],"%9.2f")+" ("+string(a[5,1],"%9.2f")+", "+string(a[6,1],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight 1500-2499g"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight 2500-2999g"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight>=4000g"
local HR_CI=string(a[1,5],"%9.2f")+" ("+string(a[5,5],"%9.2f")+", "+string(a[6,5],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\bw_5g_19822010.dta",clear
export excel using bw_5g_19822010,firstrow(variables) replace

/*birth weight for gestational age*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\sga_new_cohort19822010.log",replace
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",clear
keep if birth_year<2011
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
stcox b0.sga_new i.men i.birth_year i.country_par c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.t1d_par i.smoking1 , vce(cluster lopnr_mor)  
matrix a=r(table)

capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\sga_new_19822010",replace 
local Exposures="Small for gestational age"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Large for gestational age"
local HR_CI=string(a[1,3],"%9.2f")+" ("+string(a[5,3],"%9.2f")+", "+string(a[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\sga_new_19822010.dta",clear
export excel using sga_new_19822010,firstrow(variables) replace

/*sibling comparison in children born in 1982-2010*/
use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta", clear
keep if birth_year<2011
duplicates tag lopnr_mor, generate(withsib)
drop if withsib==0
codebook(lopnr_mor)  
/*the unique number of lopnr_mor: the number of sibling groups: 1,145,631*/
codebook(lopnr_mor) if t1d18_any==1  /*15,697 */
codebook(lopnr_mor) if t1d18_any==0  /*1,145,301 */

tempfile t1d18_any_family
keep if t1d18_any==1
gen t1d18_any_family=1
keep lopnr_mor t1d18_any_family
sort lopnr_mor 
duplicates drop 
save `t1d18_any_family'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\edu_3g_par.csv",clear /*including those meeting exclusion criteria*/
sort lopnr_barn
tempfile edu_3g_par
save `edu_3g_par'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta", clear
keep if birth_year<2011
duplicates tag lopnr_mor, generate(withsib)
drop if withsib==0
tempfile control_family
keep if t1d18_any==0
gen control_family=1
keep lopnr_mor control_family
sort lopnr_mor 
duplicates drop 
save `control_family'

use "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta", clear
sort lopnr_mor 
merge 1:1 lopnr_barn using `edu_3g_par',keepusing(edu_3g_mor edu_3g_far) keep(1 3) nogen
merge m:1 lopnr_mor using `t1d18_any_family',keep(1 3) keepusing(t1d18_any_family) nogen
merge m:1 lopnr_mor using `control_family',keep(1 3) keepusing(control_family) nogen
codebook(t1d18_any_family control_family)
keep if t1d18_any_family==1   /*41618 individuals from 15795 sibling groups discordant on t1d18 remained for sibling analysis*/
save "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_sib.dta",replace
codebook(lopnr_mor)

/*maternal age at delivery in the model turned the HR in relation to birth year from >1 to <1!*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\HR_perinatal_sib19822010.log",replace


capture postclose stats 
postfile stats str50 Exposures str50 HR_CI using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\HR_perinatal_sib19822010",replace 
stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) scale(365.25) failure(t1d18_any==1)   
/*Single mothers at childbirth*/
stcox i.mor_single_birth i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix d=r(table)

local Exposures="Single mothers at childbirth"
local HR_CI=string(d[1,2],"%9.2f")+" ("+string(d[5,2],"%9.2f")+", "+string(d[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

/*maternal BMI and smoking during pregnancy*/
stcox i.smoking1 c.bmi_mor_sd i.imputed_bmi  i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor , strata(lopnr_mor)  
matrix a=r(table)

local Exposures="Maternal BMI during pregnancy"
local HR_CI=string(a[1,4],"%9.2f")+" ("+string(a[5,4],"%9.2f")+", "+string(a[6,4],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 


local Exposures="Maternal smoking during pregnancy"
local HR_CI=string(a[1,2],"%9.2f")+" ("+string(a[5,2],"%9.2f")+", "+string(a[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

/*Maternal infection during pregnancy*/
stcox b0.infection_preg i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  

stcox b0.infection_type i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix b=r(table)

stcox b0.infection_site i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix c=r(table)

local Exposures="Maternal Virus infection"
local HR_CI=string(b[1,2],"%9.2f")+" ("+string(b[5,2],"%9.2f")+", "+string(b[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Maternal Bacteria infection"
local HR_CI=string(b[1,3],"%9.2f")+" ("+string(b[5,3],"%9.2f")+", "+string(b[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 


local Exposures="Maternal Respiratory infection"
local HR_CI=string(c[1,2],"%9.2f")+" ("+string(c[5,2],"%9.2f")+", "+string(c[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Maternal Gastrointestinal infection"
local HR_CI=string(c[1,3],"%9.2f")+" ("+string(c[5,3],"%9.2f")+", "+string(c[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Maternal Genitourinary infection"
local HR_CI=string(c[1,4],"%9.2f")+" ("+string(c[5,4],"%9.2f")+", "+string(c[6,4],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
/*Maternal mood disorders during pregnancy*/
stcox b0.mood_preg i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix d=r(table)

local Exposures="Maternal mood disorders during pregnancy"
local HR_CI=string(d[1,2],"%9.2f")+" ("+string(d[5,2],"%9.2f")+", "+string(d[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

stcox b0.caesarean1 i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix d=r(table)

local Exposures="Cesarean"
local HR_CI=string(d[1,2],"%9.2f")+" ("+string(d[5,2],"%9.2f")+", "+string(d[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
/*gestational age*/
stcox b4.ga_5g i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix f=r(table)

local Exposures="Very preterm"
local HR_CI=string(f[1,1],"%9.2f")+" ("+string(f[5,1],"%9.2f")+", "+string(f[6,1],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Preterm"
local HR_CI=string(f[1,2],"%9.2f")+" ("+string(f[5,2],"%9.2f")+", "+string(f[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Early term"
local HR_CI=string(f[1,3],"%9.2f")+" ("+string(f[5,3],"%9.2f")+", "+string(f[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Postterm"
local HR_CI=string(f[1,5],"%9.2f")+" ("+string(f[5,5],"%9.2f")+", "+string(f[6,5],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

/*small for gestational age*/
stcox b0.sga_new i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix g=r(table)
local Exposures="Small for gestational age"
local HR_CI=string(g[1,2],"%9.2f")+" ("+string(g[5,2],"%9.2f")+", "+string(g[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'")

local Exposures="Large for gestational age"
local HR_CI=string(g[1,3],"%9.2f")+" ("+string(g[5,3],"%9.2f")+", "+string(g[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'")
/*birth weight<1500g*/
stcox b4.bw_5g i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix h=r(table)
local Exposures="Birth weight<1500g"
local HR_CI=string(h[1,1],"%9.2f")+" ("+string(h[5,1],"%9.2f")+", "+string(h[6,1],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight 1500-2499g"
local HR_CI=string(h[1,2],"%9.2f")+" ("+string(h[5,2],"%9.2f")+", "+string(h[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight 2500-2999g"
local HR_CI=string(h[1,3],"%9.2f")+" ("+string(h[5,3],"%9.2f")+", "+string(h[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Birth weight 4000 or above"
local HR_CI=string(h[1,5],"%9.2f")+" ("+string(h[5,5],"%9.2f")+", "+string(h[6,5],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
/*Infection during the first year of life: the follow-up should be started from 15 months after birth*/
gen date_15month=mdy(month(date_birth_barn)+3,15,birth_year+1) if month(date_birth_barn)<=9
	replace date_15month=mdy(month(date_birth_barn)-9,15,birth_year+2) if month(date_birth_barn)>9
format date_15month %td

stset date_t1d18_any, id(lopnr_barn) origin(time date_birth_barn) enter(time date_15month) scale(365.25) failure(t1d18_any==1)   

stcox b0.infection_1st i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  

stcox b0.infection_type_1st i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix i=r(table)

stcox b0.infection_site_1st i.men i.birth_year c.age_mor_sd i.edu_post2nd_far i.edu_post2nd_mor c.bmi_mor_sd i.imputed_bmi i.smoking1 , strata(lopnr_mor)  
matrix j=r(table)

local Exposures="Virus infection"
local HR_CI=string(i[1,2],"%9.2f")+" ("+string(i[5,2],"%9.2f")+", "+string(i[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Bacteria infection"
local HR_CI=string(i[1,3],"%9.2f")+" ("+string(i[5,3],"%9.2f")+", "+string(i[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Respiratory infection"
local HR_CI=string(j[1,2],"%9.2f")+" ("+string(j[5,2],"%9.2f")+", "+string(j[6,2],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Gastrointestinal infection"
local HR_CI=string(j[1,3],"%9.2f")+" ("+string(j[5,3],"%9.2f")+", "+string(j[6,3],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 

local Exposures="Genitourinary infection"
local HR_CI=string(j[1,4],"%9.2f")+" ("+string(j[5,4],"%9.2f")+", "+string(j[6,4],"%9.2f")+")"
post stats ("`Exposures'") ("`HR_CI'") 
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\HR_perinatal_sib19822010.dta",clear
export excel using HR_perinatal_sib19822010,firstrow(variables) replace
