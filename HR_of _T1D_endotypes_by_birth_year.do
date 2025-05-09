/*Supplementary Figure 1:Hazard ratio (95% CI) of T1D at age 0-6, 6-12, and 13-18 years in each birth year (birth year for reference: 1982)*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_trend_new.csv",clear

local outcome "birth_barn t1d18_any"
foreach x of local outcome {
	gen date_`x'_new=date(date_`x',"YMD")
    format date_`x'_new %td
    drop date_`x'
    rename date_`x'_new date_`x'

}
codebook 
gen date_age7=mdy(month(date_birth_barn),15,birth_year+7) 
format date_age7 %td  /*the date before age 7 years*/

gen t1d6=1 if t1d18_any==1 & date_t1d18_any<date_age7
	replace t1d6=0 if t1d18_any!=1 | (t1d18_any==1 & date_t1d18_any>=date_age7)

gen date_t1d6=date_t1d18_any if t1d6==1
	replace date_t1d6=min(date_t1d18_any,date_age7) if t1d6==0
format date_t1d6 %td


gen date_age13=mdy(month(date_birth_barn),15,birth_year+13) 
format date_age13 %td 

gen t1d712=1 if t1d18_any==1 & (date_t1d18_any>=date_age7 & date_t1d18_any<date_age13)
	replace t1d712=0 if date_t1d18_any>=date_age7 & (t1d18_any!=1 | (t1d18_any==1 & date_t1d18_any>=date_age13))
gen date_t1d712=date_t1d18_any if t1d712==1
	replace date_t1d712=min(date_t1d18_any,date_age13) if t1d712==0 
	replace date_t1d712=date_age7+1 if t1d712==1 & date_t1d712==date_age7
format date_t1d712 %td

gen t1d1318=1 if t1d18_any==1 & date_t1d18_any>=date_age13
	replace t1d1318=0 if date_t1d18_any>=date_age13 & t1d18_any==0
gen date_t1d1318=date_t1d18_any if t1d1318==1 | t1d1318==0
format date_t1d1318 %td

/*data analysis*/
/*0-6 years*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d6_HR_birth_year.log",replace
stset date_t1d6 , id(lopnr_barn) origin(time date_birth_barn) enter(date_birth_barn) scale(365.25) failure(t1d6==1)  

stcox i.birth_year ,strata(men) vce(cluster lopnr_mor)
matrix b=r(table)

capture postclose stats 
postfile stats birth_year hr lci uci using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d6_hr_birth_year",replace 
foreach i of numlist 1/29 {
	local birth_year=`i'+1981
	local hr=b[1,`i']
	local lci=b[5,`i']
	local uci=b[6,`i']
	post stats (`birth_year') (`hr') (`lci') (`uci') 
}
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d6_hr_birth_year.dta",clear
gen Age="0-6 years"
save "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d6_hr_birth_year.dta",replace
/*export excel using t1d6_hr_birth_year,firstrow(variables) replace*/

/*7-12 years*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d712_HR_birth_year.log",replace
stset date_t1d712 , id(lopnr_barn) origin(time date_birth_barn) enter(date_age7) scale(365.25) failure(t1d712==1)  

stcox i.birth_year ,strata(men) vce(cluster lopnr_mor)
matrix b=r(table)

capture postclose stats 
postfile stats birth_year hr lci uci using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d712_hr_birth_year",replace 
foreach i of numlist 1/29 {
	local birth_year=`i'+1981
	local hr=b[1,`i']
	local lci=b[5,`i']
	local uci=b[6,`i']
	post stats (`birth_year') (`hr') (`lci') (`uci') 
}
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d712_hr_birth_year.dta",clear
replace birth_year=birth_year+0.2
gen Age="7-12 years"
save "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d712_hr_birth_year.dta",replace

/*export excel using t1d712_hr_birth_year,firstrow(variables) replace*/

/*age 13-18 years*/
capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d1318_HR_birth_year.log",replace
stset date_t1d1318 , id(lopnr_barn) origin(time date_birth_barn) enter(date_age13) scale(365.25) failure(t1d1318==1)  

stcox i.birth_year ,strata(men) vce(cluster lopnr_mor)
matrix b=r(table)

capture postclose stats 
postfile stats birth_year hr lci uci using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d1318_hr_birth_year",replace 
foreach i of numlist 1/29 {
	local birth_year=`i'+1981
	local hr=b[1,`i']
	local lci=b[5,`i']
	local uci=b[6,`i']
	post stats (`birth_year') (`hr') (`lci') (`uci') 
}
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d1318_hr_birth_year.dta",clear
replace birth_year=birth_year+0.4
gen Age="13-18 years"
save "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d1318_hr_birth_year.dta",replace

/*export excel using t1d1318_hr_birth_year,firstrow(variables) replace*/

use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d6_hr_birth_year.dta",clear
append using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d712_hr_birth_year.dta"
append using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\t1d1318_hr_birth_year.dta"
sort birth_year
export excel using Figure_HR_birth_year_finer_age,firstrow(variables) replace
