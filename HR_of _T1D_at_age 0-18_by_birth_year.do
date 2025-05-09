/*Figure 1(A):Hazard ratio (95% CI) of T1D at age 0-18 years in each birth year (birth year for reference: 1982)*/
/*basic model*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_trend_new.csv",clear

local outcome "birth_barn t1d18_any"
foreach x of local outcome {
	gen date_`x'_new=date(date_`x',"YMD")
    format date_`x'_new %td
    drop date_`x'
    rename date_`x'_new date_`x'
}

capture cap log close
log using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\HR_birth_year_t1d18_any.log",replace
stset date_t1d18_any , id(lopnr_barn) origin(time date_birth_barn) enter(date_birth_barn) scale(365.25) failure(t1d18_any==1)  

stcox i.birth_year ,strata(men) vce(cluster lopnr_mor)
matrix b=r(table)

capture postclose stats 
postfile stats birth_year hr lci uci using "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\hr_birth_year_t1d18_any",replace 
foreach i of numlist 1/29 {
    local birth_year=`i'+1981
    local hr=b[1,`i']
    local lci=b[5,`i']
    local uci=b[6,`i']
    post stats (`birth_year') (`hr') (`lci') (`uci') 
}
postclose stats
cd "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418"
use "W:\C6_Carlsson\Yuxia Wei\results\familial\trend_240418\hr_birth_year_t1d18_any.dta",clear
export excel using hr_birth_year_t1d18_any,firstrow(variables) replace
