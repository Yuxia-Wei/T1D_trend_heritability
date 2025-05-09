/*import and merge datasets*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\country_infection_t1dpar_export.csv",clear /*for those born between 1982 and 2014, including those meeting exclusion criteria*/
sort lopnr_barn
tab infection_preg infection_type,missing
tab infection_preg infection_site,missing
tab infection_preg_sv infection_type_sv,missing
tab infection_preg_sv infection_site_sv,missing
tab infection_preg_ov infection_type_ov,missing
tab infection_preg_ov infection_site_ov,missing
tempfile country_infection_t1dpar_export
save `country_infection_t1dpar_export'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mor_infection_priorpregnancy.csv",clear /*birth_year 1973-2014; including those meeting exclusion criteria*/
sort lopnr_barn  /*unique lopnr_barn*/
tab infection_prior infection_type_prior,missing
tab infection_prior infection_site_prior,missing
tab infection_prior_sv infection_type_prior_sv,missing
tab infection_prior_sv infection_site_prior_sv,missing
tab infection_prior_ov infection_type_prior_ov,missing
tab infection_prior_ov infection_site_prior_ov,missing
tempfile mor_infection_priorpregnancy
save `mor_infection_priorpregnancy'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\child_infection_1st.csv",clear /*for children born in 1973-2014, including those meeting exclusion criteri*/
sort lopnr_barn /*unique lopnr_barn*/
tab infection_1st infection_type_1st,missing
tab infection_1st infection_site_1st,missing
tab infection_1st_sv infection_type_1st_sv,missing
tab infection_1st_sv infection_site_1st_sv,missing
tab infection_1st_ov infection_type_1st_ov,missing
tab infection_1st_ov infection_site_1st_ov,missing
tempfile child_infection_1st
save `child_infection_1st'


import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mor_mood_pregnancy.csv",clear /*for children born in 1973-2014, including those meeting exclusion criteria*/
sort lopnr_barn /*unique lopnr_barn*/
tab mood_preg ,missing
tab mood_preg_sv ,missing
tab mood_preg_ov ,missing
tempfile mor_mood_pregnancy
save `mor_mood_pregnancy'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\famsit_birth.csv",clear 
codebook
gen mor_single_birth=1 if famsit!=. & famsit!=1
replace mor_single_birth=0 if famsit==1
replace mor_single_birth=9 if famsit==.
tab mor_single_birth
sort lopnr_barn /*unique lopnr_barn*/
tempfile famsit_birth
save `famsit_birth'

import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.csv",clear /*this dataset has been updated according to t1d18_any*/
local outcome "birth_barn t1d18_any" 
foreach x of local outcome {
	gen date_`x'_new=date(date_`x',"YMD")
	format date_`x'_new %td
	drop date_`x'
	rename date_`x'_new date_`x'
}
drop t1d_par

sort lopnr_barn /*unique lopnr_barn*/
merge 1:1 lopnr_barn using `country_infection_t1dpar_export',keepusing (country_par t1d_par infection*) keep(1 3) nogen
/*3,3471,310 matched; 0 not matched*/
merge 1:1 lopnr_barn using `mor_infection_priorpregnancy',keepusing (infection*) keep(1 3) nogen
/* 141,160 matched*/
merge 1:1 lopnr_barn using `child_infection_1st',keepusing (infection*) keep(1 3) nogen
/*648,218 matched*/
merge 1:1 lopnr_barn using `mor_mood_pregnancy',keepusing (mood*) keep(1 3) nogen
/*16715*/
tab t1d_par,missing  
merge 1:1 lopnr_barn using `famsit_birth',keepusing(mor_single_birth) keep(1 3) nogen
/*3,371,310  matched*/
tab mor_single_birth


tab infection_preg,missing
tab infection_type,missing
tab infection_site,missing

tab infection_preg_sv,missing
tab infection_type_sv,missing
tab infection_site_sv,missing

tab infection_preg_ov,missing
tab infection_type_ov,missing
tab infection_site_ov,missing

tab infection_prior,missing
tab infection_type_prior,missing
tab infection_site_prior,missing

tab infection_prior_sv,missing
tab infection_type_prior_sv,missing
tab infection_site_prior_sv,missing

tab infection_prior_ov,missing
tab infection_type_prior_ov,missing
tab infection_site_prior_ov,missing

tab infection_1st,missing
tab infection_type_1st,missing
tab infection_site_1st,missing

tab infection_1st_sv,missing
tab infection_type_1st_sv,missing
tab infection_site_1st_sv,missing

tab infection_1st_ov,missing
tab infection_type_1st_ov,missing
tab infection_site_1st_ov,missing

tab mood_preg,missing
tab mood_preg_sv,missing
tab mood_preg_ov,missing
/*0 missing: 0.45% of parent T1D; 0.46% of potential parental T1D*/

replace infection_preg=0 if infection_preg==.
replace infection_type=0 if infection_type==.
replace infection_site=0 if infection_site==.

replace infection_preg_sv=0 if infection_preg_sv==.
replace infection_type_sv=0 if infection_type_sv==.
replace infection_site_sv=0 if infection_site_sv==.

replace infection_preg_ov=0 if infection_preg_ov==.
replace infection_type_ov=0 if infection_type_ov==.
replace infection_site_ov=0 if infection_site_ov==.

replace infection_prior=0 if infection_prior==.
replace infection_type_prior=0 if infection_type_prior==.
replace infection_site_prior=0 if infection_site_prior==.

replace infection_prior_sv=0 if infection_prior_sv==.
replace infection_type_prior_sv=0 if infection_type_prior_sv==.
replace infection_site_prior_sv=0 if infection_site_prior_sv==.

replace infection_prior_ov=0 if infection_prior_ov==.
replace infection_type_prior_ov=0 if infection_type_prior_ov==.
replace infection_site_prior_ov=0 if infection_site_prior_ov==.

replace infection_1st=0 if infection_1st==.
replace infection_type_1st=0 if infection_type_1st==.
replace infection_site_1st=0 if infection_site_1st==.

replace infection_1st_sv=0 if infection_1st_sv==.
replace infection_type_1st_sv=0 if infection_type_1st_sv==.
replace infection_site_1st_sv=0 if infection_site_1st_sv==.

replace infection_1st_ov=0 if infection_1st_ov==.
replace infection_type_1st_ov=0 if infection_type_1st_ov==.
replace infection_site_1st_ov=0 if infection_site_1st_ov==.

replace mood_preg=0 if mood_preg==.
replace mood_preg_sv=0 if mood_preg_sv==.
replace mood_preg_ov=0 if mood_preg_ov==.

codebook(date_t1d18_any) if t1d18_any==1 
codebook(date_t1d18_any) if t1d18_any==0 

save "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta", replace
