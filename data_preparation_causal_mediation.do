/*create the dataset for causal mediation analysis*/
/*cumulative incidence for all birth years*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_trend_new.csv",clear
local outcome "birth_barn t1d18_any"
foreach x of local outcome {
  gen date_`x'_new=date(date_`x',"YMD")
    format date_`x'_new %td
    drop date_`x'
    rename date_`x'_new date_`x'
}

gen age_stop=(date_t1d18_any-date_birth_barn)/365.25
codebook( age_stop) if t1d18_any==0  /*until age 18.997946 years!*/
codebook( age_stop) if t1d18_any==1  /*Until age 18.997946 years*/

keep birth_year age_stop t1d18_any
export delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_cumulative_inc.csv",nolabel replace

/*stata*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\edu_3g_par.csv",clear
sort lopnr_barn
tempfile edu_3g_par
save `edu_3g_par'
/*'*/
import delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_trend_new.csv",clear
local outcome "birth_barn t1d18_any"
foreach x of local outcome {
  gen date_`x'_new=date(date_`x',"YMD")
    format date_`x'_new %td
    drop date_`x'
    rename date_`x'_new date_`x'
}

sort lopnr_barn
merge 1:1 lopnr_barn using `edu_3g_par',keepusing(edu_3g_mor) keep(1 3) nogen
merge 1:1 lopnr_barn using "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_perinatal_cohort.dta",keepusing(smoking1 age_mor ga_5g bw_5g infection_type infection_site caesarean1) keep(1 3) nogen
/*'*/
gen age_stop=(date_t1d18_any-date_birth_barn)/365.25
gen year2000_3g=1 if birth_year==1982
replace year2000_3g=2 if birth_year==1983
replace year2000_3g=3 if birth_year==2000

gen year2000=1 if birth_year==2000
	replace year2000=0 if birth_year==1982

gen year2000_smok=1 if birth_year==2000
	replace year2000_smok=0 if birth_year==1983

/*calculate sample size in 192-2000 without missing*/
count if age_mor!=. & birth_year>=1982 & birth_year<=2000
count if edu_3g_mor!=9 & birth_year>=1982 & birth_year<=2000

count if smoking1!=9 & birth_year!=1982
count if smoking1!=9 & birth_year>1982 & birth_year<=2000

count if infection_type!=9 & birth_year>=1982 & birth_year<=2000
count if caesarean1!=9 & birth_year>=1982 & birth_year<=2000
count if ga_5g!=9 & birth_year>=1982 & birth_year<=2000
count if bw_5g!=9 & birth_year>=1982 & birth_year<=2000

codebook(year2000)
tab smoking1,missing
tab edu_3g_mor,missing
hist age_mor
drop if year2000_3g==. /*125487 individuals remained*/
keep smoking1 age_mor edu_3g_mor infection_site infection_type ga_5g bw_5g caesarean1 year2000_3g year2000 men age_stop t1d18_any
tab t1d18_any if year2000==1  /*0.88%*/
tab t1d18_any if year2000==0 /*0.38%*/
export delimited "W:\C6_Carlsson\Yuxia Wei\data\familial\mbr_mediation.csv",nolabel replace
