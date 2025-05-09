/*prepare sibling pairs used in the calculation of heritability at age 0-18 years*/
/*two datasets are created: (1) all possible full sibling pairs from each family; (2) one full sibling pair from each family*/
libname familial "W:\C6_Carlsson\Yuxia Wei\data\familial";
data mbr_t1d18_any;
set familial.mbr_merged;
   if contra_sex=1 or contra_lopnrmor=1 or mis_lopnrmor=1 or unreliable_lopnrfar=1 
        or stillbirth^="" or year(date_t1d18_any)<year(date_birth_barn)
        or year(date_t1d18_any)=year(date_birth_barn) and month(date_t1d18_any)<=month(date_birth_barn)
        or index^=1 then exc_trend_t1d_any=1;
where birth_year<2011 and birth_year>=1982;
run;
data mbr_t1d18_any1;
set mbr_t1d18_any;
where exc_trend_t1d_any^=1 and multi_birth="1";
run; /*2852600*/

proc freq data=mbr_t1d18_any1;
tables t1d18_any;
run; /*19642: not 20,086?*/
data a;
set mbr_t1d18_any1;
where lopnr_barn=. or lopnrindex=.;
run; /*0*/
/*index persons*/
data index_t1d(keep=lopnrindex lopnr_mor date_birth_index men_index t1d18_index);
set mbr_t1d18_any1;
	if t1d18_any=1 then t1d18_index=1;
	if t1d18_any=0 then t1d18_index=0;
rename lopnr_barn=lopnrindex;
rename date_birth_barn=date_birth_index;
rename men=men_index;
run; 
/*T1D cases diagnosed at different ages in the full cohort (n=2852600)*/

/*siblings*/
data lopnr_sibling(keep=lopnr date_birth t1d18 men);
set mbr_t1d18_any1;
	if t1d18_any=1  then t1d18=1;
	if t1d18_any=0  then t1d18=0;
rename lopnr_barn=lopnr;
rename date_birth_barn=date_birth;
run; /*2852600*/

data sibling_mgr(keep=lopnr lopnrindex relation);
set familial.mgr;
where relation=:"0010";
run; /*full siblings*/

proc sql;
create table sibling_merged
	as select a.*, b.lopnrindex,b.relation
	from lopnr_sibling as a left join  sibling_mgr as b
	on a.lopnr=b.lopnr;
quit; /*lopnr of siblings: lopnr*/

data sibling_merged1;
set sibling_merged;
where relation^="";
run;

proc sql;
create table sibling_match
	as select a.*, b.*
	from index_t1d as a left join  sibling_merged1 as b
	on a.lopnrindex=b.lopnrindex;
quit;

data familial.fulsib_t1d18_any_all;
set sibling_match;
where relation^="";
run; /*fewer than familial.fulsib_t1d_all included children born <1982 and >2010*/


/********************************one sibling pair from each mother and father (children with missing lopnr_far excluded)********************************************/
/*to get lopnr_far*/
data lopnr_far(keep=lopnrindex lopnr relation);
set familial.mgr_sibling_cousin;
where relation=:"1000";
run;
proc sort data=lopnr_far nodupkey;
by lopnrindex lopnr;
run; /*4241381*/
proc sort data=lopnr_far nodupkey;
by lopnrindex;
run; /*4241381*/


proc sql;
create table fulsib_t1d_all_lopnrfar
	as select a.*, b.lopnr as lopnr_far_index
	from familial.fulsib_t1d18_any_all as a left join  lopnr_far as b
	on a.lopnrindex=b.lopnrindex;
quit;

proc sql;
create table fulsib_t1d_all_lopnrfar1
	as select a.*, b.lopnr as lopnr_far
	from fulsib_t1d_all_lopnrfar as a left join  lopnr_far as b
	on a.lopnr=b.lopnrindex;
quit;
proc freq data=fulsib_t1d_all_lopnrfar1;
tables relation/missing;
where lopnr_far_index=. or lopnr_far=.;
run;

/*any birth cohort before 2010*/
data fulsib_t1d18_any(keep=lopnrindex lopnr_mor lopnr_far lopnr men_index men t1d18_index t1d18 date_birth_index date_birth relation);
set fulsib_t1d_all_lopnrfar1;
	where year(date_birth_index)<2011 and year(date_birth)<2011 and lopnr_far_index^=. and lopnr_far^=.;
run;

proc sort data=fulsib_t1d18_any;
by lopnr_mor date_birth_index date_birth;
run; /*3098736 pairs of siblings. For each mother, only one pair of sibligns should be selected*/
proc sort data=fulsib_t1d18_any out=fulsib_t1d18_any0 nodupkey;
by lopnr_mor;
run; /*898038 pairs of siblings were selected by lopnr_mor*/
proc sort data=fulsib_t1d18_any0;
by lopnr_far date_birth_index date_birth;
run; 
proc sort data=fulsib_t1d18_any0 out=fulsib_t1d18_any1 dupout=duplicated_lopnr nodupkey;
by lopnr_far;
run;  /*884975 pairs remained; 13063 pairs deleted*/

data fulsib_t1d18_any;
set fulsib_t1d18_any1;
	birth_year_index=year(date_birth_index);
	birth_year=year(date_birth);
run; 
/*to be exported as fulsib_t1d18_any.csv*/
proc freq data=fulsib_t1d18_any;
tables t1d18_index*t1d18;
run; /*457 concordant pairs, 5627+5860 discordant pairs */
/***********************all possible pairs of full siblings****************/
/*to get lopnr_far*/
data lopnr_far(keep=lopnrindex lopnr relation);
set familial.mgr_sibling_cousin;
where relation=:"1000";
run;
proc sort data=lopnr_far nodupkey;
by lopnrindex lopnr;
run; /*4241381*/
proc sort data=lopnr_far nodupkey;
by lopnrindex;
run; /*4241381*/

proc sql;
create table fulsib_t1d_all_lopnrfar
	as select a.*, b.lopnr as lopnr_far_index
	from familial.fulsib_t1d18_any_all as a left join  lopnr_far as b
	on a.lopnrindex=b.lopnrindex;
quit;

proc sql;
create table fulsib_t1d_all_lopnrfar1
	as select a.*, b.lopnr as lopnr_far
	from fulsib_t1d_all_lopnrfar as a left join  lopnr_far as b
	on a.lopnr=b.lopnrindex;
quit;
proc freq data=fulsib_t1d_all_lopnrfar1;
tables relation/missing;
where lopnr_far_index=. or lopnr_far=.;
run;

/*exclude children with missing lopnr_far*/
data fulsib_t1d18_lopnrfar2;
set fulsib_t1d_all_lopnrfar1;
where lopnr_far_index^=. and lopnr_far^=.;
run; /*one pair of siblings appear twice in this dataset*/

data fulsib_t1d18_any_allpair1;
set fulsib_t1d18_lopnrfar2;
	where date_birth_index<date_birth and 
		year(date_birth_index)<2011 and year(date_birth)<2011
		and year(date_birth_index)>1981 and year(date_birth)>1981;
run; /*1549362 for t1d18*/ 
data fulsib_t1d18_any_allpair2;
set fulsib_t1d18_lopnrfar2;
	where date_birth_index>date_birth and 
		year(date_birth_index)<2011 and year(date_birth)<2011
		and year(date_birth_index)>1981 and year(date_birth)>1981;
run; /*1549362 for t1d18*/

data fulsib_t1d18_any_allpair3;
set fulsib_t1d18_any_allpair2;
	rename lopnr=lopnr2;
	rename lopnr_far=lopnr_far2;
	rename date_birth=date_birth2;
	rename men=men2;
	rename t1d18=t1d182;
run; 

data fulsib_t1d18_any_allpair4;
set fulsib_t1d18_any_allpair3;
	rename lopnrindex=lopnr;
	rename lopnr_far_index=lopnr_far;
	rename men_index=men;
	rename date_birth_index=date_birth;
	rename t1d18_index=t1d18;
run;
data fulsib_t1d18_any_allpair5;
set fulsib_t1d18_any_allpair4;
	rename lopnr2=lopnrindex;
	rename lopnr_far2=lopnr_far_index;
	rename date_birth2=date_birth_index;
	rename men2=men_index;
	rename t1d182=t1d18_index;
run;
data fulsib_t1d18_any_allpair;
set fulsib_t1d18_any_allpair1 fulsib_t1d18_any_allpair5;
run;
data fulsib_t1d18_any_allpair1;
set fulsib_t1d18_any_allpair;
	birth_year_index=year(date_birth_index);
	birth_year=year(date_birth);
run;
proc sort data=fulsib_t1d18_any_allpair1;
by lopnr_mor date_birth_index date_birth;
run;
proc sort data=fulsib_t1d18_any_allpair1 nodupkey dupout=dup_exc;
by lopnr_mor date_birth_index date_birth;
run; /*1549357 pairs remained for t1d18(if twins are matched to different their siblings, only one person in the twin paris remained)*/
proc freq data=fulsib_t1d18_any_allpair1;
tables t1d18_index*t1d18;
where relation=:"0010";
run; /*841 concordant pairs, 9635+10382 discordant paris for t1d18*/


/*to get familyid*/
data lopnr_far_index_mor(keep=lopnr_far_index lopnr_mor);
set fulsib_t1d18_any_allpair1;
run;

proc sort data=lopnr_far_index_mor nodupkey;
by lopnr_far_index lopnr_mor;
run;

proc sort data=lopnr_far_index_mor out=lopnr_far_index_mor1 nodupkey;
by lopnr_far_index;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_t1d18_any_allpair2
	as select a.*, b.lopnr_mor as pre_id
	from fulsib_t1d18_any_allpair1 as a left join  lopnr_far_index_mor1 as b
	on a.lopnr_far_index=b.lopnr_far_index;
quit; 

data lopnr_mor_pre_id(keep=lopnr_mor pre_id);
set fulsib_t1d18_any_allpair2;
run;

proc sort data=lopnr_mor_pre_id nodupkey;
by lopnr_mor pre_id;
run;

proc sort data=lopnr_mor_pre_id out=lopnr_mor_pre_id1 nodupkey;
by lopnr_mor;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_t1d18_any_allpair3
	as select a.*, b.pre_id as id
	from fulsib_t1d18_any_allpair2 as a left join  lopnr_mor_pre_id1 as b
	on a.lopnr_mor=b.lopnr_mor;
quit; 

data pre_id_id(keep=pre_id id);
set fulsib_t1d18_any_allpair3;
run;

proc sort data=pre_id_id nodupkey;
by pre_id id;
run;

proc sort data=pre_id_id out=pre_id_id1 nodupkey;
by pre_id;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_t1d18_any_allpair4
	as select a.*, b.id as familyid
	from fulsib_t1d18_any_allpair3 as a left join  pre_id_id1 as b
	on a.pre_id=b.pre_id;
quit; 

data fulsib_t1d18_any_allpair5(drop=pre_id id);
set fulsib_t1d18_any_allpair4;
run; /*to be exported as fulsib_t1d18_any_allpair.csv*/