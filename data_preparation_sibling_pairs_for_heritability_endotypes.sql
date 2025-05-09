/*prepare sibling pairs used in the calculation of heritability at age 0-6, 7-12, and 13-18 years*/
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
run; /*2852600: 15126 exc_trend_t1d_any=1; 77263 multi_brith^="1"*/
/*index=1: reused_lopnr^=1 and dup_lopnr^=1*/ 

/*0-6 years*/
data index_t1d(keep=lopnrindex lopnr_mor date_birth_index men_index t1d6_index);
set mbr_t1d18_any1;
	if t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25<7 then t1d6_index=1;
	if t1d18_any^=1 or t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25>=7 then t1d6_index=0;
	rename lopnr_barn=lopnrindex;
	rename date_birth_barn=date_birth_index;
	rename men=men_index;
run;
proc freq data=index_t1d;
tables t1d6_index;
run; /*n=2852600, t1d6=5838*/
proc freq data=index_t1d;
tables men_index date_birth_index;
run;

data lopnr_sibling(keep=lopnr date_birth men t1d6);
set mbr_t1d18_any1;
	if t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25<7 then t1d6=1;
	if t1d18_any^=1 or t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25>=7 then t1d6=0;
	rename lopnr_barn=lopnr;
	rename date_birth_barn=date_birth;
run;
proc freq data=lopnr_sibling;
tables t1d6;
run; /*n=2852600, t1d6=5838*/
proc freq data=lopnr_sibling;
tables men date_birth;
run;

*restart sas*
/*7-12 years*/
data index_t1d(keep=lopnrindex lopnr_mor date_birth_index men_index t1d712_index);
set mbr_t1d18_any1;
	if t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25<13 then t1d712_index=1;
	if t1d18_any^=1 or t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25>=13 then t1d712_index=0;
	rename lopnr_barn=lopnrindex;
	rename date_birth_barn=date_birth_index;
	rename men=men_index;
where (date_t1d18_any-date_birth_barn)/365.25>=7;
run;
proc freq data=index_t1d;
tables t1d712_index;
run; /*n=2801146, 8470 t1d712=1*/

data lopnr_sibling(keep=lopnr date_birth men t1d712);
set mbr_t1d18_any1;
	if t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25<13 then t1d712=1;
	if t1d18_any^=1 or t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25>=13 then t1d712=0;
	rename lopnr_barn=lopnr;
	rename date_birth_barn=date_birth;
where (date_t1d18_any-date_birth_barn)/365.25>=7;
run;
proc freq data=lopnr_sibling;
tables t1d712;
run; /*n=2799431, 7790 t1d712=1*/
/*repeat the macro steps for age 0-6 years and export data as fulsib_t1d712_trend.csv*/

*restart sas*
/*13-18 years*/
data index_t1d(keep=lopnrindex lopnr_mor date_birth_index men_index t1d1318_index);
set mbr_t1d18_any1;
	if t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25<19 then t1d1318_index=1;
	if t1d18_any^=1 or t1d18_any=1 and (date_t1d18_any-date_birth_barn)/365.25>=19 then t1d1318_index=0;
	rename lopnr_barn=lopnrindex;
	rename date_birth_barn=date_birth_index;
	rename men=men_index;
where (date_t1d18_any-date_birth_barn)/365.25>=13;
run;
proc freq data=index_t1d;
tables t1d1318_index;
run; /*n=2403067, 5334 t1d1318*/

data lopnr_sibling(keep=lopnr date_birth men t1d1318);
set mbr_t1d18_any1;
	if t1d18_any=1 then t1d1318=1;
	if t1d18_any^=1 then t1d1318=0;
	rename lopnr_barn=lopnr;
	rename date_birth_barn=date_birth;
where (date_t1d18_any-date_birth_barn)/365.25>=13;
run;
proc freq data=lopnr_sibling;
tables t1d1318;
run;  /*n=2403067, 5334 t1d1318*/
/*for different age groups*/
%macro outcome(var1);
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
proc freq data=sibling_merged;
tables relation/missing;
run;

data sibling_merged1;
set sibling_merged;
where relation^="";
run; /*3378280 for 0-6 years; 3331523 for 7-12 years; 2884547 for age 13-18 years*/

proc sql;
create table sibling_match
	as select a.*, b.*
	from index_t1d as a left join  sibling_merged1 as b
	on a.lopnrindex=b.lopnrindex;
quit;
proc freq data=sibling_match;
tables relation/missing;
run; /*nonmissing relationship: 3102882 for age 0-6 years; 3032930 for age 7-12 years; 2551232 non-missing for age 13-18 years*/

data fulsib_&var1.; /*to be used*/
set sibling_match;
where relation^="";
run;
proc freq data=fulsib_&var1.;
tables relation/missing;
run; /*3102882 for age 0-6 years; 3032930 for age 7-12 years; 2551232 for age 13-18 years*/

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
run; 
proc freq data=lopnr_far;
tables relation;
run; /*4241381*/

proc sql;
create table fulsib_&var1._lopnrfar
	as select a.*, b.lopnr as lopnr_far_index
	from fulsib_&var1. as a left join  lopnr_far as b
	on a.lopnrindex=b.lopnrindex;
quit;

proc sql;
create table fulsib_&var1._lopnrfar1
	as select a.*, b.lopnr as lopnr_far
	from fulsib_&var1._lopnrfar as a left join  lopnr_far as b
	on a.lopnr=b.lopnrindex;
quit;
proc freq data=fulsib_&var1._lopnrfar1;
tables relation/missing;
where lopnr_far_index=. or lopnr_far=.;
run; /*4146 missing paternal lopnr for age 0-6 years, 2992 missing paternal lopnr for 7-12 years,
and 1908 missing for 13-18 years*/

/**********one sibling pair from each mother**********/
data fulsib_&var1._any(keep=lopnrindex lopnr_mor lopnr_far lopnr men_index men &var1._index &var1. date_birth_index date_birth relation);
set fulsib_&var1._lopnrfar1;
	where year(date_birth_index)<2011 and year(date_birth)<2011 and lopnr_far_index^=. and lopnr_far^=.;
run; 

proc sort data=fulsib_&var1._any;
by lopnr_mor date_birth_index date_birth;
run; /*3098736 (0-6), 3029938 (7-12), 2549324 (13-18) pairs of siblings. For each mother, only one pair of sibligns should be selected*/
proc sort data=fulsib_&var1._any out=fulsib_&var1._any0 nodupkey;
by lopnr_mor;
run; /*898038 (0-6), 884761 (7-12), 751249(13-18) pairs of siblings were selected by lopnr_mor*/
proc sort data=fulsib_&var1._any0;
by lopnr_far date_birth_index date_birth;
run; 
proc sort data=fulsib_&var1._any0 out=fulsib_&var1._any1 dupout=duplicated_lopnr nodupkey;
by lopnr_far;
run;  /*884975 (0-6), 872085 (7-12), 742220 (13-18) pairs remained*/

data fulsib_&var1._any;
set fulsib_&var1._any1;
	birth_year_index=year(date_birth_index);
	birth_year=year(date_birth);
run; 
proc freq data=fulsib_&var1._any;
tables &var1._index*&var1.;
run; /*84 concordant T1D pairs, 1680+1786 discordant pairs for age 0-6 years;
90 concordant T1D pairs, 2380+2660 discordant pairs for age 7-12 years;
40 concordant T1D pairs, 1603+1656 discordant pairs for age 13-18 years*/
/*to be exported as fulsib_&var1._any.csv*/
proc export data=fulsib_&var1._any 
outfile="W:\C6_Carlsson\Yuxia Wei\data\familial\fulsib_&var1._any.csv" 
dbms=csv 
replace;
run;  
/********************all possible siblings pairs************************************/
/*exclude children with missing lopnr_far*/
data fulsib_&var1._lopnrfar2;
set fulsib_&var1._lopnrfar1;
where lopnr_far_index^=. and lopnr_far^=.;
run; /*one pair of siblings appear twice in this dataset*/

data fulsib_&var1._any_allpair1;
set fulsib_&var1._lopnrfar2;
	where date_birth_index<date_birth and 
		year(date_birth_index)<2011 and year(date_birth)<2011
		and year(date_birth_index)>1981 and year(date_birth)>1981;
run; /*1549362 for t1d6; 1514963 for t1d712; 1274660 for t1d1318*/ 
data fulsib_&var1._any_allpair2;
set fulsib_&var1._lopnrfar2;
	where date_birth_index>date_birth and 
		year(date_birth_index)<2011 and year(date_birth)<2011
		and year(date_birth_index)>1981 and year(date_birth)>1981;
run; /*1549362 for t1d6; 1514963 for t1d712; 1274660 for t1d1318*/

data fulsib_&var1._any_allpair3;
set fulsib_&var1._any_allpair2;
	rename lopnr=lopnr2;
	rename lopnr_far=lopnr_far2;
	rename date_birth=date_birth2;
	rename men=men2;
	rename &var1.=&var1.2;
run; 

data fulsib_&var1._any_allpair4;
set fulsib_&var1._any_allpair3;
	rename lopnrindex=lopnr;
	rename lopnr_far_index=lopnr_far;
	rename men_index=men;
	rename date_birth_index=date_birth;
	rename &var1._index=&var1.;
run;
data fulsib_&var1._any_allpair5;
set fulsib_&var1._any_allpair4;
	rename lopnr2=lopnrindex;
	rename lopnr_far2=lopnr_far_index;
	rename date_birth2=date_birth_index;
	rename men2=men_index;
	rename &var1.2=&var1._index;
run;
data fulsib_&var1._any_allpair;
set fulsib_&var1._any_allpair1 fulsib_&var1._any_allpair5;
run;
data fulsib_&var1._any_allpair1;
set fulsib_&var1._any_allpair;
	birth_year_index=year(date_birth_index);
	birth_year=year(date_birth);
run;
proc sort data=fulsib_&var1._any_allpair1;
by lopnr_mor date_birth_index date_birth;
run;
proc sort data=fulsib_&var1._any_allpair1 nodupkey dupout=dup_exc;
by lopnr_mor date_birth_index date_birth;
run; /*1549357 pairs remained for t1d6(if twins are matched to different their siblings, only one person in the twin paris remained);
		1514958 pairs remained for t1d712; 1274660 pairs for t1d1318*/
proc freq data=fulsib_&var1._any_allpair1;
tables &var1._index*&var1.;
where relation=:"0010";
run; /*146 concordant pairs, 2785+3084 discordant paris for t1d6;
150 concordant pairs, 4000+4710 discordant pairs for t1d712;
65 concordant pairs, 2788+2838 discordant pairs for t1d1318*/

/*to get familyid*/
data lopnr_far_index_mor(keep=lopnr_far_index lopnr_mor);
set fulsib_&var1._any_allpair1;
run;

proc sort data=lopnr_far_index_mor nodupkey;
by lopnr_far_index lopnr_mor;
run;

proc sort data=lopnr_far_index_mor out=lopnr_far_index_mor1 nodupkey;
by lopnr_far_index;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_&var1._any_allpair2
	as select a.*, b.lopnr_mor as pre_id
	from fulsib_&var1._any_allpair1 as a left join  lopnr_far_index_mor1 as b
	on a.lopnr_far_index=b.lopnr_far_index;
quit; 

data lopnr_mor_pre_id(keep=lopnr_mor pre_id);
set fulsib_&var1._any_allpair2;
run;

proc sort data=lopnr_mor_pre_id nodupkey;
by lopnr_mor pre_id;
run;

proc sort data=lopnr_mor_pre_id out=lopnr_mor_pre_id1 nodupkey;
by lopnr_mor;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_&var1._any_allpair3
	as select a.*, b.pre_id as id
	from fulsib_&var1._any_allpair2 as a left join  lopnr_mor_pre_id1 as b
	on a.lopnr_mor=b.lopnr_mor;
quit; 

data pre_id_id(keep=pre_id id);
set fulsib_&var1._any_allpair3;
run;

proc sort data=pre_id_id nodupkey;
by pre_id id;
run;

proc sort data=pre_id_id out=pre_id_id1 nodupkey;
by pre_id;
run; /*all children with the same father have the same pre_id (first lopnr_mor in each lopnr_far)*/

proc sql;
create table fulsib_&var1._any_allpair4
	as select a.*, b.id as familyid
	from fulsib_&var1._any_allpair3 as a left join  pre_id_id1 as b
	on a.pre_id=b.pre_id;
quit; 

data fulsib_&var1._any_allpair5(drop=pre_id id);
set fulsib_&var1._any_allpair4;
run; /*to be exported as fulsib_&var1._any_allpair.csv*/
proc export data=fulsib_&var1._any_allpair5
outfile="W:\C6_Carlsson\Yuxia Wei\data\familial\fulsib_&var1._any_allpair.csv" 
dbms=csv 
replace;
run;  /*"any" means any childhood-onset t1d regardless of conflicting types of diagnosis (t1d18_any)*/

%mend;
%outcome(t1d6);
%outcome(t1d712);
%outcome(t1d1318);



