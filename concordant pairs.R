#Supplementary Table 2: to calculate concordant and discordant pairs
datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d18_any_allpair.csv', sep=";" , header=TRUE)
table(datWide$t1d18_index,datWide$t1d18) 

datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d6_any_allpair.csv', sep="," , header=TRUE)
table(datWide$t1d6_index,datWide$t1d6) 

datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d712_any_allpair.csv', sep="," , header=TRUE)
table(datWide$t1d712_index,datWide$t1d712) 

datWide <- read.csv('W:/C6_Carlsson/Yuxia Wei/data/familial/fulsib_t1d1318_any_allpair.csv', sep="," , header=TRUE)
table(datWide$t1d1318_index,datWide$t1d1318) 
