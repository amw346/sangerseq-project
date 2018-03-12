#Testing NewSequences

library(stringi)

ks1f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-1F kdrFL-R7 s kdrFL.ab1"
compareHap(ks1f, newCombinedNoGapAllAdded)
#o rows 0 col


ks2f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-2F kdrFL-R7 s kdrFL.ab1"
compareHap(ks2f, newCombinedNoGapAllAdded)
#6913 match found 
#run from 6904-7021
#kdr1_(super-kdr3).seq ad kdr1_(super-kdr3).seq

ks3f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-3F kdrFL-R7 s kdrFL.ab1"
compareHap(ks3f, newCombinedNoGapAllAdded)
#6915 MATCH FOUND 
#run from 6904-7021
#kdr3_(super-kdr2).seq kdr3_(super-kdr2).seq

ks5f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-5F kdrFL-R7 s kdrFL.ab1"
compareHap(ks5f, newCombinedNoGapAllAdded)
#6915 MATCH FOUND 
#run from 6904-7021
#kdr3_(super-kdr2).seq kdr3_(super-kdr2).seq

