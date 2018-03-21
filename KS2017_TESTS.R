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

ks8f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-8F kdrFL-R7 s kdrFL.ab1"
compareHap(ks8f, newCombinedNoGapAllAdded)
#6915 MATCH FOUND 
#run from 6904-7021
#kdr3_(super-kdr2).seq kdr3_(super-kdr2).seq

ks9f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-9F kdrFL-R7 s kdrFL.ab1"
compareHap(ks9f, newCombinedNoGapAllAdded)
testMatch(ks9f,6913)
#0 rows 0 col #due to messy beginning

ks10f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-10F kdrFL-R7 s kdrFL.ab1"
compareHap(ks10f, newCombinedNoGapAllAdded)
testMatch(ks10f,6915)
#0rows 0 col

ks11f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-11F kdrFL-R7 s kdrFL.ab1"
compareHap(ks11f, newCombinedNoGapAllAdded)
#6913 match found 
#run from 6904-7021
#kdr1_(super-kdr3).seq ad kdr1_(super-kdr3)



testMatch <- function(file,indexmain) {
  sangerobj <- readsangerseq(ks10f) #read in file
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  b = addPriSec2(cutpri,cutsec)
  a = newCombinedNoGapAllAdded[6915,1]
  print(a)
  print(b)
  z = pairwiseAlignment(a,b)
  writePairwiseAlignments(z)
}

#Round 1
ks4f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-4F kdrFL-R7 s kdrFL.ab1"
ks6f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-6F kdrFL-R7 s kdrFL.ab1"
ks7f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-7F kdrFL-R7 s kdrFL.ab1"
ks12f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-12F kdrFL-R7 s kdrFL.ab1"
ks14f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-14F kdrFL-R7 s kdrFL.ab1"
ks15f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-15F kdrFL-R7 s kdrFL.ab1"
ks16f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-16F kdrFL-R7 s kdrFL.ab1"

ks18m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-18M kdrFL-R7 s kdrFL.ab1"
ks19m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-19M kdrFL-R7 s kdrFL.ab1"
ks20m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-20M kdrFL-R7 s kdrFL.ab1"
ks21m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-21M kdrFL-R7 s kdrFL.ab1"
ks22m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-22M kdrFL-R7 s kdrFL.ab1"
ks24m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-24M kdrFL-R7 s kdrFL.ab1"
ks25m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-25M kdrFL-R7 s kdrFL.ab1"
ks26m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-26M kdrFL-R7 s kdrFL.ab1"
ks30m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-30M kdrFL-R7 s kdrFL.ab1"
ks31m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-31M kdrFL-R7 s kdrFL.ab1"


sink("C:/Users/amw346/Desktop/compareHapOutputKansas.txt")
compareHap(ks4f, newCombinedNoGapAllAdded)
compareHap(ks6f, newCombinedNoGapAllAdded)
compareHap(ks7f, newCombinedNoGapAllAdded)
compareHap(ks12f, newCombinedNoGapAllAdded)
compareHap(ks14f, newCombinedNoGapAllAdded)
compareHap(ks15f, newCombinedNoGapAllAdded)
compareHap(ks16f, newCombinedNoGapAllAdded)
compareHap(ks18m, newCombinedNoGapAllAdded)
compareHap(ks19m, newCombinedNoGapAllAdded)
compareHap(ks20m, newCombinedNoGapAllAdded)
compareHap(ks21m, newCombinedNoGapAllAdded)
compareHap(ks22m, newCombinedNoGapAllAdded)
compareHap(ks24m, newCombinedNoGapAllAdded)
compareHap(ks25m, newCombinedNoGapAllAdded)
compareHap(ks30m, newCombinedNoGapAllAdded)
compareHap(ks31m, newCombinedNoGapAllAdded)

#round two
ks13f ="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-13F kdrFL-R7 s kdrFL.ab1"
ks23m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-23M kdrFL-R7 s kdrFL.ab1"
ks26m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-26M kdrFL-R7 s kdrFL.ab1"
ks27m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-27M kdrFL-R7 s kdrFL.ab1"
ks28m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-28M kdrFL-R7 s kdrFL.ab1"
ks29m="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-29M kdrFL-R7 s kdrFL.ab1"

sink("C:/Users/amw346/Desktop/compareHapOutputKansasRound2.txt")
compareHap(ks13f, newCombinedNoGapAllAdded)
compareHap(ks23m, newCombinedNoGapAllAdded)
compareHap(ks26m, newCombinedNoGapAllAdded)
compareHap(ks27m, newCombinedNoGapAllAdded)
compareHap(ks28m, newCombinedNoGapAllAdded)
compareHap(ks29m, newCombinedNoGapAllAdded)
sink()

