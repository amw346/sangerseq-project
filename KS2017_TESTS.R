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
#kdr1_(super-kdr3).seq ad kdr1_(super-kdr3).seq

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