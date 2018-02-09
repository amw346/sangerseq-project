#Testing compareHap(file, newCombinedNoGap)
#compareHap function located in compareHapFinalVersion.R



#1FLF1kdr1 kdr1/kdr-his1

file="C:/Users/amw346/Desktop/1FLF1kdr.ab1"
compareHap(file, newCombinedNoGap)
b = newCombinedNoGap[9,1]
a=combined
x=pairwiseAlignment(a,b)
writePairwiseAlignments(x)
#174, 181, 189 on chromatogram

#ERROR
#1FLF2 V20/KDR4
file2="C:/Users/amw346/Desktop/1FLF2kdr1.ab1"
compareHap(file2, newCombinedNoGap)
g =combined
b = newCombinedNoGap[1357,1]
z=pairwiseAlignment(g,b)
writePairwiseAlignments(z)

file3="C:/Users/amw346/Desktop/1FLF1kdrother.ab1"
z=pairwiseAlignment(combined3,b)
writePairwiseAlignments(z)

#2F1F1kdr kdr2 v 39
file4="C:/Users/amw346/Desktop/2FLF1kdr.ab1"
z=pairwiseAlignment(combined4,h)
b = newCombinedNoGap[1166,1]
writePairwiseAlignments(z)


# 10FLF2KDR V40 KDRHIS1
file5="C:/Users/amw346/Desktop/10FLF2kdr.ab1"
z=pairwiseAlignment(combined5,newCombinedNoGap[53,1])
writePairwiseAlignments(z)
#205 A vs W

WTF(file5,53)
WTF <- function(file,index) {
  sangerobj <- readsangerseq(file) #read in file
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  b = addPriSec2(cutpri,cutsec)
  a = newCombinedNoGap[index,1]
  a
  b
  z = pairwiseAlignment(a,b)
  writePairwiseAlignments(z)
}






