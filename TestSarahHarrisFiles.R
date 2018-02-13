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


testMatch <- function(file,indexmain) {
  sangerobj <- readsangerseq(file8) #read in file
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj, ratio = 0.1)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  b = addPriSec2(cutpri,cutsec)
  a = newCombinedNoGap[indexmain,1]
  print(a)
  print(b)
  z = pairwiseAlignment(a,b)
  writePairwiseAlignments(z)
}

#
file5="C:/Users/amw346/Desktop/1KSF1kdr.ab1"
compareHap(file5, newCombinedNoGap)
testMatch(file5, 1152)

file6="C:/Users/amw346/Desktop/Batch 1/1KSF1kdr.ab1"
testMatch(file6,1152)

file7 = "C:/Users/amw346/Desktop/Batch 1/2KSF1kdr.ab1"
#doont have answer testMatch(file7)

file8 = "C:/Users/amw346/Desktop/Batch 1/3KSF1kdr.ab1"
testMatch(file8,354)

file9 = "C:/Users/amw346/Desktop/Batch 1/4KSF1kdr.ab1"
testMatch(file9,1275)

file10 = "C:/Users/amw346/Desktop/Batch 1/5KSF1kdr.ab1"
testMatch(file10, 242)

file11 = "C:/Users/amw346/Desktop/Batch 1/8KSF1kdr.ab1"
testMatch(file11, 282)
compareHap(file10,newCombinedNoGap)

file12 = "C:/Users/amw346/Desktop/10KSF2kdr.ab1"
testMatch(file12, 4279)

file13 = "C:/Users/amw346/Desktop/Batch 1/11KSF1kdr.ab1"
testMatch(file13, 1145)

#test old file 
file14 = "C:/Users/amw346/Desktop/aa.ab1"
testMatch(file16,3999)

file14="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"
file16="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1"


compareHap(file16,newCombinedNoGap)


