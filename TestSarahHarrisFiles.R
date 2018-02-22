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


testMatch <- function(file,indexmain=354) {
  sangerobj <- readsangerseq(file8) #read in file
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  b = addPriSec2(cutpri,cutsec)
  a = newCombinedNoGapAllAdded[indexmain,1]
  print(a)
  print(b)
  z = pairwiseAlignment(a,b)
  writePairwiseAlignments(z)
}
nchar(b)
#from kansas batch 1 folder
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

#test old file # v3v5 and
file16="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1" 
testMatch(file16,3999)

file14="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"
file16="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1"


compareHap(file14,newCombinedNoGapAllAdded)


#fuzzymatrix test
sangerobj <- readsangerseq(file16) #read in file
index = clipIndex(sangerobj) #cut off all Ns at the end

#cut and add combined string
basecalls <- makeBaseCalls(sangerobj)
pri <- primarySeq(basecalls, string = 'TRUE')
sec <- secondarySeq(basecalls, string = 'TRUE')
cutpri = substr(pri,20,index) 
cutsec = substr(sec,20,index)
b = addPriSec2(cutpri,cutsec)
a = newCombinedNoGap[3999,1]
d=pairwiseAlignment(b,a)
writePairwiseAlignments(d)


nmismatch(d)
nmatch(d)
mapping <- diag(4)
dimnames(mapping) <- list(DNA_BASES, DNA_BASES)
mapping["C", "T"] <- mapping["T", "C"] <- 1
mapping["G", "A"] <- mapping["A", "G"] <- 1


#ratio 1/3 

#testcase1
#v26/kdr2
file6="C:/Users/amw346/Desktop/Batch 1/1KSF1kdr.ab1"
#mod4 range: 1100-1200
compareHap(file6, newCombinedNoGapAllAdded)
#0 rows 0 cols
testMatch(file6,1152)
#273 on chromatagram

#testcase2
#super-kdr2/superkdr3
file8 = "C:/Users/amw346/Desktop/Batch 1/3KSF1kdr.ab1"
#mod4range: 350-400
compareHap(file8, newCombinedNoGapAllAdded)
#0 rows 0 col
testMatch(file8,354)
#167, 181, 189

#testcase3
#v41/superkdr2
file9 = "C:/Users/amw346/Desktop/Batch 1/4KSF1kdr.ab1"
#mod4 range: 1200-1400
compareHap(file9, newCombinedNoGapAllAdded)
#match found 1275
testMatch(file9,1275)

#testcase4 
#kdr-his3/super-kdr2
file10 = "C:/Users/amw346/Desktop/Batch 1/5KSF1kdr.ab1"
#mod4 range: 200-400
compareHap(file10, newCombinedNoGapAllAdded)
#0 rows 0 columns
testMatch(file10, 242)
#208

#testcase5
#v39/ kdr-his3
file11 = "C:/Users/amw346/Desktop/Batch 1/8KSF1kdr.ab1"
#mod4: 250-300
compareHap(file11, newCombinedNoGapAllAdded)
#kdr-his3/v30 AND kdr-his3/v39
testMatch(file11, 282)

#testcase6
#v33/v37
file12 = "C:/Users/amw346/Desktop/10KSF2kdr.ab1"
#mod4: 4200,4300
compareHap(file12, newCombinedNoGapAllAdded)
#0 rows 0 col
testMatch(file12, 4279)
#203

#testcase7
#super-kdr2/super-kdr2
file13 = "C:/Users/amw346/Desktop/Batch 1/11KSF1.ab1"
#mod4: 1140,1160
compareHap(file13, newCombinedNoGapAllAdded)
#0 row 0 col
testMatch(file13, 1145)


#0 rows 0 col
testMatch(file12, 4279)
#203

