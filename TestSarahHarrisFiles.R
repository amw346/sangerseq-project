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

#








