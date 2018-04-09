#Test sequences for match with 7029
library(stringi)

KS17_18Jul17_12F_colony1_kdrFL_R7 = "TCGCTTCAAGGACCATGAATTACCGCGCTGGAATTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTACTCCCCCGCTCCCAGGATGGAGGCTTCATCCGTAATATACATAAATTTGACATTTATCTCTCTCTTTCTCTCTCCCAACTTTATTCTCTCCACTGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"

compareHap4string<- function(combined, master) {
  #takes string and compares to sequence in master- must already be single sequence
  print(combined)
  
  #look for matches within the master list supplied to function
  match = MODcheckMasterList6(combined,master)
  return (match)
}

MODcheckMasterList6<- function(newseq, master) {
  matchesList = data.frame()
  len = dim(master)[1]
  count = 1
  for (i in 6903:7021) {
    print(i)
    found = FALSE
    overlap = overlapIsTrue(newseq,master[i,1])
    
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE) |stri_detect_fixed(newseq, master[i,1], case_insensitive = TRUE) | overlap )
    
    if (found == TRUE) {
      matchesList[count,1]= master[i,1]
      matchesList[count,2]= master[i,2]
      matchesList[count,3]= master[i,3]
      print(i)
      print("match found")
      
    }
  }
  count = 1+count
  
  return(matchesList)
}

addPriSeq3 <- function(pri,sec) {
  #addPriSeq inputs two strings of same length and outputs a string that is the combined version according to IUPAC codes
  
  #checking if same length otherwise error message
  len = nchar(pri)
  len2 = nchar(sec)
  if (len != len2) {return ("input sequences different lengths")}
  
  #adding the elements
  combined = ""
  for (i in 1:len) {
    newchar = addbases2(substr(pri,i,i), substr(sec,i,i))
    combined= paste0(combined, newchar)
  }
  
  #returning the added sequence
  return (AddedSeq = combined)
}

#tests of 14f 12 f sequences
#test1
compareHap4string(KS17_18Jul17_12F_colony1_kdrFL_R7, newCombinedNoGapAllAdded7021)
#match found
#kdr3_(super-kdr2).seq

#test 2
KS17_18Jul17_12F_colony2_kdrFL_R7 = "TCGCTTCAAGGACCATGAATTACCGCGCTGGAATTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTACTCCCCCGCTCCCAGGATGGAGGCTTCATCCGTAATATACATAAATTTGACATTTATCTCTCTCTTTCTCTCTCCCAACTTTATTCTCTCCACTGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"
compareHap4string(KS17_18Jul17_12F_colony2_kdrFL_R7, newCombinedNoGapAllAdded7021)
#matchfound
#kdr3_(superkdr2)

#test3
KS17_18Jul17_12F_colony3_kdrFL_R7 = "TCGCTCAAGGACCATGAATTACCGCGCTGGAACTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGTATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTGCTCCCAGGATGGAGGCTTCTGATGGCCAATTAAAAAAAATTAAATCAACCTCTCTCTTTCTCTCTCTCTCAACTTTATTCCGTCCATCCGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"
compareHap4string(KS17_18Jul17_12F_colony3_kdrFL_R7, newCombinedNoGapAllAdded7021)
#matchfound
#kdr1_(super-kdr3)

#test 4 
KS17_18Jul17_14F_colony1_kdrFL_R7 = "TCGCTTCAAGGACCATGAATTACCGCGCTGGAATTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTACTCCCCCGCTCCCAGGATGGAGGCTTCATCCGTAATATACATAAATTTGACATTTATCTCTCTCTTTCTCTCTCCCAACTTTATTCTCTCCACTGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"
compareHap4string(KS17_18Jul17_14F_colony1_kdrFL_R7, newCombinedNoGapAllAdded7021)
#matchfound
#kdr3_(superkdr2)

#test 5 
KS17_18Jul17_14F_colony2_kdrFL_R7 = "TCGCTTCAAGGACCATGAATTACCGCGCTGGAATTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTACTCCCCCGCTCCCAGGATGGAGGCTTCATCCGTAATATACATAAATTTGACATTTATCTCTCTCTTTCTCTCTCCCAACTTTATTCTCTCCACTGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"
compareHap4string(KS17_18Jul17_14F_colony2_kdrFL_R7, newCombinedNoGapAllAdded7021)
#match found
#kdr3_superkdr2

#test6
KS17_18Jul17_14F_colony3_kdrFL_R7 = "TCGCTTCAAGGACCATGAATTACCGCGCTGGAACTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGTATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTGCTCCCAGGATGGGGGCTTCTGATGGCCAATTAAAAAAAATTAAATCAACCTCTCTCTTTCTCTCTCTCTCAACTTTATTCCGTCCATCCGTTGCAGGTTCTTAATCTTTTCTTAGCTTTGCTTTTGTCCAACTTCGGTTCATCTAGTTTATCAGCCCCGACTGCCGACAATGATACCA"
compareHap4string(KS17_18Jul17_14F_colony3_kdrFL_R7, newCombinedNoGapAllAdded7021)
#nomatch


#observing why nomatch was found for these sequences
#add12F colony 1 and 12F colony 3

nchar(KS17_18Jul17_12F_colony1_kdrFL_R7)
nchar(KS17_18Jul17_12F_colony3_kdrFL_R7)

z1 = pairwiseAlignment(KS17_18Jul17_12F_colony1_kdrFL_R7,KS17_18Jul17_12F_colony3_kdrFL_R7)

writePairwiseAlignments(z1)

#get them to be the same length
KS17_18Jul17_12F_colony1_kdrFL_R7_short = substr(KS17_18Jul17_12F_colony1_kdrFL_R7,1,378)
combinedcolony1_3 = addPriSeq3(KS17_18Jul17_12F_colony1_kdrFL_R7_short, KS17_18Jul17_12F_colony3_kdrFL_R7)
nchar(combinedcolony1_3)

#generated sequence from 12F chromatogramfile
ks12f="C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-12F kdrFL-R7 s kdrFL.ab1"

sangerobj <- readsangerseq(ks12f) #read in file
index = clipIndex(sangerobj) #cut off all Ns at the end

#cut and add combined string
basecalls <- makeBaseCalls(sangerobj)
pri <- primarySeq(basecalls, string = 'TRUE')
sec <- secondarySeq(basecalls, string = 'TRUE')
cutpri = substr(pri,20,index) 
cutsec = substr(sec,20,index)
combined = addPriSeq3(cutpri,cutsec)
print(combined)

z = pairwiseAlignment(combined, newCombinedNoGapAllAdded7021[1019,1])
writePairwiseAlignments(z)

a = "CGAATNCGACTTCNGCACAGCATTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTGCTCCCNCGNTNNAGGNATGGNAGGNTCANCNNAAAAANANAAAA"
nchar(a)


superkdr2_superkdr3_7029 = newCombinedNoGapAllAdded7021[1019,1]
SUPER2 = toupper(super2)

super3 = toupper(newCombinedNoGapAllAdded7021[1019,4])
super2 = newCombinedNoGapAllAdded7021[1019,5]

w4 = pairwiseAlignment(SUPER2, KS17_18Jul17_12F_colony1_kdrFL_R7)
writePairwiseAlignments(w4)

w6 = pairwiseAlignment(super3, KS17_18Jul17_12F_colony3_kdrFL_R7)
writePairwiseAlignments(w6)
#1019

w2 = pairwiseAlignment(superkdr2_superkdr3_7029,combinedcolony1_3)
writePairwiseAlignments(w2)

w5 = pairwiseAlignment(SUPER2, super3)
writePairwiseAlignments(w5)

nchar(super3)
nchar(SUPER2)
nchar(superkdr2_superkdr3_7029)

substr(SUPER2,1,276)

w7 = pairwiseAlignment(KS17_18Jul17_12F_colony1_kdrFL_R7, KS17_18Jul17_12F_colony3_kdrFL_R7)
writePairwiseAlignments(w7)

super2short = substr(SUPER2,8,276)
super3short = substr(super3,8,276)

shorter2and3 = addPriSeq3(super2short,super3short)
w10 = pairwiseAlignment(shorter2and3,newCombinedNoGapAllAdded7021[1019,1])
writePairwiseAlignments(w10)


filenamesforMD <- list.files("C:/Users/amw346/Desktop/mdfiels", full.name = TRUE)
sink("C:/Users/amw346/Desktop/comparehapmdseq.txt")
for (i in 1:32) {
  compareHap(filenamesforMD[i],newCombinedNoGapAllAdded7021)
}







