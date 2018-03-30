#change threshold of compare hap to 2 mismatches
library(stringi)

file="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"
#examples of how to use function and arguments
file="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1"
compareHap(file, newCombinedNoGap)


#warnings
#addbases requires capital letters
#you should run library(stringi)before hand otherwise stri_detect_fixed will throw error
compareHap<- function(file, master) {
  #inputs a chromatogram file and outputs the list of matches with the master list and the sequence that was a match
  sangerobj <- readsangerseq(file) #read in file
  print(file)
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  combined = addPriSec2(cutpri,cutsec)
  print(combined)
  #look for matches within the master list supplied to function
  match = MODcheckMasterList4(combined,master)
  return (match)
}



#HELPER FUNCTIONS FOR REFERENCE
clipIndex <- function(sangob) { 
  pri = primarySeq(sangob,string = TRUE) 
  len = nchar(pri) 
  num = 0 
  found = FALSE 
  for (i in 1:len) { 
    c1= substr(pri,i,i) 
    c2 = substr(pri,i+1,i+1) 
    c3 = substr(pri,i+2,i+2) 
    if (found) { 
      return(num-2) 
    } 
    if ((c1 == "N") & (c2 == "N") & (c3 == "N")) { 
      found = TRUE 
    } 
    num = i +1 
  } 
  return (num-2)
}


MODcheckMasterList4<- function(newseq, master) {
  matchesList = data.frame()
  len = dim(master)[1]
  count = 1
  for (i in 200:300) {
    found = FALSE
    overlap = overlapIsTrue(newseq,master[i,1])
    
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE) |stri_detect_fixed(newseq, master[i,1], case_insensitive = TRUE) | overlap )
    print(i)
    if (found == TRUE) {
      matchesList[count,1]= master[i,1]
      matchesList[count,2]= master[i,2]
      matchesList[count,3]= master[i,3]
      count = 1+count
      print(i)
      print("match found")
      
    }
  }
  
  return(matchesList)
}



MODcheckMasterList5<- function(newseq, master) {
  matchesList = data.frame()
  count = 1
  for (i in 1:7021) {
    found = FALSE
    alignseq = pairwiseAlignment(newseq,master[i,1])
    sum = summary(alignseq)
    found = (sum@nmismatch < 3) & (sum@ninsertion[2] == 0) & (sum@ndeletion[2] == 0)
   
    if (found == TRUE) {
      matchesList[count,1]= master[i,1]
      matchesList[count,2]= master[i,2]
      matchesList[count,3]= master[i,3]
      count = 1+count
      print(i)
      print("match found")
      
    }
  }
  
  return(matchesList)
}


w = "AATTTCCTTTCCATCTATCTATCTATCATMMDMDMDMDMDMDDDMDM"
C = "AATTTCCTTTCCATCTATCTATCTATCATYYYYYYYYYYYYYYYYYY"
u = pairwiseAlignment(w,C)
b8 = writePairwiseAlignments(u)
g =summary(u)
g@nmismatch

addPriSeq2 <- function(pri,sec) {
  #addPriSeq inputs two strings of same length and outputs a string that is the combined version according to IUPAC codes
  
  #checking if same length otherwise error message
  len = nchar(cutpri)
  len2 = nchar(cutsec)
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

#add bases for capital letters
addbases2<- function(a,b) {
  #addbases inputs two single character strings and outputs a single letter string according to IUPAC codes
  if (a == "A") {
    if (b == "A") {return ("A")}
    if (b == "C") {return ("M")}
    if (b == "G") {return ("R")}
    if (b == "T") {return ("W")}
    if (b == "N") {return ("N")}
    if (b == "R") {return ("R")}
    if (b == "Y") {return ("H")}
    if (b == "S") {return ("V")}
    if (b == "W") {return ("W")}
    if (b == "K") {return ("D")}
    if (b == "M") {return ("M")}
  }
  if (a == "C") {
    if (b == "A") {return ("M")}
    if (b == "C") {return ("C")}
    if (b == "G") {return ("S")}
    if (b == "T") {return ("Y")}
    if (b == "N") {return ("N")}
    if (b == "R") {return ("V")}
    if (b == "Y") {return ("Y")}
    if (b == "S") {return ("S")}
    if (b == "W") {return ("H")}
    if (b == "K") {return ("B")}
    if (b == "M") {return ("M")}
  }
  if (a == "G") {
    if (b == "A") {return ("R")}
    if (b == "C") {return ("S")}
    if (b == "G") {return ("G")}
    if (b == "T") {return ("K")}
    if (b == "N") {return ("N")}
    if (b == "R") {return ("R")}
    if (b == "Y") {return ("B")}
    if (b == "S") {return ("S")}
    if (b == "W") {return ("D")}
    if (b == "K") {return ("K")}
    if (b == "M") {return ("V")}
  }
  if (a == "T") {
    if (b == "A") {return ("W")}
    if (b == "C") {return ("Y")}
    if (b == "G") {return ("K")}
    if (b == "T") {return ("T")}
    if (b == "N") {return ("N")}
    if (b == "R") {return ("D")}
    if (b == "Y") {return ("Y")}
    if (b == "S") {return ("B")}
    if (b == "W") {return ("W")}
    if (b == "K") {return ("K")}
    if (b == "M") {return ("H")}
  }
  if (a == "N") {return("N")}
  return ("input strings not valid")	
}

#add bases for lowercase situation
addbases1<- function(a,b) {
  #addbases inputs two single character strings and outputs a single letter string according to IUPAC codes
  if (a == "a") {
    if (b == "a") {return ("A")}
    if (b == "c") {return ("M")}
    if (b == "g") {return ("R")}
    if (b == "t") {return ("W")}
    if (b == "n") {return ("N")}
    if (b == "r") {return ("R")}
    if (b == "y") {return ("H")}
    if (b == "s") {return ("V")}
    if (b == "w") {return ("W")}
    if (b == "k") {return ("D")}
    if (b == "m") {return ("M")}
  }
  if (a == "c") {
    if (b == "a") {return ("M")}
    if (b == "c") {return ("C")}
    if (b == "g") {return ("S")}
    if (b == "t") {return ("Y")}
    if (b == "n") {return ("N")}
    if (b == "r") {return ("V")}
    if (b == "y") {return ("Y")}
    if (b == "s") {return ("S")}
    if (b == "w") {return ("H")}
    if (b == "k") {return ("B")}
    if (b == "m") {return ("M")}
  }
  if (a == "g") {
    if (b == "a") {return ("R")}
    if (b == "c") {return ("S")}
    if (b == "g") {return ("G")}
    if (b == "t") {return ("K")}
    if (b == "n") {return ("N")}
    if (b == "r") {return ("R")}
    if (b == "y") {return ("B")}
    if (b == "s") {return ("S")}
    if (b == "w") {return ("D")}
    if (b == "k") {return ("K")}
    if (b == "m") {return ("V")}
  }
  if (a == "t") {
    if (b == "a") {return ("W")}
    if (b == "c") {return ("Y")}
    if (b == "g") {return ("K")}
    if (b == "t") {return ("T")}
    if (b == "n") {return ("N")}
    if (b == "r") {return ("D")}
    if (b == "y") {return ("Y")}
    if (b == "s") {return ("B")}
    if (b == "w") {return ("W")}
    if (b == "k") {return ("K")}
    if (b == "m") {return ("H")}
  }
  if (a == "n") {return("N")}
  return ("input strings not valid")	
}

overlapIsTrue <- function(long, short) {
  #function inputs two strings and checks to see if they have a matching overlapping portion
  #returns true if match is found
  
  t = pairwiseAlignment(long, short) # align strings to see where they overlap
  
  shortindex = t@subject@range@start
  longindex = t@pattern@range@start
  
  if (shortindex < longindex) {
    cutshort = substr(short,1,nchar(long)-longindex +1)
    cutlong = substr(long,longindex, nchar(long))
  } 
  
  if (longindex == shortindex) {
    if (nchar(long)> nchar(short)) {
      cutshort = substr(short,1,nchar(short))
      cutlong = substr(long,1, nchar(short))
    } else {
      cutshort = substr(short,1,nchar(long))
      cutlong = substr(long,1, nchar(long))
    }
  }
  if (longindex < shortindex) {
    cutshort = substr(short,shortindex,nchar(short))
    cutlong = substr(long,1,nchar(short)-shortindex+1)
  }
  if (cutlong ==cutshort) {
    return (TRUE)
  }
  return (FALSE)
}


#trying to learn how agrep works
a = "TCATGATYGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCMTGATKGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTTCTTCTTGGCCACGGTCGTGATCGGCAATYTTGTGGTAAGTTGACGTGGCCGAAACTRCTCCCCCGCTCCCAGGATGGAGGCTTCTGATGGCCWATTAAAAAAMMWTWWATYWWMMWYYMWCYTYYYYYYTYYYYYMWCTYYMWWCYKTMYWYYCKYYRCDGKTKCWKRWTCTTWWYYTWTYYTTRSYTTTGYYYHWSTYCRRMTTCGGATTTT"
b = "YTTGTGGTAAGTTGACGTGGCCGAAACTRCTCCCCCGCTCCCAGGATGGAGGCTTCTGATGGCCWATTAAAAAAMMWTWWATYWWMMWYYMWCYTYYYYYYTYYYYYMWCTYYMWWCYKTMYWYYCKYYRCWGKTKCWK"
c = "ACTRCTCCCCCGCTCCCAGGATGGAGGCTTCTGATGGCCWATTAAAAAAMMWTWWATYWWMMWYYMWCYTYYYYYYTYYYYYMWCTYYMWWCYKTMYWYYCKYYRCDGKTKCWKRWTCTTWWYYTWTYYTTRSYTTTGYYYHWSTYCRRMTTCGGATTTT"
agrep(c, b, max.distance = list(all = 2))
agrep("aaa", "xaab",max.distance = list(sub = 2))


testMatch <- function(file,indexmain) {
  sangerobj <- readsangerseq(file16) #read in file
  index = clipIndex(sangerobj) #cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) 
  cutsec = substr(sec,20,index)
  b = addPriSec2(cutpri,cutsec)
  a = newCombinedNoGapAllAdded[1126,1]
  print(a)
  print(b)
  z = pairwiseAlignment(a,b)
  writePairwiseAlignments(z)
  str(z)
  z[patter]
  z@pattern@mismatch
  z@subject@indel
  z@subject@unaligned@quality
}

z@nmatch
s = summary(z)
s@nmatch

nchar(a)
nchar(b)
str(z)
str(summary(z))
atrr()
s@ninsertion

b="TCATGTATYGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCMTGTTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTRCTCCCCCGCTCCCAGGATGGAGGCTTCWKMYGKMMWATWMAWAAAWWTKAMATYWAYCTCTCTCTTTCTCTCTCYCWMMWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
a2 = "GAATTTCACCGACTTCATGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATTTTGTGGTAAGTTGACGTGGCCGAAACTRCTCCCCCGCTCCCAGGATGGAGGCTTCWKMYGKMMWATWMAWAAAWWTKAMATYWAYCTCTCTCTTTCTCTC"
r = pairwiseAlignment(b,a2)
writePairwiseAlignments(r)
y = summary(r)
y@ninsertion
y@ndeletion
y@nmismatch

r2 = pairwiseAlignment(a2,b)
writePairwiseAlignments(r2)
y2 = summary(r2)
y2@ninsertion


testfile = "C:/Users/amw346/Desktop/8KSF1kdr2.ab1"
sangerobj <- readsangerseq(testfile) #read in file
print(file)
index = clipIndex(sangerobj) #cut off all Ns at the end
#221
seqcutoff = "TNCGACTNAGCACAGCTTCATGATTGTGGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCANGNTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTTCTTCTTGGCCACGGTCGTGATCGGCAATCNTGTGGTAAGTTGACGTGGCCGAAACTGCTCCCCCGCTCCCAGGATGGAGGCTTCAGANGNCCAANTAAANANANTTAANNCNTNNN"
nchar(seqcutoff)
#run again its still 221
basecalls <- makeBaseCalls(sangerobj)
pri <- primarySeq(basecalls, string = 'TRUE')
sec <- secondarySeq(basecalls, string = 'TRUE')
cutpri = substr(pri,20,index) 
cutsec = substr(sec,20,index)
combined = addPriSec2(cutpri,cutsec)
print(combined)
sangerobj <- readsangerseq(testfile) #read in file
print(file)
index = clipIndex(sangerobj)
#221

testfile3 = "C:/Users/amw346/Desktop/Batch 1/4KSF1kdr.ab1"
sangerobj <- readsangerseq(testfile3) #read in file
index = clipIndex(sangerobj)
basecalls <- makeBaseCalls(sangerobj)
pri <- primarySeq(basecalls, string = 'TRUE')
sec <- secondarySeq(basecalls, string = 'TRUE')
cutpri = substr(pri,20,index) 
cutsec = substr(sec,20,index)
combined = addPriSec2(cutpri,cutsec)
index

#index seems to remain the same, problem remains
compareHap(testfile3, newCombinedNoGapAllAdded)

sewe = combined
alignseq2 = pairwiseAlignment(sewe,newCombinedNoGapAllAdded[1275,1])
writePairwiseAlignments(alignseq2)
sum2 = summary(alignseq2)
found = (sum2@nmismatch < 3) & (sum2@ninsertion[2] == 0) & (sum2@ndeletion[2] == 0)
