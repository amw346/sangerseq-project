#COMPARE A FILE TO THE MASTERLIST AND FIND POSSIBLE COMBOS
library(stringi)
file="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"

file="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1"

#addbases requires capital letters
compareHap<- function(file, master) {
  sangerobj <- readsangerseq(file)
  index = clipIndex(sangerobj)
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index)
  cutsec = substr(sec,20,index)
  combined = addPriSec2(cutpri,cutsec)
  match = MODcheckMasterList4(combined,master)
  return (match)
}

compareHap("C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1",TESTcombinedtypes)
compareHap(file, newCombinedNoGap)



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
  for (i in 3000:4100) {
    found = FALSE
    overlap1 = overlapIsTrue(newseq,master[i,1])
    overlap2 = overlapIsTrue(master[i,1],newseq)
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE) |stri_detect_fixed(newseq, master[i,1], case_insensitive = TRUE) | overlap1 | overlap2)
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

file1="C:/Users/amw346/Desktop/aa.ab1"
aabysseq <- readsangerseq(file1)
abseq = "ATCCGACTTCAGCACAGCTTCATGATTGTGTTCCGAGTGCTGTGCGGAGAGTGGATCGAGTCCATGTGGGACTGCATGTATGTGGGCGATGTCAGCTGTATACCCTTCTTCTTGGCCACGGTCGTGATCGGCAATCTTGTGGTAAGTTGACGTGGCCGAAACTGCTCYCCCGCTCCCAGGATGGAGGCTTCWKMYGKMMWATWMAAAAAWWWKWMAWTYAWCYYYYYYYTTYYCYCYYCYMWCTYWAYTYTMTYCMCWGYWKCMKKTKCWTRWTCTTWWYYTWKYYTTRCYTTTGYYYWWSTYCRRYTCMTSTASWTYWWCWKYMYCRACYSCCRAYRMYGAYAMYRATACCAA"
v2v3combo=  "CTTGTGGTAAGTTGACGTGGCCGAAACTGCTCYCCCGCTCCCAGGATGGAGGCTTCWKMYGKMMWATWMAAAAAWWWKWMAWTYAWCYYYYYYYTTYYYYYYYCYMWCTYWAYTYTMTYCMSWSYWKCMK"
v= pairwiseAlignment(abseq,v2v3combo)
writePairwiseAlignments(v)

v3v5combo = newCombinedNoGap[3999,1]
v= pairwiseAlignment(abseq,v3v5combo)
writePairwiseAlignments(v)

inputseqaabys = combined

v= pairwiseAlignment(inputseqaabys,v3v5combo)
writePairwiseAlignments(v)

v= pairwiseAlignment(inputseqaabys,v2v3combo)
writePairwiseAlignments(v)

found = (stri_detect_fixed(v3v5combo,inputseqaabys, case_insensitive = TRUE) |stri_detect_fixed(inputseqaabys, v3v5combo, case_insensitive = TRUE))

(stri_detect_fixed(v2v3combo,inputseqaabys, case_insensitive = TRUE) |stri_detect_fixed(inputseqaabys, v2v3combo, case_insensitive = TRUE))

k = newCombinedNoGap[4231,1]
v= pairwiseAlignment(inputseqaabys,k)
writePairwiseAlignments(v)

overlapIsTrue <- function(combined, MScandidate) {
  t=pairwiseAlignment(combined, MScandidate)
  
  MScanindex = t@subject@range@start
  comboindex = t@pattern@range@start
  cutcombined = substr(combined,comboindex,nchar(combined))
  cutMScan = substr(MScandidate,1,nchar(combined)-comboindex+1)
  if (cutcombined ==cutMScan ) {
    return (TRUE)
  }
  return (FALSE)
}

overlapIsTrue(v3v5combo,inputseqaabys)
overlapIsTrue(inputseqaabys, v3v5combo)

combined = v3v5combo
MScandidate=inputseqaabys
c = substr(combined,comboindex,nchar(combined))
nchar(combined)-comboindex +1
nchar(combined)-28
