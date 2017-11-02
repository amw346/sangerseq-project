#10-17-2017 Final Version

library(seqinr)

#Generate dataframe of desired files [118 haplotype files]
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
for (i in 1:118) {
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
}

addseq <- function(seq1,seq2, index, CutSeq1) { 
  #input is two string sequences of DNA, an index integer where alligning should start and a boolean cutseq1 indicating which sequence to cut
  #cutseq1 = TRUE means seq1 should be cut
  #output is a combined sequence starting at index 
  
  short = seq1 
  long = seq2 
  if (CutSeq1) { 
    short = seq2 
    long = seq1 
  } 
  len = nchar(long) 
  #cut long string 
  newlong = substr(long,index,len) 
  
  lenlong = nchar(newlong) 
  lenshort = nchar(short) 
  length = lenlong 
  if (lenlong > lenshort) { 
    length = lenshort 
  } 
  
  #adding the elements 
  combined = "" 
  for (i in 1:length) { 
    newchar = addbases(substr(short,i,i), substr(long,i,i)) 
    combined= paste0(combined, newchar) 
  } 
  
  #returning the added sequence 
  return (combined) 
}

addbases <- function(a,b) {
  #addbases inputs two single character strings and outputs a single letter string according to IUPAC codes
  if (a == "a") {
    if (b == "a") {return ("A")}
    if (b == "c") {return ("M")}
    if (b == "g") {return ("R")}
    if (b == "t") {return ("W")}
  }
  
  if (a == "c") {
    if (b == "a") {return ("M")}
    if (b == "c") {return ("C")}
    if (b == "g") {return ("S")}
    if (b == "t") {return ("Y")}
  }
  
  if (a == "g") {
    if (b == "a") {return ("R")}
    if (b == "c") {return ("S")}
    if (b == "g") {return ("G")}
    if (b == "t") {return ("K")}
  }
  if (a == "t") {
    if (b == "a") {return ("W")}
    if (b == "c") {return ("Y")}
    if (b == "g") {return ("K")}
    if (b == "t") {return ("T")}
  }
  if (a == "N") {return("N")}
  return ("input strings not valid")	
}


#Code to generate the matrix of added pairwise combinations
combinedtypes = matrix() 

#two loops to make ((n-1)*n)/2 pairwise comparisons
for (d in 1:117) { 
  y = 118-d 
  for (i in 1:y) { 
    
    #initialize two sequences
    seq1 = seq[d,1] 
    seq2 = seq[i+d,1] 
    
    #alligning them
    allign = pairwiseAlignment(seq1,seq2)
    s = summary(allign) 
    
    #define start index
    x = "NaN" 
    index = 0 
    while (x == "NaN") { 
      index = 1 +index 
      x = s@mismatchSummary$pattern$position[index,3] 
    } 
    
    #define cutseq1 boolean
    #cutseq1 = TRUE means the first sequence inputted into pairwiseAllign needs to be cut
    cutSeq1 = TRUE
    if  (allign@pattern@range@start == 1) {
      cutSeq1 == FALSE
    }
    
    #combine seq1&2 and cut them to size
    combined =addseq(seq1,seq2, index, cutSeq1) 
    combinedtypes = rbind(combinedtypes,combined)
  } 
}



#Testcases 10182017

seqA =  "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttctgat"
seqB =  "gaatttcaccgacttcatgcacagcttcatgattgtgttcc"
seqC = "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactctctctttctctctctctcaactttattccgaccatcctttgcag"
seqD = "gtaagttgccaggatggaggcttctgat"
seqE = "ccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactc"
seqF = "gcacagcttcatgattgtgttcc"
seqG = "cccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactc"
seqH = "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttctgat"

seq = data.frame(c(1:8))
seq[1,1] = seqA
seq[2,1] = seqB
seq[3,1] = seqC
seq[4,1] = seqD
seq[5,1] = seqE
seq[6,1] = seqF
seq[7,1] = seqG
seq[8,1] = seqH

TESTcombinedtypes3 = matrix() 
#two loops to make ((n-1)*n)/2 pairwise comparisons
for (d in 1:7) { 
  y = 8-d 
  for (i in 1:y) { 
    
    #initialize two sequences
    seq1 = seq[d,1] 
    seq2 = seq[i+d,1] 
    
    #alligning them
    allign = pairwiseAlignment(seq1,seq2)
    s = summary(allign) 
    
    #define start index
    x = "NaN" 
    index = 0 
    while (x == "NaN") { 
      index = 1 +index 
      x = s@mismatchSummary$pattern$position[index,3] 
    } 
    
    #define cutseq1 boolean
    #cutseq1 = TRUE means the first sequence inputted into pairwiseAllign needs to be cut
    cutSeq1 = TRUE
    if  (allign@pattern@range@start == 1) {
      cutSeq1 == FALSE
    }
    
    #combine seq1&2 and cut them to size
    combined =addseq(seq1,seq2, index, cutSeq1) 
    TESTcombinedtypes3 = rbind(TESTcombinedtypes3,combined)
  } 
}


#10-19-2017 edits
ab = pairwiseAlignment(seqA,seqH)
ab = pairwiseAlignment(seqA,seqB)
writePairwiseAlignments(ab)

TESTcombinedtypes[2]
D = DNAString(TESTcombinedtypes[2])
pairwiseAlignment(TESTcombinedtypes[3],TESTcombinedtypes[4])


for (d in 1:7) { 
  y = 8-d 
  for (i in (d+1):8) { 
        #initialize two sequences
        seq1 = seq[d,1] 
        seq2 = seq[i+d,1] 
        
        #alligning them
        allign = pairwiseAlignment(seq1,seq2)
        s = summary(allign) 
        
        #define start index
        x = "NaN" 
        index = 0 
        while (x == "NaN") { 
          index = 1 +index 
          x = s@mismatchSummary$pattern$position[index,3] 
        } 
        
        #define cutseq1 boolean
        #cutseq1 = TRUE means the first sequence inputted into pairwiseAllign needs to be cut
        cutSeq1 = TRUE
        if  (allign@pattern@range@start == 1) {
          cutSeq1 == FALSE
        }
        
        #combine seq1&2 and cut them to size
        combined =addseq(seq1,seq2, index, cutSeq1) 
        combinedtypes = rbind(combinedtypes,combined)
      } 
    } 
    

#test for 8 sequenves
for (d in 1:7) { 
  y = 8-d 
  for (i in 1:y) { 
    c = c(d, i+d)
    t = i +d
    a = pairwiseAlignment(seq[d,1],seq[t,1])
    writePairwiseAlignments(a)
    print(c(c, seq[d,1],seq[t,1]))
  }}

seq[5,1]


# compare
getPriSec <- function(file) {
  #getPriSec inputs a file path ie "/Users/aliciawilliams/Desktop/ab1/practicefile2.ab1" and outputs a list containing the primary and secondary sequences as strings
  sangerobj <- readsangerseq(file)
  basecalls <- makeBaseCalls(sangerobj)
  primary <- primarySeq(basecalls, string = 'TRUE')
  secondary <- secondarySeq(basecalls, string = 'TRUE')
  return((list(PrimarySeq = primary, SecondarySeq = secondary)))
}

library(sangerseqR)
sangerobj <- readsangerseq("C:/Users/amw346/Desktop/aabys11May16-2FkdrFL-R7skdrFL.ab1")
getPriSec()

#original
clipToString <- function(sangob) { 
  pri = primarySeq(sangob,string = TRUE) 
  len = nchar(pri) 
  num = 0 
  found = FALSE 
  for (i in 1:len) { 
    c1= substr(pri,i,i) 
    c2 = substr(pri,i+1,i+1) 
    c3 = substr(pri,i+2,i+2) 
    if (found) { 
      return(substr(pri,20,num-2)) 
    } 
    if ((c1 == "N") & c2 == "N" & c3 == "N") { 
      found = TRUE 
    } 
    num = i +1 
  } 
  return (num)
} 


#modified for compareHAP
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
    if ((c1 == "N") & c2 == "N" & c3 == "N") { 
      found = TRUE 
    } 
    num = i +1 
  } 
  return (num-2)
} 

#10202017 more edits
compareHap<- function(file) {
  sangerobj <- readsangerseq(file)
  index = clipIndex(sangerobj)
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index)
  cutsec = substr(sec,20,index)
  combined = addPriSec2(cutpri,cutsec)
  match = checkMasterList(combined,combinedtypes)
}

checkMasterListORIGINAL<- function(seq, master) {
  len = dim(master)[1]
  found = FALSE
  
  for (i in 2:len) {
  
    master[i]
    if (found = TRUE) {
      print("match found")
      return(list[i])
    }
  }
  print("no match found") 
}
  
p = pairwiseAlignment(combined,seqA)
writePairwiseAlignments(p)

addPriSec2 <- function(pri,sec) {
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

addbases2 <- function(a,b) {
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

#20171024 with edits from 20171025
library(stringi)

#MOST UP TO DATE CHECKMASTERLIST
checkMasterList<- function(newseq, master) {
  len = dim(master)[1]
  found = FALSE
  
  for (i in 2:len) {
    found = (stri_detect_fixed(newseq,master[i,1], case_insensitive = TRUE) | stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE))
    if (found == TRUE) {
      print("match found")
      return(master[i,1])
    }
  }
  print("no match found") 
  return(found)
}

#RANDOM CRAP USED TO TEST CHECK MASTER LIST AND COMPAREHAP
newseq =  "gtaagttgacgtggccgaaactgctcccccgctcccaggatggag"

seq[1,2]

stri_detect_fixed(seq[1,1],newseq)
stri_detect_fixed("aaa", "AAA", case_insensitive = TRUE)
stri_detect_fixed("AAA", "A")



#NChis11May16-1FkdrFL-R7skdrFL

file="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"
  compareHap<- function(file, master) {
  sangerobj <- readsangerseq(file)
  index = clipIndex(sangerobj)
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index)
  cutsec = substr(sec,20,index)
  combined = addPriSec2(cutpri,cutsec)
  match = checkMasterList(combined,master)
  return (match)
  }
  
compareHap("C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1")

#20171026 Testing of comparehap
file="C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1"
compareHap<- function(file, master) {
  sangerobj <- readsangerseq(file)
  index = clipIndex(sangerobj)
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index)
  cutsec = substr(sec,20,index)
  combined = addPriSec2(cutpri,cutsec)
  match = checkMasterList(combined,master)
  return (match)
}
len
compareHap("C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1",TESTcombinedtypes)

TESTcombinedtypes


MODcheckMasterList<- function(newseq, master) {
  len = dim(TESTcombinedtypes3)[1]
  num= 0
  matchseq= matrix()
  for (i in 2:len) {
    found = FALSE
    found = (stri_detect_fixed(TESTcombinedtypes3[i,1],newseq4, case_insensitive = TRUE))
    if (found == TRUE) {
      print(i)
      print("match found")
      num = num +1
      print(num)
      matchseq= rbind(matchseq,TESTcombinedtypes3[i,1])
    }
  }
  
  return(c(match = matchseq, nummatches = num))
}

TESTcombinedtypes3[2,1]

matchseq= rbind(matchseq,TESTcombinedtypes3[2])
matchseq= rbind(matchseq,TESTcombinedtypes3[3])

matchseq= data.frame(matches = character())
matchseq= rbind(matchseq,c(matches = newseq))

matchseq= rbind(matchseq,TESTcombinedtypes3[3])
newseq =  "SYRARWYKRCKYSSCCGMWMCYRSKMYSSMGSYYYCASSMKKSRGRMTWMWKAT"
newseq2  = "SYRARWYKRCKYSSCCGMWMCYRSKMYSSMGSYYYCASSM"
newseq3= "SYRARWYKRCKYSSCCGMWMCYRSK"
newseq4 = "SYRARWYKSCWSSMYSGMKSCYWSKRWK"
newlist = MODcheckMasterList(newseq3, TESTcombinedtypes3)

stri_detect_fixed("A","AAA")
stri_detect_fixed("AAA","AA")

data = TESTcombinedtypes3[-1]
len= dim(TESTcombinedtypes3[1])[1]
for (i in 2:len) {
   MODcheckMasterList(TESTcombinedtypes3[i],TESTcombinedtypes3) 
}


str(TESTcombinedtypes3)

#10312017
seqA =  "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttctgat"
seqB =  "gaatttcaccgacttcatgcacagcttcatgattgtgttcc"
seqC = "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactctctctttctctctctctcaactttattccgaccatcctttgcag"
seqD = "gtaagttgccaggatggaggcttctgat"
seqE = "ccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactc"
seqF = "gcacagcttcatgattgtgttcc"
seqG = "cccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactc"
seqH = "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttctgat"

seq = data.frame(c(1:8))
seq[1,1] = seqA
seq[2,1] = seqB
seq[3,1] = seqC
seq[4,1] = seqD
seq[5,1] = seqE
seq[6,1] = seqF
seq[7,1] = seqG
seq[8,1] = seqH

TESTcombinedtypes4 = matrix() 
#two loops to make ((n-1)*n)/2 pairwise comparisons
for (d in 1:7) { 
  y = 8-d 
  for (i in 1:y) { 
    
    #initialize two sequences
    seq1 = seq[d,1] 
    seq2 = seq[i+d,1] 
    
    #alligning them
    allign = pairwiseAlignment(seq1,seq2)
    s = summary(allign) 
    
    #define start index
    x = "NaN" 
    index = 0 
    while (x == "NaN") { 
      index = 1 +index 
      x = s@mismatchSummary$pattern$position[index,3] 
    } 
    
    #define cutseq1 boolean
    #cutseq1 = TRUE means the first sequence inputted into pairwiseAllign needs to be cut
    cutSeq1 = TRUE
    if  (allign@pattern@range@start == 1) {
      cutSeq1 == FALSE
    }
    
    #combine seq1&2 and cut them to size
    combined =addseq(seq1,seq2, index, cutSeq1) 
    TESTcombinedtypes4 = rbind(TESTcombinedtypes4,combined)
  } 
}

#code to input your own 118
len5 = dim(seq)[1]
for (i in 1:len5) {
  print(i)
  seq1 = seq[i,1] 
  TESTcombinedtypes4 = rbind(TESTcombinedtypes4,seq1)
}

#how many am i inside
newseq4 = "SYRARWYKRCKYSSCCGMWMCYRSKMYSSMGSYYYCASSM" 
MODcheckMasterList2<- function(newseq, master,curI) {
  len = dim(master)[1]
  matchseq= matrix()
  matchnum= matrix()
  curnum = matrix()
  for (i in 2:len) {
    
    
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE))
    if (found == TRUE & i !=curI) {
      print(i)
      print("match found")
      matchseq= rbind(matchseq,master[i,1])
      matchnum= rbind(matchnum,i)
      #curnum = rbind(curnum,curI)
    }
  }
  newr= data.frame(seq=matchseq[,1],index = matchnum[,1])
  #, current = curnum[,1]
  if(nrow(newr) > 1) {
    return(newr)
  } else {
    empty = TRUE
    return(empty)
  }
 
}

#run through all the sequences in master
len6 = dim(TESTcombinedtypes4)[1] 
dataf = data.frame(seq= c("string"),index= c(1))
for (j in 2:len6) {
  newr = MODcheckMasterList2(TESTcombinedtypes4[j,1],TESTcombinedtypes4,j)
  if (!is.logical(newr)) {
    dataf = rbind(dataf,newr)
  } 
}


data = MODcheckMasterList2(newseq4, TESTcombinedtypes4, 1)


#11022017
#how many am i inside
newseq4 = "SYRARWYKRCKYSSCCGMWMCYRSKMYSSMGSYYYCASSM" 
MODcheckMasterList2<- function(newseq, master,curI) {
  len = dim(master)[1]
  matchseq= matrix()
  matchnum= matrix()
  curnum = matrix()
  for (i in 2:len) {
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE))
    if (found == TRUE & i !=curI) {
      print(i)
      print(curI)
      print("match found")
      matchseq= rbind(matchseq,master[i,1])
      matchnum= rbind(matchnum,i)
      curnum = rbind(curnum,curI)
    }
  }
  newr= data.frame(seq=matchseq[,1],index = matchnum[,1],original = curnum[,1])
  #, current = curnum[,1]
  if(nrow(newr) > 1) {
    return(newr)
  } else {
    empty = TRUE
    return(empty)
  }
  
}

#run through all the sequences in master
len6 = dim(TESTcombinedtypes4)[1] 
dataf = data.frame(seq= c("string"),index= c(1),original = c(1))
for (j in 2:len6) {
  newr = MODcheckMasterList2(TESTcombinedtypes4[j,1],TESTcombinedtypes4,j)
  if (!is.logical(newr)) {
    dataf = rbind(dataf,newr)
  } 
}

#code for combinedtypes
len7 = dim(combinedtypes)[1]
dataf = data.frame(seq= c("string"),index= c(1),original = c(1))
for (h in 2:len7) {
  newr = MODcheckMasterList2(combinedtypes[h,1],combinedtypes,h)
  if (!is.logical(newr)) {
    dataf = rbind(dataf,newr)
  } 
}

combinedtypes[3]
combinedtypes[56]
combinedtypes[300]
combinedtypes[4964]
