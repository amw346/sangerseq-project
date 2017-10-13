#Test 10132017 Generating Combined Haplotypes


library(seqinr)
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
namelist = data.frame(n)
name <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas")
for (i in 1:118) {
  namelist[i,1] = filenames[i]
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
}

seqA =  "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttctgat"
seqB = "gaatttcaccgacttcatgcacagcttcatgattgtgttcc"
seqC = "gtaagttgacgtggccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactctctctttctctctctctcaactttattccgaccatcctttgcag"
seqD = "gtaagttgccaggatggaggcttctgat"
seqE = "gtggccgaaactgctcccccgctcccaggatggaggcttcacccgtcggaataaatatattttacacttatcactc"
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
combinedtypes[4,1]

combinedtypes = data.frame() 
for (d in 1:7) { 
  y = 8-d 
  
  for (i in 1:y) { 
    
    seq1 = seq[d,1] 
    seq2 = seq[i+d,1] 
    
    allign = pairwiseAlignment(seq1,seq2) 
    s = summary(allign) 
    
    #define start index
    x = "NaN" 
    index = 0 
    while (x == "NaN") { 
      index = 1 +index 
      x = s@mismatchSummary$pattern$position[index,3] 
    } 
    
    cutSeq1 = TRUE
    #define cutseq1 boolean
    if  (allign@pattern@range@start == 1) {
      cutSeq1 == FALSE
    }
    
    combined =addseq(seq1,seq2, index, cutSeq1) 
    name = paste(d, "and", i) 
    combinedtypes = rbind(combinedtypes, c(name,combined))
  } 
}

addseq <- function(seq1,seq2, index, CutSeq1) { 
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
