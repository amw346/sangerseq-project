getPriSec <- function(file) {
#getPriSec inputs a file path ie "/Users/aliciawilliams/Desktop/ab1/practicefile2.ab1" and outputs a list containing the primary and secondary sequences as strings
	sangerobj <- readsangerseq(file)
	basecalls <- makeBaseCalls(sangerobj)
		primary <- primarySeq(basecalls, string = 'TRUE')
		secondary <- secondarySeq(basecalls, string = 'TRUE')
	return((list(PrimarySeq = primary, SecondarySeq = secondary)))
}

addbases <- function(a,b) {
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

addPriSec <- function(pri,sec) {
#addPriSeq inputs two strings of same length and outputs a string that is the combined version according to IUPAC codes

	 #checking if same length otherwise error message
	 len = nchar(pri)
	 len2 = nchar(sec)
	 if (len != len2) {return ("input sequences different lengths")}
	 
	 #adding the elements
	 combined = ""
	 for (i in 1:len) {
	 	newchar = addbases(substr(pri,i,i), substr(sec,i,i))
	 	combined= paste0(combined, newchar)
	 }
	 
	 #returning the added sequence
	 return (AddedSeq = combined)
}

library(seqinr)
file2 = read.fasta(file="Z:/Shared Documents/Alicia Williams/r/test1.fas", as.string = TRUE)
view(file)
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
namelist = data.frame(n)
name <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas")
for (i in 1:118) {
  namelist[i,1] = filenames[i]
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)
}

seq1 = seq[1,1]
seq2 = seq[2,1]

allign = pairwiseAllignment(seq1, seq2)

library(sangerseqR)
read.fasta("Z:/Shared Documents/Alicia Williams/haplotypes/kdr-his6 sans exons.fas")
readsangerseq("/Users/amw346/Desktop/test2.ab1")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

install.packages("sangerseqR", lib= "C:/Users/amw346/Documents/R/win-library/3.4")


library(seqinr)
file2 = read.fasta(file="Z:/Shared Documents/Alicia Williams/r/test1.fas", as.string = TRUE)
view(file)
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
namelist = data.frame(n)
name <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas")
for (i in 1:118) {
  namelist[i,1] = filenames[i]
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)
}
library(sangerseqR)
read.fasta("Z:/Shared Documents/Alicia Williams/haplotypes/kdr-his6 sans exons.fas")

seq2 = "cgacttcatgcacagcttcatgattgtgttccgagtgctgtgcggagagtggatcgagtccatgtgggactgcatgtatgtgggcgatgtcagctgtatacccttcttcttggccacggtcgtgataggcaatcttgtggtaagttgacgtggccgaaactactcccccgctcccaggatggaggcttctgatggccaattaaaaaagaacatttattttaaatcaacctctctctttctctctctctctcaactttattccgtccatcctttgcaggttcttaat"

seq1 = "gaatttcaccgacttcatgcacagcttcatgattgtgttccgagtgctgtgcggagagtggatcgagtccatgtgggactgcatgtatgttggcgatgtgagctgtatacccttcttcttggccacggtcgtgatcggcaatcatgtggtaagttgacgtggccgaaactactcccccgctcccaggagagaggcttgatccgtaatatacaaaaatttgacatttatctctctctttctctctcccaactttattctctccactgttgcaggttcttaat"
seq3 = "tcgagtccatgtgggactgcatgtatgttggcgatgtgagctgtatacccttcttcttggccacggtcgtgatcggcaatcatgtggtaagttgacgtggccgaaactactcccccgctcccaggagagaggcttgatccgtaatatacaaaaatttgacatttatctctctctttctctctcccaactttattctctccactgttgcaggttcttaat"
allign = pairwiseAlignment(seq1,seq2)
writePairwiseAlignments(allign)

combinedtypes = matrix[]
#TODo nest loops: for i 1:117, y = 118 -x , for i in 1:y
n = 118
for (i in 1:n) {
  x = FALSE
  seq1 = seq[i,1]
  seq2 = seq[i+1,1]

  allign = pairwiseAlignment(seq1,seq2)
  allign2 = pairwiseAlignment(seq1,seq3)
  writePairwiseAlignments(allign)
  writePairwiseAlignments(allign2)
  s = summary(allign2)
  s@position
  allign@start
  
  allign@ranges
  p = pattern(allign)
  str(allign)
  str(p)
  allign@pattern@unaligned@ranges
  x= mismatchTable(allign)
  bstring = aligned(allign)
  str(bstring)
  bstring@ranges
  des@start 
  start = allign@start
  start = compareStrings(allign)
  indel(allign)
  deletion(allign)
  nedit(allign)
  DNAString(seq1)
  DNAString(seq2)
  s = summary(allign)
  a  = summary(x)
  str(s)
  count = 1
  s@mismatchSummary$pattern$position[count,3]
  
  
  if (startseq1>1){
    x  = TRUE
  }
  
  combined =addseq(seq1,seq2, start, x)
  combinedtypes.append(combined)
}
s = 
index = 1
x = "NaN"
while (x = "NaN") {
  index = 1 +index
  s = summary(allign)
  x = s@mismatchSummary$pattern$position[count,3]
}



addseq <- function(seq1,seq2, start,x) {
  short = seq2
  long = seq1
  if (x) {
    short = seq2
    long = seq1
  }
  len = nchar(long)
  #cut long string
  newlong = substr(long,start,len)
  len2 = nchar(sec)
  
  #adding the elements
  #TODO: diff lengths issue
  combined = ""
  for (i in 1:len) {
    newchar = addbases(substr(short,i,i), substr(long,i,i))
    combined= paste0(combined, newchar)
  }
  
  #returning the added sequence
  return (combined)
}


addbases <- function(a,b) {
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




