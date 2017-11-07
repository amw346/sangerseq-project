library(seqinr)

#Generate dataframe of desired files [118 haplotype files]
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
for (i in 1:118) {
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
  seq[i,2] = names(read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE))
  }


library(stringi)
for (i in 1:118) {
  current.seq = seq[i,1]
  his = FALSE
  v= FALSE
if(stri_detect_fixed(substr(current.seq,1,7),"gtaagtt", case_insensitive = TRUE)) {
#if(charmatch("gtaagtt", substr(current.seq,1,7))==1) {
  name = seq[i,2]
  if (substr(name,5,7) == "his") {
    his = TRUE
    modseq = paste("catgtg", current.seq, sep="")
    seq[i,3] = modseq
    stop= nchar(name)-4
    newfilename = paste(substr(name,1,stop),".fasta", sep = "2")
    write.fasta(sequences = modseq, file.out = newfilename, names = name)
  }
  if (substr(name,1,1) == "v") {
    v = TRUE
    modseq = paste("cttgtg", current.seq, sep="")
    seq[i,3] = modseq
    stop= nchar(name)-4
    newfilename = paste(substr(name,1,stop),".fasta", sep = "2")
    write.fasta(sequences = modseq, file.out = newfilename, names = name)
  }
  #if (substr(name,1,1) != "v" & substr(name,5,7) != "his") {
  if (!his & !v) {
    modseq = paste("tttgtg", current.seq, sep="")
    seq[i,3] = modseq
    stop= nchar(name)-4
    newfilename = paste(substr(name,1,stop),".fasta", sep = "7")
    write.fasta(sequences = modseq, file.out = newfilename, names = name)
  }
}
}

for (d in 1:118) {
  if (is.na(seq[d,3])) {
    seq[d,3] = seq[d,1]}
}

#Code to generate the matrix of added pairwise combinations
newCombined = data.frame() 
count = 1
#two loops to make ((n-1)*n)/2 pairwise comparisons
for (d in 1:117) { 
  y = 118-d 
  for (i in 1:y) { 
   
    #initialize two sequences
    seq1 = seq[d,3]
    name = seq[d,2]
    seq2 = seq[i+d,3] 
    name2 = seq[i+d,2]
    
    #alligning them
    allign = pairwiseAlignment(seq1,seq2)
    
 
    index = allign@pattern@range@start
    
    allign2 = pairwiseAlignment(seq2,seq1)
   

    index2 = allign2@pattern@range@start
    
    #define cutseq1 boolean
    #cutseq1 = TRUE means the first sequence inputted into pairwiseAllign needs to be cut
    cutSeq1 = TRUE
    if  (allign@pattern@range@start == 1) {
      cutSeq1 = FALSE
      index = index2
    }
    
    #combine seq1&2 and cut them to size
    combined =addseq(seq1,seq2, index, cutSeq1) 
    newCombined[count,1] = combined
    newCombined[count,2] = name
    newCombined[count,3] = name2

    count = count +1
  } 
}

#CORRECTED ADDSEQ
addseq <- function(seq1,seq2, index, cutSeq1) { 
  #input is two string sequences of DNA, an index integer where alligning should start and a boolean cutseq1 indicating which sequence to cut
  #cutseq1 = TRUE means seq1 should be cut
  #output is a combined sequence starting at index 
  
  short = seq1 
  long = seq2 
  if (cutSeq1) { 
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
    newchar = addbases(substr(short,i,i), substr(newlong,i,i)) 
    combined= paste0(combined, newchar) 
  } 
  
  #returning the added sequence 
  return (combined) 
}

library(stringi)
newseq = newCombined[1,1]
master = newCombined
#comparing against itself
MODcheckMasterList2<- function(newseq, master,curI) {
  len = dim(master)[1]
  matchseq= matrix()
  matchnum= matrix()
  curnum = matrix()
  for (i in 1:len) {
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
len8 = dim(newCombined)[1] 
dataf2 = data.frame(seq= c("string"),index= c(1),original = c(1))
for (j in 1:len8) {
  newr = MODcheckMasterList2(newCombined[j,1],newCombined,j)
  if (!is.logical(newr)) {
    dataf2 = rbind(dataf2,newr)
  } 
}



#random code used to test
seq[1,1]
seq[1,2]
seq1 = seq[1,3]
seq2 =seq[2,3]
newCombined[1,1]

allign = pairwiseAlignment(seq[1,3],seq[2,3])
writePairwiseAlignments(allign)



allign2 = pairwiseAlignment(seq[2,3],seq[1,3])
writePairwiseAlignments(allign2)
index2 = allign2@pattern@range@start
combined =addseq(seq[1,3],seq[2,3], index, cutSeq1) 

newchar = addbases(substr(short,1,1), substr(newlong,1,1)) 
combined= paste0(combined, newchar)




