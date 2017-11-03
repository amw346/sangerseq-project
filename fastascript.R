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

#two loops to make ((n-1)*n)/2 pairwise comparisons
for (d in 1:117) { 
  y = 118-d 
  for (i in 1:y) { 
    d=1
    i=1
    #initialize two sequences
    seq1 = seq[d,3] 
    seq2 = seq[i+d,3] 
   
    #alligning them
    allign = pairwiseAlignment(seq1,seq2)
    writePairwiseAlignments(allign)
    s = summary(allign) 
    index = allign@pattern@range@start
    
    allign2 = pairwiseAlignment(seq2,seq1)
    writePairwiseAlignments(allign2)
    s = summary(allign2) 
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
    newCombined[i+d,1] = combined
  } 
}






