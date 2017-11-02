library(seqinr)

#Generate dataframe of desired files [118 haplotype files]
n = c(1:118)
seq= data.frame(n)
filenames <- list.files("Z:/Shared Documents/Alicia Williams/haplotypes", pattern= "*.fas", full.name = TRUE)
for (i in 1:118) {
  seq[i,1] = read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
  seq[i,2] = names(read.fasta(file=filenames[i], as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE))
  }


seq[1,2]
current.seq = seq[1,1]
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


  paste("")  

write.fasta
read.fasta()






