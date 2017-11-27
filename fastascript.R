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
      matchseq= rbind(matchseq, master[i,1])
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
newCombined[10,1]
newCombined[11,1]
newCombined[21,1]
newCombined[41,1]
newCombined[52,1]
newCombined[96,1]

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

#how many in 118 are exactly the same
library(stringi)
count = 1
matchesog = data.frame()
for (d in 1:117) { 
  y = 118-d 
  for (i in 1:y) { 
    print("ggg")
    print(i+d)
    print(d)
    #initialize two sequences
    seq1 = seq[d,3]
    name = seq[d,2]
    seq2 = seq[i+d,3] 
    name2 = seq[i+d,2]
    if (stri_detect_fixed(seq1,seq2, case_insensitive = TRUE) & stri_detect_fixed(seq2,seq1, case_insensitive = TRUE)) {
      matchesog[count, 1] = name
      matchesog[count,2] = name2
      count = count+1
    }
   
  }}


#how many matches match because of concidence gives index in newCombined
count = 1
matchesCOMBO = data.frame()
for (d in 1:6903) { 
  y = 6903-d 
  for (i in 1:y) { 
  
    #initialize two sequences
    seq1 = newCombined[d,1]
    name = d 
    seq2 = newCombined[i+d,1] 
    name2 = i+d
    if (stri_detect_fixed(seq1,seq2, case_insensitive = TRUE) & stri_detect_fixed(seq2,seq1, case_insensitive = TRUE)) {
      matchesCOMBO[count, 1] = name
      matchesCOMBO[count,2] = name2
      count = count+1
      print(d)
    }
    
  }}
seq1 = "CATGTGGTAAGTTGACGTGGCCGAAACTGCTCCCCCGCTCCCAGGATGGAGGCTTCTGATGRCCAATTAAAAAAAATTAAATCAACCTCTCTCTTTCTCTCTCTCTCWMMWYTWTWYYSYSYMYMYSYTKYRS"
seq2="CATGTGGTAAGTTGACGTGGCCGAAACTRCTCCCCCGCTCCCAGGAKRGAGGCTTSWKMYGKMMWATWMAAAAAWWTKAMATYWAYCTCTCTCTTTCTCTCTCYCWMMWYTWTWYTCYSYMCWKYYGYWGSWK"
d = 3
i = 119-d


newCombined[3148,1]
newCombined[5458,1]

newCombined[3148,2:3]
newCombined[5458,1:3]


newCombined[10,1:3]
newCombined[96,1:3]

#test what the heck is happening with the 587 matches
q = seq[11,3]
z = seq[97,3]
h = pairwiseAlignment(q,z)
writePairwiseAlignments(h)

s = seq[1,3]

newcombined[10]

g = pairwiseAlignment(q,s)
writePairwiseAlignments(g)

k = pairwiseAlignment(z,s)
writePairwiseAlignments(k)
h = pairwiseAlignment(seq[64,3],seq[104,3])
writePairwiseAlignments(h)
h@pattern@range@start
seq[64,3]
seq[104,3]
g = pairwiseAlignment(seq[104,3],seq[31,3])
writePairwiseAlignments(g)
seq4 = substr(seq[64,3],143,nchar(seq[64,3]))

g = pairwiseAlignment(seq[104,3],seq[31,3])
writePairwiseAlignments(g
        )

#11142017 tests
newCombined[10,1:3]
newCombined[11,1:3]

a = seq[11,3]
 b =seq[1,3]
c =seq[12,3]

g = pairwiseAlignment(a,b)
writePairwiseAlignments(g)

k = pairwiseAlignment(b,c, type="local")
pairwiseAlignment(b,c, type="local")
writePairwiseAlignments(k)

j = pairwiseAlignment(a,c)
writePairwiseAlignments(j)

samepairs = data.frame()
count = 1
for (i in 1:6903) {
  if (newCombined[i,2] == "v101_sans_exons.seq" | newCombined[i,3] == "v101_sans_exons.seq" | newCombined[i,3] == "v70_sans_exons.seq" | newCombined[i,2] == "v70_sans_exons.seq") {
    samepairs[count,1] = newCombined[i,1]
    samepairs[count,2] = newCombined[i,2]
    samepairs[count,3] = newCombined[i,3]
    samepairs[count,4] = i
    count = count +1
  }
}

for (i in 1:6903) {
  if ((newCombined[i,2] == "v101_sans_exons.seq" & newCombined[i,3] == "v70_sans_exons.seq") | (newCombined[i,3] == "v101_sans_exons.seq")) {
    indexpair = i
  }
}

for (i in 1:6903) {
  if ((newCombined[i,2] == "v16_sans_exons.seq" & newCombined[i,3] == "v91_sans_exons.seq") | (newCombined[i,3] == "v16_sans_exons.seq")) {
    indexpair = i
  }
}

for (i in 1:6903) {
  if ((newCombined[i,2] == "v45_sans_exons.seq" & newCombined[i,3] == "v96_sans_exons.seq") | (newCombined[i,3] == "v45_sans_exons.seq")) {
    indexpair = i
  }
}


for (j in 1:587) {
  if (matchesCOMBO[j,1] == 5189 | matchesCOMBO[j,2] == 5189) {
    found = TRUE
    nindex = j
  }
}

#add 1 if matchecombo is result of v101 v 71 match
for (v in 1:233) {
  index = samepairs[v,4]
  for (c in 1:587) {
    if (matchesCOMBO[c,1] == index | matchesCOMBO[c,2] == index) {
      matchesCOMBO2[c,3]= 1
    }
  }
} 

#pair dataframe for v16 and v91
samepairs1691 = data.frame()
count = 1
for (i in 1:6903) {
  if (newCombined[i,2] == "v16_sans_exons.seq" | newCombined[i,3] == "v16_sans_exons.seq" | newCombined[i,3] == "v91_sans_exons.seq" | newCombined[i,2] == "v91_sans_exons.seq") {
    samepairs1691[count,1] = newCombined[i,1]
    samepairs1691[count,2] = newCombined[i,2]
    samepairs1691[count,3] = newCombined[i,3]
    samepairs1691[count,4] = i
    count = count +1
  }
}

for (v in 1:233) {
  index = samepairs1691[v,4]
  for (c in 1:587) {
    if (matchesCOMBO[c,1] == index | matchesCOMBO[c,2] == index) {
      matchesCOMBO2[c,3]= 2
    }
  }
} 


#pair dataframe for v45 and v96
samepairs4596 = data.frame()
count = 1
for (i in 1:6903) {
  if (newCombined[i,2] == "v45_sans_exons.seq" | newCombined[i,3] == "v45_sans_exons.seq" | newCombined[i,3] == "v96_sans_exons.seq" | newCombined[i,2] == "v96_sans_exons.seq") {
    samepairs4596[count,1] = newCombined[i,1]
    samepairs4596[count,2] = newCombined[i,2]
    samepairs4596[count,3] = newCombined[i,3]
    samepairs4596[count,4] = i
    count = count +1
  }
}

for (v in 1:233) {
  index = samepairs1691[v,4]
  for (c in 1:587) {
    if (matchesCOMBO[c,1] == index | matchesCOMBO[c,2] == index) {
      matchesCOMBO2[c,3]= 3
    }
  }
}

count= 1
problems = data.frame()
for (e in (1:587)) {
  if (is.na(matchesCOMBO2[e,3])) {
    problems[count,1] = matchesCOMBO2[e,1]
    problems[count,2] = matchesCOMBO2[e,2]
    count = 1+ count
  }
}
newCombined[10,1:3]
newCombined[11,1:3]

#kdr2
a = seq[11,3]
#his1
b =seq[1,3]

#kdr3
c =seq[12,3]

g = pairwiseAlignment(a,b)
writePairwiseAlignments(g)


k = pairwiseAlignment(a,b, type="local")
writePairwiseAlignments(k)
str(k)
k@pattern@range@start
k@subject@range@start

k = pairwiseAlignment(b,a, type="local")
writePairwiseAlignments(k)
str(k)
k@pattern@range@start
k@subject@range@start
newCombined[18,2:3]


for (x in 1:329) {
pair1 = problems[x,1]
seq11= 
seq12 = 
pair2 = problems[x,2]
k = pairwiseAlignment(a,b, type="local")  
longstartindex = k@pattern@range@start
shortstartindex = k@subject@range@start
if (min(k@pattern@range@start,k@subject@range@start) != 1) {
    longstartindex = k@subject@range@start - (k@pattern@range@start-1) 
    if (k@pattern@range@start > k@subject@range@start) {
      longstartindex = k@pattern@range@start - (k@subject@range@start-1) 
    }
}
}

#REDO with no aliignment gap problem

#kdr2
a = seq[11,3]
#his1
b =seq[1,3]

#kdr3
c =seq[12,3]

g = pairwiseAlignment(a,b)
writePairwiseAlignments(g)


k = pairwiseAlignment(a,b, type="local")
writePairwiseAlignments(k)


#Code to generate the matrix of added pairwise combinations with allignment problem fixed
newCombinedNoGaptest = data.frame() 
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
    allign = pairwiseAlignment(seq1,seq2,type = "local")
    writePairwiseAlignments(allign)
    #get two start indexes
    indexseq1 = allign@pattern@range@start
    indexseq2 = allign@subject@range@start
  
    
    #in event neither seq starts from 1, find new start index for the long strand  
    if  (min(indexseq1, indexseq2) !=1) {
      #assume seq1 is longer and find start index for seq1
      seq2islong = FALSE
      longstartindex = indexseq1 - (indexseq2-1) 
      #check if seq2 is actually longer and change index if necessary
      if (indexseq2 > indexseq1) {
        seq2islong = TRUE
        longstartindex = indexseq2 - (indexseq1-1) 
      }
      #cut beginning of sequences
      if (seq2islong) {
        len = nchar(seq2)
        seq2cut = substr(seq2,longstartindex,len)
        seq1cut = seq1
      } else {
        len = nchar(seq1)
        seq1cut = substr(seq1,longstartindex,len)
        seq2cut = seq2
        
  #otherwise either 1 starts from 1 or they both do
   }} else {
     #if seq2 is longer cut beginning of seq2
     if (indexseq2 > indexseq1) {
       len = nchar(seq2)
       seq2cut = substr(seq2,indexseq2,len) 
       seq1cut= seq1
     #if seq1 is longer or both starting at one
     } else {
       len = nchar(seq1)
       seq1cut = substr(seq1,indexseq1,len) 
       seq2cut= seq2
     }
   }
  
    len1 = nchar(seq1cut)
    len2= nchar(seq2cut)
    length = min(len1,len2)
    
    combo = "" 
    for (i in 1:length) { 
      newchar = addbases(substr(seq1cut,i,i), substr(seq2cut,i,i)) 
      combo= paste0(combo, newchar) 
    } 
    
    newCombinedNoGaptest[count,1] = combo
    newCombinedNoGaptest[count,2] = name
    newCombinedNoGaptest[count,3] = name2
    newCombinedNoGaptest[count,4] = seq1
    newCombinedNoGaptest[count,5] = seq2
    count = count +1
  } 
}

#Code to generate dataframe of exact matches
#how many matches match because of concidence gives index in newCombined
library(stringi)
count = 1
matchesCOMBONoGap = data.frame()
for (d in 1:6903) { 
  y = 6903-d 
  for (i in 1:y) { 
    
    #initialize two sequences
    seqcombo1 = newCombinedNoGap[d,1]
    name = d 
    seqcombo2 = newCombinedNoGap[i+d,1] 
    name2 = i+d
    if (stri_detect_fixed(seqcombo1,seqcombo2, case_insensitive = TRUE) & stri_detect_fixed(seqcombo2,seqcombo1, case_insensitive = TRUE)) {
      matchesCOMBONoGap[count, 1] = name
      matchesCOMBONoGap[count, 2] = newCombinedNoGap[d,2]
      matchesCOMBONoGap[count, 3] = newCombinedNoGap[d,3]
      matchesCOMBONoGap[count, 4] = newCombinedNoGap[d,4]
      matchesCOMBONoGap[count, 5] = newCombinedNoGap[d,5]
      matchesCOMBONoGap[count,6] = name2
      matchesCOMBONoGap[count, 7] = newCombinedNoGap[i+d,2]
      matchesCOMBONoGap[count, 8] = newCombinedNoGap[i+d,3]
      matchesCOMBONoGap[count, 9] = newCombinedNoGap[i+d,4]
      matchesCOMBONoGap[count, 10] = newCombinedNoGap[i+d,5]
      count = count+1
      print(d)
    }
    
  }}

#testing above two
f = seq[1,3]
g = seq[11,3]

v = pairwiseAlignment(f,g,type = "local")
writePairwiseAlignments(v)
newCombinedNoGap[10,1]
substr(g,143,270)

f = seq[1,3]
d = seq[15,3]
j = pairwiseAlignment(f,d,type = "local")
writePairwiseAlignments(j)

#code to locate which ones are matching as a result of adding identical sequence pairs (ie v101 and v70) to the same sequence
#if you rerun this code because workspace is lost just change the 14 to 11 there is no overlap bewtten matches, 116 exist for all three
matchesCOMBONoGap2 = matchesCOMBONoGap

#pairs for v45 and v96
count = 1
for (i in 1:505) {
  sit1 = (matchesCOMBONoGap[i,2] == "v45_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v96_sans_exons.seq")
  sit2 = (matchesCOMBONoGap[i,2] == "v45_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v96_sans_exons.seq")
  sit3 = (matchesCOMBONoGap[i,3] == "v45_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v96_sans_exons.seq")
  sit4 = (matchesCOMBONoGap[i,3] == "v45_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v96_sans_exons.seq")
  
  sit5 = (matchesCOMBONoGap[i,2] == "v96_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v45_sans_exons.seq")
  sit6 = (matchesCOMBONoGap[i,2] == "v96_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v45_sans_exons.seq")
  sit7 = (matchesCOMBONoGap[i,3] == "v96_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v45_sans_exons.seq")
  sit8 = (matchesCOMBONoGap[i,3] == "v96_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v45_sans_exons.seq")
  if (sit1 | sit2 | sit3 | sit4 | sit5 | sit6 | sit7 | sit8) {
    matchesCOMBONoGap2[i,14] = 1
    count = count +1
  }
}

#pairs for v101 and v70
count2 = 0
for (i in 1:505) {
  sit1 = (matchesCOMBONoGap[i,2] == "v101_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v70_sans_exons.seq")
  sit2 = (matchesCOMBONoGap[i,2] == "v101_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v70_sans_exons.seq")
  sit3 = (matchesCOMBONoGap[i,3] == "v101_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v70_sans_exons.seq")
  sit4 = (matchesCOMBONoGap[i,3] == "v101_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v70_sans_exons.seq")
  
  sit5 = (matchesCOMBONoGap[i,2] == "v70_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v101_sans_exons.seq")
  sit6 = (matchesCOMBONoGap[i,2] == "v70_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v101_sans_exons.seq")
  sit7 = (matchesCOMBONoGap[i,3] == "v70_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v101_sans_exons.seq")
  sit8 = (matchesCOMBONoGap[i,3] == "v70_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v101_sans_exons.seq")
  if (sit1 | sit2 | sit3 | sit4 | sit5 | sit6 | sit7 | sit8) {
    matchesCOMBONoGap2[i,14] = 1
    count2 = count2 +1
  }
}

#pairs for v16 and v91
count3 = 0
for (i in 1:505) {
  sit1 = (matchesCOMBONoGap[i,2] == "v16_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v91_sans_exons.seq")
  sit2 = (matchesCOMBONoGap[i,2] == "v16_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v91_sans_exons.seq")
  sit3 = (matchesCOMBONoGap[i,3] == "v16_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v91_sans_exons.seq")
  sit4 = (matchesCOMBONoGap[i,3] == "v16_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v91_sans_exons.seq")
  
  sit5 = (matchesCOMBONoGap[i,2] == "v91_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v16_sans_exons.seq")
  sit6 = (matchesCOMBONoGap[i,2] == "v91_sans_exons.seq") & (matchesCOMBONoGap[i,3] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v16_sans_exons.seq")
  sit7 = (matchesCOMBONoGap[i,3] == "v91_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,7]) & (matchesCOMBONoGap[i,8] == "v16_sans_exons.seq")
  sit8 = (matchesCOMBONoGap[i,3] == "v91_sans_exons.seq") & (matchesCOMBONoGap[i,2] == matchesCOMBONoGap[i,8]) & (matchesCOMBONoGap[i,7] == "v16_sans_exons.seq")
  if (sit1 | sit2 | sit3 | sit4 | sit5 | sit6 | sit7 | sit8) {
    matchesCOMBONoGap2[i,14] = 1
    count3 = count3 +1
  }
}

#leftover problem dataset
#problems2 contains 
problems2 = data.frame()
count4 = 1
for (i in 1:505) {
  if (is.na(matchesCOMBONoGap2[i,14])) {
    problems2[count4,1]=  matchesCOMBONoGap2[i,2]
    problems2[count4,2]=  matchesCOMBONoGap2[i,3]
    problems2[count4,3]=  matchesCOMBONoGap2[i,7]
    problems2[count4,4]=  matchesCOMBONoGap2[i,8]
    problems2[count4,5]= i
    count4 = count4+1
  }
}

#code to look at what is happening with problems2
problemCompare <- function(index) {
  #index is for matchesCOMBoNoGap2 problems
  seq1 = matchesCOMBONoGap2[index,4]
  seq2 = matchesCOMBONoGap2[index,5]
  
  newcombinedindex = matchesCOMBONoGap2[index,1]
  
  print(seq1)
  print(seq2)
  allign1 = pairwiseAlignment(seq1,seq2,type = "local")
 
  print("second match")
  seq3 = matchesCOMBONoGap2[index,9]
  seq4 = matchesCOMBONoGap2[index,10]
  
  print(seq3)
  print(seq4)
  allign2 = pairwiseAlignment(seq3,seq4,type = "local")
  
  writePairwiseAlignments(allign1)
  writePairwiseAlignments(allign2)
  print(newCombinedNoGap[newcombinedindex,1])
  
}

#usign the function to test
problemCompare(1)
c =pairwiseAlignment(seq2,seq4)
writePairwiseAlignments(c)
# index =1 example of different seq added to same seq but starting after difference occurs

problemCompare(4)
# index = 4 example of adding a and g in one and g and a in the other so is 5
r= matchesCOMBONoGap2[4,4]
t = matchesCOMBONoGap2[4,9]

m= matchesCOMBONoGap2[4,5]
n = matchesCOMBONoGap2[4,10]

g=pairwiseAlignment(r,m, type="local")
writePairwiseAlignments(g)

problemCompare(7)
