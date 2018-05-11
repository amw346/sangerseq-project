#File contains the function compareHap() which compares a file to a database  
#also contains all the dependant functions called in compareHap()
#made in R version 3.4.2

#necessary packages: sangerseqR, Biostrings, BioCGenerics, parallel
library(stringi)

#examples of how to use function and arguments
#example file="C:/Users/amw346/Desktop/aabys11May16-3F kdrFL-R7 s kdrFL.ab1"
#example database: newCombinedNoGapAllAdded
compareHap(file, newCombinedNoGapAllAdded)

#warnings
#addbases2 requires capital letters
#you should run library(stringi)before hand otherwise stri_detect_fixed will throw error


compareHap<- function(file, master) {
  #inputs a chromatogram file and outputs the list of matches from database and the pairs of haplotypes that was a match
  sangerobj <- readsangerseq(file) #read in file
  index = clipIndex2(sangerobj) #finds the index to cut off all Ns at the end
  
  #cut and add combined string
  basecalls <- makeBaseCalls(sangerobj)
  pri <- primarySeq(basecalls, string = 'TRUE')
  sec <- secondarySeq(basecalls, string = 'TRUE')
  cutpri = substr(pri,20,index) #cuts string from 20 to the index
  cutsec = substr(sec,20,index)
  combined = addPriSeq3(cutpri,cutsec)
  print(combined)
  
  #look for matches within the master list supplied to function
  match = MODcheckMasterList4(combined,master)

  return (match)
}

clipIndex2 <- function(sangob) { 
  #input a sangerobject from sangerseqR package and will output the index of the first occurance of three NNNs.
  pri = primarySeq(sangob,string = TRUE) 
  len = nchar(pri) 
  num = 0 
  print(num)
  found = FALSE 
  for (i in 1:len) { 
    c1= substr(pri,i,i) 
    c2 = substr(pri,i+1,i+1) 
    c3 = substr(pri,i+2,i+2) 
    if (found) { 
      print(i)
      print("found was true")
      return(num-2) 
    } 
    if ((c1 == "N") & (c2 == "N") & (c3 == "N")) { 
      found = TRUE 
    } 
    num = i +1 
  } 
  print("theres definitely an error coming from clipindex2")
  return (num-2)
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


MODcheckMasterList4<- function(newseq, master) {
  #input is a sequence in the form of a string and output is dataframe of matches
  
  #make dataframe for the matches to return in
  matchesList = data.frame()
  len = dim(master)[1]
  count = 1
  #loop through all possible sequences in database
  for (i in 1:7021) { #modify this loop range to run through subset of database
    #reset varible to not found and then check if it matches
    print(i)
    found = FALSE
    overlap = overlapIsTrue(newseq,master[i,1])
    
    #looks for overlap or one inside the other
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE) |stri_detect_fixed(newseq, master[i,1], case_insensitive = TRUE) | overlap )
    
    if (found == TRUE) {
      matchesList[count,1]= master[i,1]
      matchesList[count,2]= master[i,2]
      matchesList[count,3]= master[i,3]
      print(i)
      print("match found")
      count = 1+count
    }
  }
  
  
  return(matchesList)
}


overlapIsTrue <- function(long, short) {
  #function inputs two strings and checks to see if they have a matching overlapping portion
  #returns true if match is found
  print(long)
  print(short)
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





sink("C:/Users/amw346/Desktop/testerrror.txt")
compareHap("C:/Users/amw346/Desktop/KS17 kdr sequences/KS2017_18Jul17-26M kdrFL-R7 s kdrFL.ab1", newCombinedNoGapAllAdded)
sink()
