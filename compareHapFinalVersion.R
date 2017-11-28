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
  match = checkMasterList(combined,master)
  return (match)
}

compareHap("C:/Users/amw346/Desktop/NChis11May16-1FkdrFL-R7skdrFL.ab1",TESTcombinedtypes)
compareHap(file,)



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
    if ((c1 == "N") & c2 == "N" & c3 == "N") { 
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
  for (i in 1:6903) {
    found = FALSE
    found = (stri_detect_fixed(master[i,1],newseq, case_insensitive = TRUE) |stri_detect_fixed(newseq, master[i,1], case_insensitive = TRUE))
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