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
	 return (sang@primarySeq)
}

newPrimary <- function(sangob, string) {
# A function to input a sanger object and a string and reassign the primary sequence of a sangerobject to that string as a DNAString. Example code to put in command line: sang <- newPrimary(sang,sangstring) where sang is a sanger object and sang string is the desired new primary string.
	sangob@primarySeq<- DNAString(string)
	return(sangob)
}

fourpeaks <- function(sangob, ratio) {
	Apeaks = sangob@traceMatrix[,1]
	Cpeaks = sangob@traceMatrix[,2]
	Gpeaks = sangob@traceMatrix[,3]
	Tpeaks = sangob@traceMatrix[,4]
	len = dim(sangob@traceMatrix)[1]
	
	maxpeak = matrix(nrow = len, ncol = 1)
	for (i in 1:len) {
		maxpeak[i] =max(sangob@traceMatrix[i,])
	}
	
	Aratio = Apeaks/maxpeak
	Cratio = Cpeaks/maxpeak
	Gratio = Tpeaks/maxpeak
	Tratio = Tpeaks/maxpeak
	
	Acol = matrix(0L,nrow= len, ncol = 1)
	Ccol = matrix(0L, nrow= len, ncol = 1)
	Gcol = matrix(0L, nrow= len, ncol = 1)
	Tcol = matrix(0L, nrow= len, ncol = 1)
	
	primary = matrix(nrow= len)
	secondary= matrix(nrow= len) 
	warningposition = matrix()
	for (i in len) {
	 if (maxpeak[i] == 0) {
	 	primary[i] = "N"
	 	secondary [i] = "N"
	 }
	 
	 warning = 0
	 if (Aratio[i] > ratio) {
	 	Acol[i] = Aratio[i]
	 	warning = warning +1
	 }
	 
	 if (Cratio[i] > ratio) {
	 	Ccol[i] = Cratio[i]
	 	warning = warning +1
	 }
	 
	 if (Gratio[i] > ratio) {
	 	Gcol[i] = Gratio[i]
	 	warning = warning +1
	 }
	 
	 if (Tratio[i] > ratio) {
	 	Tcol[i] = Tratio[i]
	 	warning = warning +1
	 }
	 if (warning == 3) | (warning == 4) {
	  warningposition = cbind(warningposition,i)
	 }
	 ACGT = cbind(Aratio[i],Cratio[i],Gratio[i],Tratio[i])
	 max = maxcol(Aratio[i],Cratio[i],Gratio[i],Tratio[i])
	 if max == 1 {
	 	
	 }
	 
	  
	}
	print("More than two bases above threshold check following positions:")
	print(warningposition)
	return() 	
}

addDNACol <- function(Aratio,Cratio,Gratio,Tratio) {
	len = dim(Aratio)[1]
	for i in 
}


addBase <- function(nuc1,nuc2) {
	if 
}


test1 <- function(sangob, len) {
	maxpeak = matrix(nrow = len, ncol = 1)
	for (i in 1:10) {
		maxpeak[i] =max(sangob@traceMatrix[1,])
		print ( "There was error at ")
		print(maxpeak[i])
	}
	return(maxpeak[1:100])
}


cutpoint <- function(sangob) {
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
	 	num = 1 + i
	 	}
	 return(num-2) 
}


test7 <- function(sangob, ratio, start) {
	pri = primarySeq(sangob,string = TRUE)
	len = nchar(pri)
	cutpoint = cutpoint(sangob)
	hetcalls = makeBaseCalls(sangob, ratio)
	chromatogram(hetcalls, width = 75, height = 4, showcalls = "both", 
	              trim5 = start, trim3 =  len - cutpoint)
	newpri = substr(primarySeq(hetcalls, string = TRUE), start, cutpoint)
	sec = secondarySeq(hetcalls, string =TRUE)
    len = nchar(newpri)
    warningpos = matrix()
    for (i in 1:len) {
    	b = substr(sec,i, i)
    	if ( (b != "A") & (b != "C") & (b != "G") & (b != "T")) {
    		warningpos = cbind(warningpos,i)
    	}
    }
    if (dim(warningpos)[2] > 1 ) {
    	print("Warning: More than two peaks above threshold at following positions:")
    	return(warningpos)	
    }
    return ("No warnings")
}






