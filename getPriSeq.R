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
	for (i in len) {
		maxpeak[i] = max(sang@traceMatrix[1,])
	}
	
	Aratio = Apeaks/maxpeak
	Cratio = Cpeaks/maxpeak
	Gratio = Tpeaks/maxpeak
	Tratio = Tpeaks/maxpeak
	
	Acol = matrix(nrow= len, ncol = 1)
	for (i in len) {
	 if (Aratio[i] > ratio) {
	 	Acol[i] = "A"
	 }
	 if (Cratio[i] > ratio) {
	 	Ccol[i] = "C"
	 }
	 if (Gratio[i] > ratio) {
	 	Gcol[i] = "G"
	 }
	 if (Tratio[i] > ratio) {
	 	Tcol[i] = "T"
	 }
	 
	}
	SeqString = addDNACol(Aratio, Cratio, Gratio, Tratio)
	return( SeqString)
	
}

addDNACol <- function(Aratio,Cratio,Gratio,Tratio) {
	len = dim(Aratio)[1]
	for i in 
}


addBase <- function(nuc1,nuc2) {
	if 
}
