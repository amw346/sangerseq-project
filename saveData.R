#save important dataframes
newCombinedNoGapAllAdded7021 = newCombinedNoGapAllAdded
write.csv(newCombinedNoGapAllAdded7021, "C:/Users/amw346/Desktop/newCombinedNoGapAllAdded7021.csv")
write.csv(seq, "C:/Users/amw346/Desktop/seq118.csv")

#contains index of pair in NewCombinedNoGap pair 1 information, index of pair in newCOmbinedNOGap, pair 2 info and then last 4 columns corespond to which pair its from, last col is counting if match is due to 3 pairs that are the same
write.csv(matchesCOMBONoGap2, "C:/Users/amw346/Desktop/matchesCOMBONoGap505.csv")

write.csv(matchesog, "C:/Users/amw346/Desktop/matchesog3.csv")

#code for saving the evolution of newCombined
write.csv(combinedtypes,"C:/Users/amw346/Desktop/evo/combinedtypes6903archive.csv")
write.csv(newCombined, "C:/Users/amw346/Desktop/evo/newCombined6903archive.csv")
write.csv(newCombinedNoGap, "C:/Users/amw346/Desktop/evo/newCombinedNoGap6903archive.csv")
write.csv(newCombinedNoGapAllAdded, "C:/Users/amw346/Desktop/evo/newCombinedNoGapAllAdded7021archive.csv")

write.csv(namelist,"C:/Users/amw346/Desktop/evo/namelist118.csv")

#investigation of matches
#dataframe that is all the sequences that contain v101 or v70
#and will therefore have a match with at least one other
#contains combined sequences, names, index in newCombinedNewGap
write.csv(samepairs,"C:/Users/amw346/Desktop/evo/samepairs10170_233.csv")
write.csv(samepairs1691,"C:/Users/amw346/Desktop/evo/samepairs1691_233.csv")

#dataframe of indexes in NewCombinedNewGap of both pairs that are matches
#original 587 matches 
write.csv(matchesCOMBO,"C:/Users/amw346/Desktop/evo/matchesCOMBOog587_archive.csv")
#reduced to 505 after dealing with gap / added beginning part
write.csv(matchesCOMBONoGap,"C:/Users/amw346/Desktop/evo/matchesCOMBONoGap505_archive.csv")