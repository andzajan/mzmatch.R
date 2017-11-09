PeakML.RankProducts <- function (sampleGroups, inputTable, groupA, groupB, numberOfPermutations=1000, calculateProduct=TRUE, MinNumOfValidPairs=NA, na.rm=FALSE, RandomPairs=NA)
{
	## Class labels
	HITS <- c(which(sampleGroups==groupA), which(sampleGroups==groupB))

	sampleGroups <- sampleGroups[c(HITS)]
	IntensitiesData <- t(inputTable[c(HITS),])
	cl <- as.numeric(as.factor(sampleGroups))-1

	out <- RankProducts (data=IntensitiesData, cl=cl, num.perm=1000, logged=FALSE, na.rm=na.rm, gene.names=NULL, plot=FALSE, rand = NULL, calculateProduct=calculateProduct, MinNumOfValidPairs=MinNumOfValidPairs, RandomPairs=RandomPairs)
	
	out$sampleindex <- HITS
	
	out
}
