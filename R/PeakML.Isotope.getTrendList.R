PeakML.Isotope.getTrendList <- function(intList, sampleGroups, useArea){
	# PRE:
	#	intList <- List of intensities. See PeakML.Isotope.getChromData
	#	sampleGroups <- vector of sample groups
	#	useArea <- TRUE will sum up all intensities to give area under the curve, FALSE gives maximum intensities
	# POST:
	#	Mean of the sum or max of internsities of all replicates in the form trendList[[peakGroup]][[sampleGroup]][[isotop]] <- meanIntRep
	
	trendList <- vector("list", length(intList))
	
	for (peakGroup in 1:length(intList)){
		trendList[[peakGroup]]<-vector("list",  length(intList[[peakGroup]])) # This to store the values for the trend plot

		for (sampleGroup in 1:length(intList[[peakGroup]])){
			numIsotopes <- length(intList[[peakGroup]][[sampleGroup]])
			trendList[[peakGroup]][[sampleGroup]] <- vector("list", numIsotopes)

			for (isotop in 1:numIsotopes){
				numReplicates <- length(intList[[peakGroup]][[sampleGroup]][[isotop]])
				intensities <- rep(NA,numReplicates)

				for (replicate in 1:numReplicates){
					if(!is.null(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])){
						if (useArea==TRUE){
							intensities[[replicate]] <- sum(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])[1]
						} else {
							intensities[[replicate]] <- max(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])[1]
						}
					}
				}
				meanIntRep <- mean(intensities,na.rm=TRUE)
				if (!is.na(meanIntRep)){
					trendList[[peakGroup]][[sampleGroup]][[isotop]] <- meanIntRep
				}
			}
		}
	}
	trendList
}
