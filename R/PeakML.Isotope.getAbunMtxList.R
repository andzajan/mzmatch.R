PeakML.Isotope.getAbunMtxList <- function(isotopeChroms, sampleGroups, useArea){
	# PRE:
	#	isotopeChroms <- the chrom data for isotopes see. PeakML.Isotope.getChromData
	#	sampleGroup <- the sample groups
	#	useArea -> TRUE will sum up all intensities to give area under the curve, FALSE gives maximum intensities
	# POST:
	#	Generates matrix containing the relative abundance of labelled vs unlabelled isotopes based on the the mean of the peaks in each groups

	getAbunMtx <- function(abunList, sampleGroups, numCarbons){
		cat("sampleGroups : ", sampleGroups, "\n")
		abunMtx <- matrix(nrow = length(sampleGroups), ncol = numCarbons)
		#row.names(abunMtx) <- sampleGroups
		names(abunMtx) <- c("AI", 1:(numCarbons-1))
	
		for (sg in 1:length(sampleGroups)){
			maxSignals <- unlist(abunList[[sg]][[1]])
			for (itop in 2:numCarbons){
				numReplicates <- length(abunList[[sg]][[itop]])
				sigAbunts <- vector()
				labIntesity <- vector()
				for (repl in 1:numReplicates){
					if (!is.null(abunList[[sg]][[itop]][[repl]])){
						if (!is.null(abunList[[sg]][[1]][[repl]])){
							sigAbunts[repl] <- abunList[[sg]][[itop]][[repl]] / abunList[[sg]][[1]][[repl]]
						} else{
							labIntesity <- c(labIntesity, abunList[[sg]][[itop]][[repl]]) # this is in case no unlabelled peaks were found
						}
						maxSignals <- c(maxSignals, abunList[[sg]][[itop]][[repl]])
					} else{
						sigAbunts[repl] <- NA
						maxSignals <- c(maxSignals, NA)
					}
				}
				meanAbunts <- mean(sigAbunts, na.rm=TRUE)
				if (!is.nan(meanAbunts)){
					abunMtx[sg, itop] <- round(meanAbunts, 2)
				}
				if (length(labIntesity)==numReplicates){
					abunMtx[sg, itop] <- round(mean(labIntesity), 2)
				}
			}
			abunMtx[sg, 1] <- round(mean(maxSignals, na.rm=TRUE), 2)
		}
		abunMtx
	}
	
	intList <- isotopeChroms[[2]]
	abunMtxList <- vector("list", length(intList))
	
	for (peakGroup in 1:length(intList)){

		abunList <- vector("list", length(intList[[peakGroup]]))

		for (sampleGroup in 1:length(intList[[peakGroup]])){

			numIsotopes <- length(intList[[peakGroup]][[sampleGroup]])
			abunList[[sampleGroup]] <- vector("list", numIsotopes)

			for (isotop in 1:numIsotopes){
				numReplicates <- length(intList[[peakGroup]][[sampleGroup]][[isotop]])
				maxintensities <- rep(NA,numReplicates)
				sumintensities <- rep(NA, numReplicates)
				abunList[[sampleGroup]][[isotop]] <- vector("list", numReplicates) 
				
				for (replicate in 1:numReplicates){
					
					if(!is.null(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])){
						maxintensities[[replicate]] <- max(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])[1]
						if (useArea==FALSE){
							abunList[[sampleGroup]][[isotop]][[replicate]] <- max(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])
						} else {
							abunList[[sampleGroup]][[isotop]][[replicate]] <- sum(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])
							sumintensities[[replicate]] <- sum(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])[1]
						}
					}
				}
				if (useArea==FALSE){
					meanIntRep <- mean(maxintensities,na.rm=TRUE)
				} else{
					meanIntRep <- mean(sumintensities,na.rm=TRUE)
				}
			}
		}
		print(abunList)
		abunMtxList[[peakGroup]] <- getAbunMtx(abunList, sampleGroups, numIsotopes)
	}
	abunMtxList
}
