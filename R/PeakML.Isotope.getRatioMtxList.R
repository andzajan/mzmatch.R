PeakML.Isotope.getRatioMtxList <- function(intList, sampleGroups, useArea, metName=""){
	# PRE:
	#	intList <- list containing intensities of isotopes
	#	sample Groups <- the sample groupes
	#	useArea <-  TRUE will sum up all intensities to give area under the curve, FALSE gives maximum intensities
	#	metName <- name of the metabolite to give the column name of the matrix.
	# POST:
	#	Returns the matrix for each peak group, containing the ratio of labelled vs unlabelled signals for each replicate normalized to one along the columns and isotopes along the columns.
	
	getSignalRatioMtx <- function(abunList, sampleGroup, numCarbons, numReplicates){

		abunMtx <- matrix(nrow = numReplicates, ncol = numCarbons) # +1 is for adding the rownames
		rnames <- c(paste(sampleGroup, as.character(c(1:numReplicates)), sep=""))
		row.names(abunMtx) <- rnames
		
		for (itop in 2:numCarbons){
			for (repl in 1:numReplicates){
				if (!is.null(abunList[[1]][[repl]])){
					abunMtx[repl, 1] <- abunList[[1]][[repl]]
				} else{
					abunMtx[repl, 1] <- 0
				}
				if (!is.null(abunList[[itop]][[repl]])){
					abunMtx[repl, itop] <- abunList[[itop]][[repl]]# / abunList[[sampleGroup]][[1]][[repl]]
				} else{
					abunMtx[repl, itop] <- NA
				}
			}
		}
		abunMtx
	}

	ratioMtxList <- vector("list", length(intList))
	for (peakGroup in 1:length(intList)){
		abunList <- NULL
		ratioMtx <- NULL

		for (sampleGroup in 1:length(intList[[peakGroup]])){
			numIsotopes <- length(intList[[peakGroup]][[sampleGroup]])
			abunList[[sampleGroup]] <- vector("list", numIsotopes)
			
			for (isotop in 1:numIsotopes){
				numReplicates <- length(intList[[peakGroup]][[sampleGroup]][[isotop]])
				abunList[[sampleGroup]][[isotop]] <- vector("list", numReplicates) 
				
				for (replicate in 1:numReplicates){
					if(!is.null(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])){
						if (useArea==FALSE){
							abunList[[sampleGroup]][[isotop]][[replicate]] <- max(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])
						} else {
							abunList[[sampleGroup]][[isotop]][[replicate]] <- sum(intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]])
						}
					}
				}
			}
			ratioMtx <- rbind(ratioMtx, getSignalRatioMtx(abunList[[sampleGroup]], sampleGroups[[sampleGroup]], numIsotopes, numReplicates))
		}
		cNames <- paste(metName, as.character(c(0:(numIsotopes-1))))
		colnames(ratioMtx) <- cNames
		ratioMtx[is.na(ratioMtx)] <- 0
		ratioMtxList[[peakGroup]] <- t(ratioMtx)
	}
	ratioMtxList
}
