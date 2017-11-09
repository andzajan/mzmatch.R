PeakML.Isotope.getChromData <- function(isotopeList, chromDataList, phenoData, sampleGroups){
	# PRE:
	#	isotopeList -> List of isotopes see PeakML.Isotope.getIsotopes
	#	chromDataList -> the chrom data list see PeakML.Methods.getChromData
	# 	phenoData -> PeakML.Methods.getPhenoData
	#	sampleGroups -> the sample groups

	# POST:
	#	Generates a list of lists containing mass, intensities and retentiontimes of all the isotopes in the isotopeList
	
	getBestPeakID <- function(chromDataList, peakIDs){
		# This function is used if in case two or more peaks are detected of a given isotop mass window, takes the one with the highest intensity.
		# May need to think about this at some point.
		cdList <- c()
		for (peak in 1:length(peakIDs)){
			cdList[peak] <- max(chromDataList[[peakIDs[peak]]][2,])
		}
		peakIDs[which(cdList==max(cdList))[1]]
	}
	
	mzList <- vector("list", length(isotopeList))
	intList <- vector("list", length(isotopeList))
	rtList <- vector("list", length(isotopeList))
	
	for (peakGroup in 1:length(isotopeList)){
		mzList[[peakGroup]] <-  vector("list", length(sampleGroups))
		intList[[peakGroup]] <-  vector("list", length(sampleGroups))
		rtList[[peakGroup]] <-  vector("list", length(sampleGroups))

		for (sampleGroup in 1:length(sampleGroups)){
			# create empty list store the intensities, rt and masses to get the window for plotting
			numIsotopes <- length(isotopeList[[peakGroup]])
			mzList[[peakGroup]][[sampleGroup]] <- vector("list", numIsotopes)
			intList[[peakGroup]][[sampleGroup]] <- vector("list", numIsotopes)
			rtList[[peakGroup]][[sampleGroup]]<- vector("list", numIsotopes)

			replicatesList <- which(phenoData == sampleGroups[[sampleGroup]])
			numReplicates <- length(replicatesList)
			for (itop in 1:numIsotopes){
				mzList[[peakGroup]][[sampleGroup]][[itop]] <- vector("list", numReplicates)
				intList[[peakGroup]][[sampleGroup]][[itop]] <- vector("list", numReplicates)
				rtList[[peakGroup]][[sampleGroup]][[itop]] <- vector("list", numReplicates)
			}

			for (isotop in 1:numIsotopes){
				for (replicate in 1:numReplicates){
					sample <- replicatesList[[replicate]]
					peakID <- isotopeList[[peakGroup]][[isotop]][[sample]] # peakid can be either the id of the peak or the gapfilled value, I know the name couses misinter..
					
					if (!is.null(peakID)){
						if (typeof(peakID) == "list" & peakID[[1]] == "gapfilled"){
							vals <- peakID[[2]]
							mzList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- vals[1,]
							intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- vals[2,]
							rtList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- vals[3,]
						} else{
							if (length(peakID)>1){ #If more than one peak exists with similar intensity and RT
								peakID <- getBestPeakID(chromDataList, peakID)
							}
							chrom <- chromDataList[[peakID]]
							mzList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- chrom[1,]
							intList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- chrom[2,]
							rtList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] <- chrom[3,]
						}
					}
				}
			}
		}
	}
	list(mzList, intList, rtList)
}
