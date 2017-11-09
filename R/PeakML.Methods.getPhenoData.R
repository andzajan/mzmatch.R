PeakML.Methods.getPhenoData <- function(sampleGroups, sampleNames, peakDataMtx){
	# This function will return the sample groups (eg: BLK, BLK, BLK, CO, CO, CO etc)	
	# PRE: 
	#	sampleClasses: list of available sample classes
	#	sampleNames: sample names
	#	peakDataMtx: the peak data matrix
	
	getSampleClass <- function (x)
	{
		phenoData <- sampleGroups[sGroups[peakDataMtx[,9]==x][1]]
		phenoData
	}

	sGroups <- peakDataMtx[,11]
	phenoData <- unlist(lapply(1:length(sampleNames),getSampleClass)) 
	phenoData
}
