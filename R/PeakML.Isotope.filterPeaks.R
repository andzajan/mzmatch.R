PeakML.Isotope.filterPeaks <- function(peakDataMtx, rtWindow, massWindow){
	# This function will identify the isotops in the mass and rt window given from any matrix
	# PRE:
	#	peakData: this is a mass filterd group essntialy within similar mass window
	#	peakDataMtx: the origianl peak data matrix
	# 	caron: the number of carbons in the isotop
	#	ppm: carbon mass window for filtering out the isotops
	# POST:
	#	list of all the isotops satisfying the above criteria

	# Filtering the isotops out
	filterHits <- which(peakDataMtx[,1] >= massWindow[[1]]	& # isotop mass 
				  peakDataMtx[,1] <= massWindow[[2]] &
				  peakDataMtx[,4] >= rtWindow[[1]]	& # min RT
				  peakDataMtx[,4] <= rtWindow[[2]])	#& # max RT
				 # peakDataMtx[,9] == peakData[9])	 # sample ids should = ipeak appears after a while
	filterHits
}

