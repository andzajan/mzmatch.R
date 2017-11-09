PeakML.Isotope.getTrendMtx <- function (trendList, sampleGroups){

	numCarbons <-  length(trendList[[1]])
	plotMtx <- matrix(nrow = numCarbons, ncol = length(sampleGroups))
	dimnames(plotMtx) <- list(c(1:numCarbons), sampleGroups)
	
	for (sam in 1:length(sampleGroups)){
		for (isotop in 1:numCarbons){
			VAL <- trendList[[sam]][[isotop]]
			if (!is.null(VAL)){
				plotMtx[isotop, sam] <- VAL
			} else {
				plotMtx[isotop, sam] <- 0
			}
		}
	}
	plotMtx
}

