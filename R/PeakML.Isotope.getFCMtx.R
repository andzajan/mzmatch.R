PeakML.Isotope.getFCMtx <- function (trendList, sampleGroups, followCarbon){

	numCarbons <-  length(trendList[[1]])
	plotMtx <- matrix(nrow = 1, ncol = length(sampleGroups))
	dimnames(plotMtx) <- list(1, sampleGroups)
	
	for (sam in 1:length(sampleGroups)){
		VAL <- trendList[[sam]][[followCarbon]]
		if (!is.null(VAL)){
			plotMtx[1, sam] <- VAL
		} else {
			plotMtx[1, sam] <- 0
		}
	}
	
	plotMtx
}
