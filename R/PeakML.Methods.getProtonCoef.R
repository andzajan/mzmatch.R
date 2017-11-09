PeakML.Methods.getProtonCoef <- function(PeakMLtree, ionisation="detect"){
	# PRE:
	# 	Java project, inonisation (NULL, detect, negative, positive)
	# POST:
	#	protonCoef, mass correction factor
	
	#Detect ionisation mode, only first scan from first sample is used, no ionisation switching is supported yet
	if (ionisation=="detect") {
		ionisation <- getNodeSet(PeakMLtree,"/peakml/header/measurements/measurement/scans/scan/polarity")
		ionisation <- tolower(unique(sapply(ionisation,xmlValue)))
		if (length(ionisation)>1)
		{
			cat ("PeakML file contains data from more than one ionisation mode. At this moment such data can't be hanled with mzmatch.R","\n")
			stop ()
		}
	}

	if (ionisation=="positive") {
		protonCoef <- 1
	} else if (ionisation=="negative") {
		protonCoef <- -1
	} else	protonCoef <- 0
	
	list(protonCoef,ionisation)
}
