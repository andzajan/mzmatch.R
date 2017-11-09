PeakML.Methods.getMassCorrection <- function(PeakMLtree=NULL, ionisation="detect",filename=NULL){
	# PRE:
	# 	the jave project, ionisation mode (see PeakML.Methods.getProtonCoef)
	# POST:
	#	returns mass correction as numeric
	if (!is.null(filename))
	{
		st <- system.time(PeakMLtree <- xmlInternalTreeParse(filename))
	}
	protonMass <- PeakML.Methods.getProtonMass()
	protonCoef <- PeakML.Methods.getProtonCoef(PeakMLtree, ionisation)
	rv <- protonMass * protonCoef[[1]]
	rv
}
