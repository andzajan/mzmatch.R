PeakML.Isotope.filteroutUnlabelled <- function(isotopeList, numCarbons, fillGaps, sampleNames, stringency=30){
	
	possibleSamples <- numCarbons * length(sampleNames) 
	stringency <- stringency * possibleSamples / 100
	delGroups <- c()

	for (grp in 1:length(isotopeList)){
		
		itops <- unlist(isotopeList[[grp]][2:numCarbons])
	
		if (is.null(itops)){
			delGroups <- c(delGroups, grp)
		} else {
			if (fillGaps == "NONE"){
				if (length(itops) < stringency){
					delGroups <- c(delGroups, grp)
				}
			} else {
				if (length(which(itops == "gapfilled")) < stringency){
					delGroups <- c(delGroups, grp)
				}
			}
		}
	}
	if (length(delGroups)>0) isotopeList <- isotopeList[-delGroups]
	isotopeList
}

