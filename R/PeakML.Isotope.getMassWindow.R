PeakML.Isotope.getMassWindow <- function(mass, nIsotope, ppm, element){
	# PRE: 
	# 	mass : unlabelled mass of the compound
	#	nIsotope: number of labelled carbons
	#	ppm: the ppm window required.
	#	element: the element whose isotopes we are interest in
	# POST: 
	# 	Mass window of isotop
	# N15 <- 15.0001088982 - N14 <-14.0030740048
	
	fMass <- 0
	if (element=="C") fMass <- 1.0033
	if (element=="N") fMass <- 0.9970
	
	IsoMass <- fMass * nIsotope
	rv <- PeakML.Methods.getPPMWindow(mass + IsoMass, ppm)
	rv
}

