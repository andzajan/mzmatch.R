PeakML.Methods.getPPMWindow <- function(mass, ppm){
	# PRE:
	#	mass of the molecule
	# 	The required ppm as a number
	# POST:
	# 	list with the mass window wrt ppm

	ppmMass <- PeakML.Methods.getPPMMass(mass, ppm)
	
	rv <- vector("list",2)
	rv[[1]] <- mass-ppmMass
	rv[[2]] <- mass+ppmMass
	rv
}
