PeakML.Methods.getPPMMass <- function(mass, ppm){
	# PRE:
	#	mass of the molecule
	# 	The required ppm as a number
	# POST:
	#	measured mass wrt given ppm
	
	rv <- mass * ppm * 1e-6
	rv
}
