PeakML.Methods.getProtonMass <- function(){
	# POST:
	#	Return the proton mass
	
	proton <- 1.00782503214
	electron <- 0.00054857990924
	protonmass <- proton-electron
	protonmass
}
