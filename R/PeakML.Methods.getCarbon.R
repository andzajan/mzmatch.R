PeakML.Methods.getCarbon <- function(formula){
	# Returns the number of carbons in a compounds from the formula
	
	carbons <- 0
	elements <- c("Ca", "Cl", "Co", "Cr", "Cu", "Ce", "Cm", "Cf")
	bonds <- c("CH", "CC", "CO", "CS", "CF", "CB", "CN","CP", "CR")
	
	formula <- as.character(formula)
	index <- grep("C", strsplit(formula,"")[[1]])
	if (!length(index)==0){
		subStr <- substr(formula, index, index + 2)
		if (!substr(subStr,1,2) %in% elements){
			carbons <- as.numeric(gsub("\\D", "", subStr))
		}
		if (substr(subStr,1,2) %in% bonds){
			carbons <- 1
		}
	}
	carbons
}



