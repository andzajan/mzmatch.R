PeakML.Methods.getElements <- function(formula, element){
	# Returns the number of carbons in a compounds from the formula
	
	numElement <- 0
	elements <- c("Ca", "Cl", "Co", "Cr", "Cu", "Ce", "Cm", "Cf", "Na", "Ni", "Ne", "Nd", "Np" , "No")
	bonds <- c("CH", "CC", "CO", "CS", "CF", "CB", "CN","CP", "CR", "NO", "NH", "NC")
	
	formula <- as.character(formula)
	index <- grep(element, strsplit(formula,"")[[1]])
	if (!length(index)==0){
		subStr <- substr(formula, index, index + 2)
		if (!substr(subStr,1,2) %in% elements){
			numElement <- as.numeric(gsub("\\D", "", subStr))
		}
		if (substr(subStr,1,2) %in% bonds){
			numElement <- 1
		}
	}
	if (is.na(numElement)){
		if (element == substr(formula, nchar(formula), nchar(formula))){
			numElement <- 1
		} else {
			numElement <- 0
		}
	}
	numElement
}
