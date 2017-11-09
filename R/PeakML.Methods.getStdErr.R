PeakML.Methods.getStdErr <- function(x){
	sd(x, na.rm=TRUE)/sqrt(length(na.omit(x)))
	}
