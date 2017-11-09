PeakML.Methods.formula2mass <- function(formula){
	java_project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	rv <- .jcall (java_project, returnSig="D", method="formulaToMass", formula)
	rv
}
