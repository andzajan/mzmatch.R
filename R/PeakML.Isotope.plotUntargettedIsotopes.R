PeakML.Isotope.plotUntargettedIsotopes <- function (peakMLFile, molFrame, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, label, exclude_from_plots){

	addMasses<- function(molFrame){

		molMasses <- NULL
		for (imol in 1:length(molFrame$formula)){
			molMasses <- c (molMasses, PeakML.Methods.formula2mass(as.character(molFrame$formula)[imol]))
		}
		molFrame$mass <- molMasses
		molFrame
	}

	dir.create (outDirectory, showWarnings = FALSE)
	# To generate the tab delimited file
	csvFile <- paste(outDirectory, "/", outFileName, ".csv", sep ="")
	cat("\n", file=csvFile)
	# To generate the pdf plots
	pdfFile <- paste(outDirectory, "/", outFileName, ".pdf", sep="")
	pdf (file=pdfFile, paper="a4", height=10, width=7)
	# Create the layout for the pdf
	layout(layoutMtx, heights=c(0.4, rep(1, nrow(layoutMtx)-1)),TRUE)
	# Reading the list of targets in the mol formula file
	molFrame <- addMasses(molFrame)
	# To save the abundance matrix if needed for later processing .
	molAbunList <- vector("list", nrow(molFrame))

	if (massCorrection < 0){
		sampleType = "NEG"
	} else if (massCorrection > 0){
		sampleType = "POS"
	} else {
		sampleType = "NONE"
	}

	element <- substr(label,1,1)

	for (i in 1:nrow(molFrame)){
		metName <- as.character(molFrame$name[i])
		metFormula <- strsplit(as.character(molFrame$formula[i]), ",")[[1]][1]   #as.character(molFrame$formula[i])
		numElements <- PeakML.Methods.getElements(metFormula, element)
		metMass <- as.numeric(molFrame$mass[i])
		stdRT <- as.numeric(molFrame$rt[i]) * 60
		if(is.na(stdRT)) stdRT <- NULL


		if (is.null(molFrame$follow[i])){
			followCarbon <-  numElements + 1
		}else{
			followCarbon <- as.numeric(molFrame$follow[i])+1
		}

		if (!is.null(molFrame$follow[i]) && is.na(molFrame$follow[i])) followCarbon <- 1

		if ('include' %in% colnames(molFrame)){
			if(as.character(molFrame$include[i]) == "") next()
		}

		cat(metName, ":\n")
		if (numElements==0){
			cat("\tThe metabolite ", metName, " does not contain the preferred element (", element,"), hence skipping. \n")
			next()
		}
		numElements <- numElements + 1 					# This is to account for the basal peaks as well.

		cat ("\tIdentifying isotopes: ")
		# get the UID of isotops
		isotopeList <- PeakML.Isotope.getIsotopes (peakDataMtx, mzXMLSrc, sampleNames, label, numElements, metMass, ppm, massCorrection, baseCorrection, stdRT, stdRTWindow, fillGaps)

		if(!is.null(filterStringency)){
			# This is where to include the isotop level filter.
			isotopeList <- PeakML.Isotope.filteroutUnlabelled(isotopeList, numElements, fillGaps, sampleNames, stringency=30)
		}
                metComment <- "NULL"

		if (!is.null(unlist(isotopeList))){
			cat ("\n\tGenerating the plots. \n")
			isotopeChroms <- PeakML.Isotope.getChromData (isotopeList, chromDataList, phenoData, sampleGroups)
			PeakML.Isotope.plotSamples(isotopeChroms, metName, metFormula, metMass, metComment, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon, label, exclude_from_plots)

			ratioMtxList <- PeakML.Isotope.getRatioMtxList(isotopeChroms[[2]], sampleGroups, useArea, metName)

			molAbunList[[metName]] <- ratioMtxList

			cat("Metabolite: ", toupper(metName), "\t Formula: ", metFormula, "\tMass: ", metMass, "\n", file=csvFile, append=TRUE)
			cat("-----------------------------------------------------------------------------\n", file=csvFile, append=TRUE)
			for (pkgrp in 1:length(ratioMtxList)){
				cat("Group: ", pkgrp, "\n", file=csvFile, append = TRUE)
				sNames <- paste(sampleNames, collapse="\t")
				cat(paste(metName, metFormula, "\t", sNames , "\n"), file=csvFile, append=TRUE)
				write.table(ratioMtxList[[pkgrp]] , sep="\t", na= " ", file=csvFile, quote=FALSE, col.names=FALSE, append=TRUE)
				cat("\n", file=csvFile, append=TRUE)
			}
			cat("\n")

		} else{
			cat("\tNo peaks found with mass:", metMass ," and its isotopes\n")
		}
	}
	save("molAbunList", file="abunList.Rdata")
	dev.off()
}
