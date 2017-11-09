PeakML.Isotope.processTargettedIsotopes <- function (molFormulaFile, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, baseCorrection, label, exclude_from_plots){

	readTargetsFromFile<- function(inputFile){
		# PRE:
		#	inputFile: a tab separated csv file that conforms to RCreateXMLDB format e.i. "id", "name", "formula" as column headings
		# POST:
		#	Contents of the input file as a dataFrame that has masses added in column "mass"
		# Load the java project where the java class is located with dummy parameters

		molFrame <- read.table(inputFile, sep="\t", header=TRUE) # read the file as a data frame
		molMasses <- NULL
		for (imol in 1:length(molFrame$formula)){
                    mass <- try(PeakML.Methods.formula2mass(as.character(molFrame$formula)[imol]), silent=TRUE)
                    if(is.numeric(mass)){
                        molMasses <- c(molMasses, mass)
                    }else{
                        molMasses <- c(molMasses, NA)
                    }
		}
		if(is.null(molFrame$mass)) molFrame$mass <- molMasses
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
	molFrame <- readTargetsFromFile(molFormulaFile) # reading the molformula file
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
                cat(metName, ":\n")

		metFormula <- as.character(molFrame$formula[i])
		numElements <- PeakML.Methods.getElements(metFormula, element)
		metMass <- as.numeric(molFrame$mass[i])
                if(is.na(metMass)){
                    cat("\tThe metabolite ", metName, " with formula (", metFormula,"), does not exist, skipping. \n")
                    next()
                }
                metComment <- as.character(molFrame$comment[i])

                readRT <- as.character(molFrame$rt[i])
                if(!is.na(readRT)){
			if(length(strsplit(readRT, "-")[[1]])>1){
                	    lb <- as.numeric(strsplit(readRT, "-")[[1]][1])
                	    ub <- as.numeric(strsplit(readRT, "-")[[1]][2])
                    	    readRT <- PeakML.Methods.getRTWindowFromString(lb, ub)

                    	    stdRT <- readRT[[1]] * 60
                            stdRTWindow <- readRT[[2]] * 60
			}else{
			    stdRT <- as.numeric(readRT) * 60
			}
                } else {
                    stdRT <- NULL
                }


		if (is.null(molFrame$follow[i])){
			followCarbon <-  numElements + 1
		}else{
			followCarbon <- as.numeric(molFrame$follow[i])+1
		}

		if (is.na(molFrame$follow[i])) followCarbon <- 1

		if ('include' %in% colnames(molFrame)){
			if(as.character(molFrame$include[i]) == "") next()
		}

		if (numElements==0){
			cat("\tThe metabolite ", metName, " does not contain the preferred element (", element,"), hence skipping. \n")
			next()
		}
		numElements <- numElements + 1 				# This is to account for the basal peaks as well.

		cat ("\tIdentifying isotopes: ")
		# get the UID of isotops
		isotopeList <- PeakML.Isotope.getIsotopes (peakDataMtx, mzXMLSrc, sampleNames, label, numElements, metMass, ppm, massCorrection, baseCorrection, stdRT, stdRTWindow, fillGaps)

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
			cat("\tNo peaks found with mass:", metFormula, " (" , metMass ,") and its isotopes\n")
		}
	}
	save("molAbunList", file="abunList.Rdata")
	dev.off()
}
