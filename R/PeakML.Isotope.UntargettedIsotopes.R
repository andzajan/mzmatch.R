PeakML.Isotope.UntargettedIsotopes <- function(
        baseDir,
	outFileName,
	mzXMLSrc=NULL,
	outDirectory = "untargettedIsotops",
	peakMLFile="final_combined_related_identified.peakml",
	analyse = "databases",
	databases = c("kegg"),
	sampleGroups = NULL,
	layoutMtx = NULL,
	ppm = 3,
	trendPlots = NULL,
	fillGaps = "ALLPEAKS",
	useArea = FALSE,
	stdRTWindow = NULL,
	filterStringency=30,
	baseCorrection=FALSE,
	numSlaves = 1,
	label=1,
        exclude_from_plots=NULL
    ){
	# PRE:
	#	peakMLFiles: the complete peakml dataset
	#	molFormulaFile: file containing the list of molecules whoes isotops has to be found out
	#	outDirectory: is the outDirectory where the output has to be saved
	#	mzXMLSrc: is the source of the mzXML file is not in the current working outDirectory
	#	include_ionisation: not necessary in this cases
	#	ionisation: set this if include_ionisation=TRUE
	#	loadSavedData: load from the saved peakml file
	# 	sampleType: the sample type eg. NEG, POS etc

	# POST:
	#	vector containing the list of isotops
	## Reads the peakml file & prepare the parameters to scan for isotops
	## --------------------------------------------------------------------

	labels <- c("C13", "N15")
	if (label <= length(labels)){
		label <- labels[label]
	} else {
		stop ("Please specify the correct isotope used for labelling")
	}
        print(label)

	cat("Indentifying isotopes in sample\n")
	setwd (baseDir)

	if (is.null(mzXMLSrc)){
		stop ("Please provide the location of the raw data (mzXML) files ")
	}

	if (file.exists("cpData.Rdata") == TRUE){
		load("cpData.Rdata")
	} else{
		chromPeakData <- PeakML.Read(peakMLFile, ionisation = "neutral", mzXMLSrc)
		save("chromPeakData", file="cpData.Rdata")
	}

	peakDataMtx <- chromPeakData$peakDataMtx
	chromDataList <- chromPeakData$chromDataList
	sampleClasses <- chromPeakData$sampleClasses
	sampleNames <- chromPeakData$sampleNames
	massCorrection <- PeakML.Methods.getMassCorrection(filename=peakMLFile)
	phenoData <- PeakML.Methods.getPhenoData(sampleClasses, sampleNames, peakDataMtx)
    	## sampleType <- sampleType

	if (is.null(sampleGroups)) sampleGroups <- unique(phenoData)		# To enable the user to change the order of the samples

	if (is.null(trendPlots)) trendPlots <- c("RATIO","TREND", "LABELLED", "TOTRATIO")

	if (length(sampleGroups)>22) {
		if (is.null(layoutMtx)){
			stop("You have more than 22 samples to plot. Please specify an appropriate layout matrix.\n")
		} else {
			plotOrder <- c(sampleGroups, trendPlots)
		}
	} else {
		numSG <- length(sampleGroups)
		if (numSG < 7){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,2, 3,4,5,6,7,8,9, 10,11,11,12,12,13,13),3,7, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 7-numSG), trendPlots)
		} else if (numSG >= 7 & numSG <=14){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,2, 3,4,5,6,7,8,9, 10,11,12,13,14,15,16, 17,18,18,19,19,20,20),4,7, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 14-numSG), trendPlots)
		} else if (numSG > 14 & numSG <=21){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,2, 3,4,5,6,7,8,9, 10,11,12,13,14,15,16, 17,18,19,20,21,22,23, 24,25,25,26,26,27,27),5,7, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 21-numSG), trendPlots)
		}
	}

        exclude <- NULL
        if (!is.null(exclude_from_plots)) {
            exclude <- which(sampleGroups %in% exclude_from_plots)
            if (length(exclude_from_plots)!=length(exclude)){
                stop("The sample group you wanted to exlude from the final output does not exist in the sample groups")
            }
        }

	PeakML.Isotope.processUntargettedIsotopes(peakMLFile, analyse, databases, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, numSlaves, label, exclude_from_plots)

}

