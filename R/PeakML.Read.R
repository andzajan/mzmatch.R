PeakML.Read <- function(filename, ionisation = "detect", Rawpath=NULL){
	# Return the peak and chrom data from the file specified
	# PRE: peakML filename
	# POST:
	#	list(peakDataMtx, chromDataList)
	version.1 <- get("version.1",envir=.GlobalEnv)

	sampleLookUP <- function (x)
	{
		phenoData <- sampleClasses[sampleGroups[peakDataMtx[,9]==x][1]]
		phenoData
	}
	
	getRTScans <- function(PeakMLtree, sampleNames){
		corRT <- vector ("list",length(sampleNames))
		rawRT <- vector ("list",length(sampleNames))
		for (measID in 1:length(sampleNames)){
			# Create a string containing complete XPATH (subsetiing indexes does not work directly from R input) command, and then avluate it in R.
			STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/header/measurements/measurement[id=",measID-1,"]/scans/scan/retentiontime\"),xmlValue)",sep="")
			corRT[[measID]] <- eval (parse(text=STR))
			STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/header/measurements/measurement[id=",measID-1,"]/scans/scan/annotations/annotation/label/text ()\"),xmlValue)",sep="")
			## At this moment only RT_raw is used, but can be extended to support al annotations from header if there are any extra annotations.
			labels <- eval (parse(text=STR))
			STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/header/measurements/measurement[id=",measID-1,"]/scans/scan/annotations/annotation/value/text ()\"),xmlValue)",sep="")
			values <- eval(parse(text=STR))
			rawRT[[measID]] <- values[labels=="RT_raw"]
		}
		list(corRT, rawRT)
	}

	st <- system.time(PeakMLtree <- xmlInternalTreeParse(filename))
	
	## If olda mzmatch library is used mass will be corrected for ionisation mode and "neutral" mass stored in peakml file.
	if (version.1==TRUE)
	{
		massCorrection <- PeakML.Methods.getMassCorrection(PeakMLtree, ionisation)
	} else
	{
		massCorrection <- 0
	}

	samPath <- PeakML.Methods.getRawDataPaths(PeakMLtree, Rawpath)
	sampleNames <- samPath[[1]]
	rawDataFullPaths <- samPath[[2]]

	system.time(chromDataList <- PeakML.Methods.getChromData(filename, PeakMLtree, massCorrection))

	peakDataMtx <- PeakML.Methods.getPeakData(PeakMLtree, chromDataList)
	
	rtScanList <- getRTScans(PeakMLtree, sampleNames)

	sampleGroups <- peakDataMtx[,11]
	sampleClasses <- sapply (getNodeSet(PeakMLtree,"/peakml/header/sets/set/id"),xmlValue)
	SetMeasurementids <- sapply(getNodeSet(PeakMLtree,"/peakml/header/sets/set/measurementids"),xmlValue)
	if (length(SetMeasurementids)!=0)
	{
		SetMeasurementids <- lapply (1:length(SetMeasurementids),function(i) {base64decode(SetMeasurementids[i], what="integer",endian="swap")+1})
		phenoData <- rep (NA,length(sampleNames))
		for (cl in 1:length(SetMeasurementids))
		{
			phenoData[SetMeasurementids[[cl]]] <- sampleClasses[cl]
		}
	} else
	{
		phenoData <- NA
	}

	rv = list()
	rv$peakDataMtx <- peakDataMtx
	rv$chromDataList <- chromDataList
	rv$sampleClasses <- sampleClasses
	rv$sampleNames <- sampleNames
	rv$massCorrection <- list(massCorrection,PeakML.Methods.getProtonCoef(PeakMLtree, ionisation)[[2]])
	rv$sampleGroups <- sampleGroups
	rv$phenoData <- phenoData
	rv$correctedRTList <- rtScanList[[1]]
	rv$rawRTList <- rtScanList[[2]]
	rv$fileName <- filename
	rv$rawDataFullPaths <- rawDataFullPaths
	if (length(SetMeasurementids)!=0)
	{
		rv$GroupAnnotations <- PeakML.Methods.getGroupAnnotations(PeakMLtree)
	} else
	{
		rv$GroupAnnotations <- NA
	}
	rv
}
