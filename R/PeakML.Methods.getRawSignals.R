PeakML.Methods.getRawSignals <- function(mzXMLSrc, sampleName, rtWindow, massWindow, massCorrection){
	# PRE:
	#	Return the raw signals for the given sample in the rt and mass window specified
	#	mzXMLSrc -> mzXML raw data source path
	#	sampleName -> Name of the sample
	#	rtWindow -> the retention time window as list [min, max]
	#	massWindow -> the mass window as list [min, max]
	#	massCorrection -> PeakML.Methods.getMassCorrection
	# POST:
	#	a table with masses,intensities,retentiontimes,scanids along the rows.
	
	curDir <- mzXMLSrc
	if (is.null(curDir)){
		curDir <- getwd()
	}
	curFile <- list.files(path=curDir, pattern=paste("^", sampleName, ".mzXML", sep=""), full.names=TRUE, recursive=TRUE)

	rawFile <- xcmsRaw(curFile)
	
	
	rtStart <- rtWindow[[1]]
	rtEnd <- rtWindow[[2]]
	if (rtEnd > max(rawFile@scantime)){
		rtEnd <- max(rawFile@scantime)
	}
	if (rtStart > max(rawFile@scantime)){
		rtStart <- max(rawFile@scantime)
	}
	if (rtEnd < min(rawFile@scantime)){
		rtEnd <- min(rawFile@scantime)
	}
	if (rtStart < min(rawFile@scantime)){
		rtStart <- min(rawFile@scantime)
	}

	massStart <- as.numeric(massWindow[[1]])+massCorrection
	massEnd <- as.numeric(massWindow[[2]])+massCorrection

	#  Extract data form raw data files
	if ((rtEnd-rtStart)<5){
		rData <- c(1,1,1)
	} else{
		rData <- rawMat (rawFile,mzrange=cbind(massStart, massEnd),rtrange=c(rtStart,rtEnd))
	}
	
	rData <- rbind (rData,NULL)
	## Removing values with repeating RT's
	## Row with largest intensity are selected
	rtRepeats <- as.numeric(names(which(table(rData[,1])>=2)))
	if (length(rtRepeats!=0)){
		for (z in 1:length(rtRepeats)){
			Csub <- which(round(rData[,1],5)==round(rtRepeats[z],5))
			Csub <- Csub[-c(which(rData[Csub,3]==max(rData[Csub,3]))[1])]
			rData <- rData[-c(Csub),]
			rData <- rbind(rData,NULL)
		}
	}

	scanids <- which(rawFile@scantime%in%rData[,1])
	## if RT correction was applied, scans should be extracted from raw RT's
	if (length(scanids)<3){
		scanids <- c(-1,-1,-1)
		retentiontimes <- c(-1,-1,-1)
		masses <- c(-1,-1,-1)
		intensities <- c(-1,-1,-1)
	} else {
		retentiontimes <- rData[,1]
		masses <- rData[,2]
		intensities <- rData[,3]
	}
	rm(rawFile, rData)
	rv <- rbind(masses,intensities,retentiontimes,scanids-1)
	rv
}
