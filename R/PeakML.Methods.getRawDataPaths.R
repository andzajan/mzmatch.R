PeakML.Methods.getRawDataPaths <- function(PeakMLtree, Rawpath = NULL){
	# PRE: 
	#	PeakMLtree, path of raw files if different from the ones referened in the peakml files
	# POST:
	#	Return the samples names and the full path to where the .mzXML files are located
	sampleNames <- sapply(getNodeSet(PeakMLtree,"/peakml/header/measurements/measurement/label"),xmlValue)
	folders <- sapply(getNodeSet(PeakMLtree,"/peakml/header/measurements/measurement/files/file/location"),xmlValue)
	filenames <- sapply(getNodeSet(PeakMLtree,"/peakml/header/measurements/measurement/files/file/name"),xmlValue)
	peakmlfiles <- grep(".peakml",filenames)
	if (length(peakmlfiles)!=0)
	{
		folders <- folders[-c(peakmlfiles)]
		filenames <- filenames[-c(peakmlfiles)]
	}

	if (!is.null(Rawpath)){
		dirContent <- dir(Rawpath, recursive = TRUE, full.names=TRUE)
		fileid <- rep(NA,length(sampleNames))
		for (filenum in 1:length(sampleNames))
		{
			hit <- grep (filenames[filenum],dirContent, fixed=TRUE)
			if (length(hit)!=0)
			{
				fileid[filenum] <- hit[1]
			}
		}
	  	if(length(which(is.na(fileid)))==0){
			rawDataPaths <- dirContent[fileid]
			cat(paste("Raw data file located at: ", rawDataPaths, "\n", sep=""))
		} else{
			cat(paste("Raw data file: ",sampleNames[which(is.na(fileid))],"cant't be found in folder: ", Rawpath, "\n", sep=""))
			rawDataPaths <- NULL
		}
	} else{
		rawDataPaths <- paste(folders,filenames,sep="/")
		mzFiles <- which(file.exists(rawDataPaths)==TRUE)
		if (length(mzFiles) == length(rawDataPaths)){
			cat(paste("Raw data file located at: ", rawDataPaths, "\n", sep=""))
		} else {
			cat(paste("Raw data file can't be read from location': ", rawDataPaths, "\n", sep=""))
			rawDataPaths <- NULL
		}
	}
	list(sampleNames,rawDataPaths)
}
