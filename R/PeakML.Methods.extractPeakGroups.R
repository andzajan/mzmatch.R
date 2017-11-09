PeakML.Methods.extractPeakGroups <- function(PeakMLData, outputfile, sets)
{
	if (is.null(PeakMLData$rawDataFullPaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}
	if (max(sets)>max(PeakMLData$peakDataMtx[,10]))
	{
		cat ("Peak group index in \"sets\" parameter is larger than number of the peak groups in peakml file. \n")
		stop ()
	}

	whattoextract <- which(PeakMLData$peakDataMtx[,10]%in%sets)

	PeakMLData$peakDataMtx <- PeakMLData$peakDataMtx[whattoextract,]
	PeakMLData$chromDataList <- PeakMLData$chromDataList[whattoextract]
	PeakMLData$sampleGroups <- PeakMLData$sampleGroups[whattoextract]

	if (!is.null(PeakMLData$GroupAnnotations))
	{
		for (ann in 1:length(PeakMLData$GroupAnnotations))
		{
			PeakMLData$GroupAnnotations[[ann]] <- PeakMLData$GroupAnnotations[[ann]][sets]
		}
	}

	PeakML.Write (peakMLdata=PeakMLData, outFileName=outputfile)	
}
