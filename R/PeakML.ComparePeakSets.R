PeakML.ComparePeakSets <- function(standard_filename, filename, stdionisation="detect",ionisation="detect", stdRawpath=NULL, Rawpath=NULL, outputfile,ppm=5,rtwin=20,checkIntensity=TRUE,PeakShapeCor=TRUE,PeakShapeCor.thr=0.7)
{
	st <- system.time (stdPeakMLdata <- PeakML.Read (standard_filename,stdionisation,Rawpath=stdRawpath))
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))

	rawdatafullpaths <- PeakMLdata$rawDataFullPaths
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	rawdatafullpaths <- stdPeakMLdata$rawDataFullPaths
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	stdPeakMLDataTable <- PeakML.Methods.getCompleteTable (stdPeakMLdata,sumintensity=FALSE)
	stdpeakMasses <- apply(stdPeakMLDataTable[[2]],2,median,na.rm=TRUE)-stdPeakMLdata$massCorrection[[1]]
	stdpeakRTs <- apply(stdPeakMLDataTable[[3]],2,median,na.rm=TRUE)

	PeakMLDataTable <- PeakML.Methods.getCompleteTable (PeakMLdata,sumintensity=FALSE)
	peakMasses <- apply(PeakMLDataTable[[2]],2,median,na.rm=TRUE)-PeakMLdata$massCorrection[[1]]
	peakRTs <- apply(PeakMLDataTable[[3]],2,median,na.rm=TRUE)
	peakIntensities <- apply(PeakMLDataTable[[1]],2,max,na.rm=TRUE)

	selectedSets <- NULL
	stdMatchedSets <- NULL
	for (setnum in 1:length(stdpeakMasses))
	{
		setMass <- unlist(PeakML.Methods.getPPMWindow(stdpeakMasses[setnum],ppm))
		setRT <- c(stdpeakRTs[setnum]-rtwin,stdpeakRTs[setnum]+rtwin)

		hits <- which(peakMasses>=setMass[1] & peakMasses<=setMass[2] & peakRTs>=setRT[1] & peakRTs<=setRT[2])
		if (PeakShapeCor==TRUE & length(hits)>0)
		{
			peak.cors <- rep(NA,length(hits))
			maxint <- function (chrom,data)
			{
				max(data[[chrom]][2,])
			}
			stdChroms <- which(stdPeakMLdata$peakDataMtx[,10]==setnum)
			maxints <- sapply (stdChroms,maxint,data=stdPeakMLdata$chromDataList)
			stdchrom <- stdChroms[which(maxints==max(maxints))[1]]
			stdchrom <- stdPeakMLdata$chromDataList[[stdchrom]]
			for (hitn in 1:length(hits))
			{
				Chroms <- which(PeakMLdata$peakDataMtx[,10]==hits[hitn])
				maxints <- sapply (Chroms,maxint,data=PeakMLdata$chromDataList)
				chrom <- Chroms[which(maxints==max(maxints))[1]]
				chrom <- PeakMLdata$chromDataList[[chrom]]
				peak.cors[hitn] <- PeakML.Methods.PeakShapeCorrelation(stdchrom,chrom)[1]
			}
			hits <- hits[which(hits>PeakShapeCor.thr)]
			if (length(hits)>1)
			{
				hits <- hits[which(hits==max(hits))[1]]
			}
		}

		if (checkIntensity==TRUE & length(hits)>1)
		{
			maxint <- max(peakIntensities[hits])
			hits <- hits[which(peakIntensities[hits]==maxint)[1]]
		}
		if (length(hits)!=0)
		{
			selectedSets <- append(selectedSets,hits)
			stdMatchedSets <- append(stdMatchedSets,setnum)
		}
	}
	selectedSets <- unique(selectedSets)
	stdMatchedSets <- unique(stdMatchedSets)

	if (length(selectedSets)>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=selectedSets)
	}

	## Write no matching
	if (length(c(1:length(peakRTs))[-c(selectedSets)])>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("Not_matched_",outputfile,sep=""), sets=c(1:length(peakRTs))[-c(selectedSets)])
	}

	# Not macthed from std file
	if (length(c(1:length(stdpeakRTs))[-c(stdMatchedSets)])>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=stdPeakMLdata, outputfile=paste("STD_not_matched_",standard_filename,sep=""), sets=c(1:length(stdpeakRTs))[-c(stdMatchedSets)])
	}
}
