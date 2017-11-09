PeakML.Normalization <- function(filename,ionisation="detect",Rawpath=NULL,outputfile,values=NULL)
{
	PeakMLdata <- PeakML.Read(filename=filename,ionisation=ionisation,Rawpath=Rawpath)

	## Append peakdata with a vector of normalisation values which should be used for every vector of intensities

	## if values=NULL use TIC's of the raw data files for normalisation
	if (is.null(values))
	{
		values <- rep(NA,length(PeakMLdata$sampleNames))
		for (file in 1:length(PeakMLdata$sampleNames))
		{
			rawfile <- xcmsRaw(PeakMLdata$rawDataFullPaths[file])
			values[file] <- sum(rawfile@tic)
		}
	}
	
	normvals <- rep(NA,nrow(PeakMLdata$peakDataMtx))
	for (nval in 1:length(PeakMLdata$sampleNames))
	{
		normvals[PeakMLdata$peakDataMtx[,9]==nval] <- values[nval]
	}


	## Normalize data
	for (chrnum in 1:length(normvals))
	{
		PeakMLdata$chromDataList[[chrnum]][2,] <- PeakMLdata$chromDataList[[chrnum]][2,]/normvals[chrnum]*10^4
	}

	
	## Write all of this out
	PeakML.Write(peakMLdata=PeakMLdata, outFileName=outputfile)
}
