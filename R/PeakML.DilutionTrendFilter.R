PeakML.DilutionTrendFilter <- function(filename,ionisation="detect",Rawpath=NULL,trendSets,p.value.thr=NULL,outputfile)
{

	version.1 <- get("version.1",envir=.GlobalEnv)

	#DilutionFactorfunc <- function (x) (1/2^x)

	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))

	ionisation <- PeakMLdata$massCorrection[[2]]
	massCorrection <- PeakMLdata$massCorrection[[1]]
	samplenames <- PeakMLdata$sampleNames
	rawdatafullpaths <- PeakMLdata$rawDataFullPaths
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	trendSetsindex <- NULL
	for (setin in 1:length(trendSets))
	{
		trendSetsindex <- append(trendSetsindex, which(PeakMLdata$sampleClasses==trendSets[setin]))
	}
		
	numOfDilPoints <- length(trendSetsindex)

	PeakSetCorrelations <- matrix(ncol=2,nrow=max(PeakMLdata$peakDataMtx[,10]),data=NA)
	## for each set caculate dilution trend.
	for (setnum in 1:max (PeakMLdata$peakDataMtx[,10]))
	{
		peaksetrows <- which(PeakMLdata$peakDataMtx[,10]==setnum)
		PeakSetData <- PeakMLdata$peakDataMtx[peaksetrows,]
		PeakSetData <- rbind(PeakSetData,NULL)
		trendintensities <- rep(NA,length(trendSetsindex))
		for (trset in 1:length(trendSetsindex))
		{
			hit <- which(PeakSetData[,11]==trendSetsindex[trset])
			if (length(hit)>0)
			{
				trendintensities[trset] <- median(PeakSetData[hit,8],na.rm=TRUE)
			}
		}

		nas <- which(is.na(trendintensities))
		nas <- append(nas,which(trendintensities==0))
		## Check if all measured values are consequtive. i.e there ar no gaps between dillution points
		if (length(nas)!=0 & length(nas)<length(trendintensities))
		{
			validvals <- c(1:length(trendintensities))[-c(nas)]
			minNas <- min(validvals)
			maxNas <- max(validvals)
			valid <- (maxNas-minNas+1)==length(validvals)
		} 
		if (length(nas)==0)
		{
			valid <- TRUE
		}

		if ((numOfDilPoints - length(nas))>=3 & valid==TRUE)
		{
			if (length(nas)==0)
			{
				RES <- cor.test(c(1:numOfDilPoints),log(trendintensities,base=2))
			} else
			{
				RES <- cor.test(c(1:numOfDilPoints)[-c(nas)],log(trendintensities[-c(nas)],base=2))
			}
			PeakSetCorrelations[setnum,] <- c(RES$estimate,RES$p.value)
		}
	}

	PeakMLdata$GroupAnnotations$dillution.corr <- PeakSetCorrelations[,1]
	PeakMLdata$GroupAnnotations$dillution.p.val <- PeakSetCorrelations[,2]
	if (is.null(p.value.thr))
	{
		PeakML.Write (peakMLdata=PeakMLdata,outFileName=outputfile)	
	} else
	{
		# Good peaksets
		HITS <- which(PeakSetCorrelations[,2]<=p.value.thr)
		if (length(HITS)>0)
		{
			PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=HITS)
		}
		# Bad peaksets
		HITS <- c(1:nrow(PeakSetCorrelations))[-c(HITS)]
		if (length(HITS)>0)
		{
			PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("discarded_",outputfile,sep=""), sets=HITS)
		}
	}
}
