PeakML.FilterPeakGroups <- function (filename, ionisation="detect", Rawpath=NULL, ppm=5, rtwin=60, outputfile)
{
	version.1 <- get("version.1",envir=.GlobalEnv)

	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	Masses <- apply(PeakTable[[2]],2,median,na.rm=TRUE)
	RTs <- apply(PeakTable[[3]],2,median,na.rm=TRUE)

	MassOrder <- order(Masses)
	Masses <- sort(Masses)
	RTs <- RTs[MassOrder]

	validSetIndex <- rep(1,length(Masses))

	for (mm in 1:length(Masses))
	{
		if(!is.na(Masses[mm]))
		{
			massWindow <- PeakML.Methods.getPPMWindow(Masses[mm],ppm)
			HITS <- which(Masses >= massWindow [[1]] & Masses<=massWindow[[2]])
			if (length(HITS)>1)
			{
				selRTs <- RTs[HITS]
				RTmax <- selRTs[selRTs==min(selRTs)[1]][1]+rtwin
				HITS <- HITS[which (selRTs<=RTmax)]
				if (length(HITS)>1)
				{
					setsToCheck <- MassOrder[HITS]
					nsamples <- rep(NA,length(setsToCheck))
					maxint <- rep(NA,length(setsToCheck))
					rtwith <- rep(NA,length(setsToCheck))
					for (setn in 1:length(setsToCheck))
					{
						#cat (setn,"\n")
						## number of samples in peak set
						nsamples[setn] <- nrow(rbind(PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[,10]==setsToCheck[setn],],NULL))
						## max signal intensity
						maxint[setn] <- max(PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[,10]==setsToCheck[setn],8])
						## delta RT between largest and smalles RT in peak set
						rtwith[setn] <- max(PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[,10]==setsToCheck[setn],6])-min(PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[,10]==setsToCheck[setn],5])
					}
					setsToKeep <- which(nsamples==max(nsamples))
					if (length(setsToKeep)>1)
					{
						setsToKeep <- which(nsamples==max(nsamples) & maxint==max(maxint))
					}
					if (length(setsToKeep)>1)
					{
						setsToKeep <- which(nsamples==max(nsamples) & rtwith==max(rtwith) & rtwith==max(rtwith))
					}
					
					validSetIndex[setsToCheck[-c(setsToKeep)]] <- 0
				}
			}
			Masses[HITS] <- NA
		}
	}

	## Write valid sets 
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=which(validSetIndex==1))

	## Write discarded sets 
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("removed_",outputfile,sep=""), sets=which(validSetIndex==0))
}



