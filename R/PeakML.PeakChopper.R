PeakML.PeakChopper <- function (filename, ionisation="detect", Rawpath=NULL, outputfile, peaktailI=95, ignoreGroups=NULL, filterOnPeakMax=FALSE)
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	numberOfpeakSets <- length(unique(PeakMLdata$peakDataMtx[,10]))
	PeakMLdataold <- PeakMLdata
	
	if (!is.null(ignoreGroups))
	{
		samplesToIgnore <- NULL
		for (ig in ignoreGroups)
			{
				samplesToIgnore <- append(samplesToIgnore, which(PeakMLdata$phenoData==ig))
			} 	
	}
	
	# single peaks to remove, as those does not match 3 scan length criteria after filtering
	toremove <- NULL
	
	# counter for full peak sets which does not match criteria after filtering
	toremoveset <- NULL
	
	system.time({
	for (nset in 1:numberOfpeakSets)
	{
		counter <- 0
		hits <- which(PeakMLdata$peakDataMtx[,10]==nset)
		
		if (!is.null(ignoreGroups))
		{
			ignorepeakset <- which(PeakMLdata$peakDataMtx[hits,9]%in%samplesToIgnore)
		}
		
		if (length(hits)>1)
		{
			minRT <- PeakMLdata$peakDataMtx[hits,5]
			maxRT <- PeakMLdata$peakDataMtx[hits,6]
			minRTm <- median(minRT)
			maxRTm <- median(maxRT)
			
			# maxint before chopping
			intmax <- rep(NA,length(hits))
			# maxint after chopping
			intmaxc <- rep(NA,length(hits))
			# int on the left after chopping
			intl <- rep(NA,length(hits))
			# int on the right tail after chopping
			intr <- rep(NA,length(hits))
			
			
			for (chrnum in 1:length(hits))
			{
				chromatogram <- PeakMLdata$chromDataList[[hits[chrnum]]]
				intmax[chrnum] <- max(chromatogram[2,])
								
				tokeep <- which(chromatogram[3,]>=minRTm & chromatogram[3,] <= maxRTm)
				
				if (length(tokeep)<=2)
				{
					toremove <- append(toremove, hits[chrnum])
					counter <- counter+1
					intmaxc[chrnum] <- intmax[chrnum] 
					
				} else
				{
					intmaxc[chrnum] <- max(chromatogram[2,tokeep])
					intl[chrnum] <- chromatogram[2,tokeep][1]
					intr[chrnum] <- chromatogram[2,tokeep][length(tokeep)]
					PeakMLdata$chromDataList[[hits[chrnum]]] <- chromatogram[,tokeep]
				}	
			}
			
			if (counter==length(hits))
			{
				toremoveset <- append(toremoveset,nset)
			}
		
			# if maxint != maxint after chopping, discard peak set. More than one maximum in the peak set.
			
			if (filterOnPeakMax==TRUE & length(ignorepeakset)==0)
			{
				intmaxc[is.na(intmaxc)] <- 1
				intmax[is.na(intmax)] <- 1
				if (!all(intmaxc==intmax))
				{
					toremoveset <- append(toremoveset, nset)
					toremove <- append(toremove, hits)
				}
			}
			
			# peaktailI exceeded discard peakset
			
			Inten <- median(intmaxc,na.rm=TRUE)
			Intlm <- median(intl, na.rm=TRUE)
			Intrm <- median(intr, na.rm=TRUE)
			if (is.na(Intlm)) Intlm <-1
			if (is.na(Intrm)) Intrm <-1
			
			if ((Intlm/Inten)*100 > peaktailI | (Intrm/Inten)*100 > peaktailI)
			{
				toremoveset <- append(toremoveset, nset)
				toremove <- append(toremove, hits)	
			}
			
			
		}
		
	}
	})
	
	toremove <- unique (toremove)
	toremoveset <- unique (toremoveset)
	
	## Write discarded sets 
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdataold, outputfile=paste("removed_",outputfile,sep=""), sets=toremoveset)
	rm (PeakMLdataold)
	
	if (!is.null(toremove))
	{
		keeppeaks <- c(1:nrow(PeakMLdata$peakDataMtx))[-c(toremove)]
		PeakMLdata$peakDataMtx <- PeakMLdata$peakDataMtx[keeppeaks,]
		PeakMLdata$chromDataList <- PeakMLdata$chromDataL[keeppeaks]
		PeakMLdata$sampleGroups <- PeakMLdata$sampleGroups[keeppeaks]
	}
	
	if (!is.null(toremoveset) & !is.null(PeakMLdata$GroupAnnotations) )
	{	
		for (ann in 1:length(PeakMLdata$GroupAnnotations))
		{
			PeakMLdata$GroupAnnotations[[ann]] <- PeakMLdata$GroupAnnotations[[ann]][-c(toremoveset)]
		}
	}
		
	PeakML.Write (peakMLdata=PeakMLdata, outFileName=outputfile)
	
}
