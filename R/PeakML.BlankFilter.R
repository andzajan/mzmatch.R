PeakML.BlankFilter <- function (filename, ionisation="detect", Rawpath=NULL, outputfile, BlankSample=NULL, IgnoreIntensity=FALSE, detectedInNumOfBlanks=NULL, rtwindow=NULL)
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	numberOfpeakSets <- ncol(PeakTable[[1]])

	## Check which samples are blanks
	blanksamples <- which (PeakMLdata$phenoData==BlankSample)
	if (length(blanksamples)<1)
	{
		stop ("No samples matching a given names for \"blanks\" are found. Please check if function parameter `BlankSample` is set correctly.")
	} 
	Intensities.blank <- PeakTable[[1]][blanksamples,]
	Intensities.blank[is.na(Intensities.blank)] <- 0
	Intensities.blank <- rbind(Intensities.blank,NULL)
	if (nrow(Intensities.blank)>1)
	{
		Intensities.blank <- apply(Intensities.blank,2,max)
	}
	Intensities.rest <- PeakTable[[1]][-c(blanksamples),]
	Intensities.rest[is.na(Intensities.rest)] <- 0
	Intensities.rest <- apply(Intensities.rest,2,max,na.rm=TRUE)
	to.remove <- which(Intensities.blank >= Intensities.rest)

	if (IgnoreIntensity==TRUE)
	{
		to.remove <- which(Intensities.blank!=0)
	}

	if (!is.null(detectedInNumOfBlanks))
	{
		Intensities.blank <- PeakTable[[1]][blanksamples,]
		ndetect <- rep(NA, ncol(PeakTable[[1]]))
		for (nd in 1:length(ndetect))
		{
			ndetect[nd] <- nrow(Intensities.blank)-length(which(is.na(Intensities.blank[,nd])))
		}

		hits <- which (ndetect >= detectedInNumOfBlanks)
		
		to.remove <- to.remove[to.remove%in%hits]	
	}



	if (!is.null(rtwindow))
	{
		rtwindow <- rtwindow*60
		rt.times <- apply (PeakTable[[3]],2,median,na.rm=TRUE)
		hits <- which(rt.times >= rtwindow[1] & rt.times <= rtwindow[2])

		to.remove <- to.remove[to.remove%in%hits]
	}

	## Write no matching
	if (length(to.remove)>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("blank_removed_",outputfile,sep=""), sets=c(1:numberOfpeakSets)[to.remove])
	}
	# Not macthed from std file
	if (length(to.remove)>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=c(1:numberOfpeakSets)[-c(to.remove)])
	} else
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=c(1:numberOfpeakSets))
	}
}
