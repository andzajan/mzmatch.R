PeakML.QCIntensityCorrection <- function (filename, ionisation="detect", Rawpath=NULL, outputfile, QCSample=NULL, RSD=0.3, writeRejected=FALSE)
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	numberOfpeakSets <- ncol(PeakTable[[1]])

	## Check which samples are QC's
	qcsamples <- which (PeakMLdata$phenoData==QCSample)
	if (length(qcsamples)<1)
	{
		stop ("No samples matching a given names for \"QC's\" are found. Please check if function parameter `QCSample` is set correctly.")
	} 

	Intensities.QC <- PeakTable[[1]][qcsamples,]
	Intensities.QC <- rbind(Intensities.QC, NULL)
	if (nrow(Intensities.QC)<2)
	{
		stop ("Only one QC sample is detected. QC RSD correction is not possible.")
	}
	
	na.vals <- which(is.na(Intensities.QC))

	QC.feature.meadians <- 

	if (length(na.vals)>0)
	{
		cat ("Warning, QC sample data contain missing values. Missing data are replaced with median value of current peak set.")
	}


	



	Intensities.all <- PeakTable[[1]]


	QC.RSD <- apply(Intensities.QC,2,var)
	QC.RSD <- sqrt(QC.RSD)/apply(Intensities.QC,2,mean)
	
	to.remove <- which(is.na(QC.RSD))
	to.remove <- append (to.remove, which(QC.RSD>RSD))
	to.remove <- unique(to.remove)

	## Write no matching
	if (length(to.remove)>0 & writeRejected==TRUE)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("QC_removed_",outputfile,sep=""), sets=c(1:numberOfpeakSets)[to.remove])
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
