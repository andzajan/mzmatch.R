PeakML.QCRSDFilter2 <- function (filename, ionisation="detect", Rawpath=NULL, outputfile, QCSample=NULL, RSD=0.3, writeRejected=FALSE, sampleFraction=0.7)
{
  
  filename="gapfilled_chop_dclean.peakml"
  ionisation="detect"
  Rawpath=NULL
  outputfile="gapfilled_chop_dclean_QCRSD.peakml"
  QCSample="QC"
  RSD=100
  writeRejected=TRUE
  sampleFraction=1
	
  
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
	Intensities.QC <- rbind(Intensities.QC, NULL)  #
	if (nrow(Intensities.QC)<2)
	{
		stop ("Only one QC sample is detected. QC RSD correction is not possible.")
	}
	
	na.count <- apply(Intensities.QC,2,function (x) length(which(is.na(x))))   
	na.count <- 1-(na.count/nrow(Intensities.QC)) 
	

	QC.RSD <- apply(Intensities.QC,2,var,na.rm=TRUE)
	QC.RSD <- sqrt(QC.RSD)/apply(Intensities.QC,2,mean,na.rm=TRUE)
	
	to.remove <- which(is.na(QC.RSD))
	to.remove <- append(to.remove, which(QC.RSD>=RSD ))
	to.remove <- append(to.remove, which(na.count < sampleFraction))   
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
