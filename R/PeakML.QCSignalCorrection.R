PeakML.QCSignalCorrection <- function(filename, outputfile, Rawpath=NULL, ionisation="detect", pbqc_groupid="PBQC", f=2/3, minIntThrSamples=NULL, scatter.plots=TRUE, sumintensity=FALSE, measurementOrder=NULL) 
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata, sumintensity=sumintensity)
	numberOfpeakSets <- ncol(PeakTable[[1]])
	
	inputdata <- PeakTable[[1]]
		
	#### order input data and extract the group	
	# extract order of samples 1,2,3,4 ... 
	# build dataframe that includes orig_order and file_time(runorder)
	orig_order <- 1:nrow(inputdata)
	
	if (is.null(measurementOrder))
	{																
		file_time <- sampleList$measurementTimeStamp											
	} else
	{
		file_time <- measurementOrder
	}
	
	file_time_order <- order(file_time)
	orig_order <- order(file_time_order)
	
	# order resulting dataframe according to Runorder = funcdata
	funcdata 	<- inputdata [file_time_order, ]
	sampClass <- PeakMLdata$phenoData[file_time_order]								
	pbqc_rows <- which(sampClass==pbqc_groupid)

	

	if (!is.null(minIntThrSamples))
	{
		minIntSamp <- which(sampClass%in%minIntThrSamples)
	} else
	{
		minIntSamp <- NULL
	}

	if (scatter.plots==TRUE)
	{
		pdf("QC_RLSC_output.pdf")
		Masses <- apply(PeakTable[[2]], 2, mean, na.rm=TRUE)
		RTs <- apply(PeakTable[[3]], 2, mean, na.rm=TRUE)
		RTminutes <- floor(RTs/60)
		RTseconds <- floor(RTs-(RTminutes*60))
		titles <- paste ("Mass: ", round(Masses,5), ", RT: ",RTminutes," min ",RTseconds," s", sep="")
	}

	QC.RLSC.Correction.Factor <- vector("list", ncol(funcdata))

	for (jj in 1:ncol(funcdata)) 										 							
	{
		#### Predict loess curve from pbqc data in each column 'jj'
		pbqc_mvals <- funcdata[pbqc_rows, jj]
		NAs <- which(is.na(pbqc_mvals))
		if (length(NAs)==0)
		{			  		 		 					
			pbqc.lowess <- lowess(x=1:length(pbqc_mvals), y=pbqc_mvals, f=f, iter=3)
		} else
		{
			pbqc.lowess <- lowess(x=c(1:length(pbqc_mvals))[-c(NAs)], y=pbqc_mvals[-c(NAs)], f=f, iter=10)
			newvals <- rep(NA_real_,length(pbqc_mvals))
			newvals[-c(NAs)] <- pbqc.lowess$y
			pbqc.lowess$y <- newvals
		}
		
		# Interpolate LOWESS fit to the total number of samples, using not-so-CV cubic spline function.
		spline.input <- rep(NA_real_,nrow(funcdata))
		spline.input[pbqc_rows] <- pbqc.lowess$y
		
		naCount <- length(which(is.na(pbqc_mvals)))
		
		if ((length(pbqc_mvals) - naCount)>5)
		{
			dat <- data.frame(x=c(1:length(spline.input)), y=spline.input)
			smoo <- with(dat[!is.na(dat$y),],smooth.spline(x,y))
			result <- with(dat,predict(smoo,x[is.na(y)]))
			dat[is.na(dat$y),] <- result
			
			correction.factor <- max(dat$y)-dat$y
		} else
		{
			correction.factor <- rep(0,nrow(funcdata))
			cat ("Warning: This peak was detected in less than 5 QC samples. You may consider using more stringent QC RSD filter settings.", "\n")
		}		
		
		if (length(minIntSamp)>0 & max(correction.factor, na.rm=TRUE)!=0 )
		{
			if (!all(is.na(funcdata[minIntSamp,jj])))
			{
				maxInt <- max(funcdata[minIntSamp,jj], na.rm=TRUE)
			} else
			{
				maxInt <- NA
			}
			if (!is.na(maxInt))
			{
				hits <- which(funcdata[,jj]<=maxInt)
				correction.factor[hits] <- 0
			}
		}
		
		## For peaks with no detection replace correction.factor with NA
		NAs <- which(is.na(funcdata[,jj]))
		if (length(NAs)>0)
		{
			correction.factor[NAs] <- NA
		}
				
		if (scatter.plots==TRUE)
		{
			layout (matrix(c(1,2),1,2,byrow=TRUE), widths=c(4,1.3))
			cols <- rep (1,nrow(funcdata))
			cols[pbqc_rows] <- 2
			pchars <- rep(1,nrow(funcdata))
			pchars[pbqc_rows] <- 2
			pchars2 <- rep(21,nrow(funcdata))
			pchars2[pbqc_rows] <- 24
			
			if (length(minIntSamp)>0)
			{
				cols[minIntSamp] <- 3
			}
			YMAX <- max(c(funcdata[,jj],funcdata[,jj]+correction.factor), na.rm=TRUE)
			YMIN <- min(c(funcdata[,jj],funcdata[,jj]+correction.factor), na.rm=TRUE)
			par(mar=c(3,3,3,0.5))
			plot (funcdata[,jj], col=cols, pch=pchars, type="p", ylab="", xlab="Run order", main=titles[jj], ylim=c(YMIN, YMAX))
			points (dat$y, type="l", col=4)
			if (max(correction.factor, na.rm=TRUE)!=0)
			{
				points (funcdata[,jj]+correction.factor, type="p", pch=pchars2, bg=cols)
			}
			par(mar=c(0,0,3,0))
			plot (1,1,pch="",xlab="",ylab="", axes=F)
			legend ("topleft",legend=c("Sample, uncorr.","QC, uncorr.","Blank, uncorr.","Sample, corr.","QC, corr.","Blank. corr."), pch=c(1,2,1,21,24,21), pt.bg=c(1,2,3,1,2,3), col=c(1,2,3,1,1,1), ncol=1)
			
		}	
		
		## Reorded correction.factor back to original order of the files in PeakML file and write it back to PeakML. Write groupa annotation as well, just in case
		correction.factor <- correction.factor[orig_order]
		
		
		## And now mess up with chromatogram data. 
		selchroms <- which(PeakMLdata$peakDataMtx[,10]==jj)
		correction.factor <- correction.factor[!is.na(correction.factor)]
		
		QC.RLSC.Correction.Factor[[jj]] <- paste(round(correction.factor,3),collapse=", ")
		
		if (length(selchroms)!=length(correction.factor))
		{
		 cat ("Oh no! Something unexpected happened! Please contact, Andris Jankevics (andris.jankevics@manchester.ac.uk) and report this error.","\n")
		 stop ()
		}
		
		for (chromn in 1:length(selchroms))
		{
			chromdata <- PeakMLdata$chromDataList[[selchroms[chromn]]]
			#maxint <- max(chromdata[2,],na.rm=TRUE)
			#chromdata[2,chromdata[2,]==maxint] <- maxint+correction.factor[chromn]
			chromdata[2,] <- chromdata[2,]+correction.factor[chromn]
			PeakMLdata$chromDataList[[selchroms[chromn]]] <- chromdata
		}
						 	 		 					
	}
	
	if (scatter.plots==TRUE)
	{
		dev.off ()
	}
	
	# Write PeakML output
	PeakMLdata$GroupAnnotations$QC.RLSC.Correction.Factor <- QC.RLSC.Correction.Factor
	PeakML.Write (peakMLdata=PeakMLdata, outFileName=outputfile)
}
