PeakML.Plot.RelatedPeaks <- function (filename, ionisation="detect", Rawpath=NULL, DBS=dir(paste(find.package("mzmatch.R"), "/dbs", sep=""),full.names=TRUE),outputfile,sampleClasses=NULL)
{
	plotPeak <- function (peakn)
	{
		peakmass <- round(Masses[which.peaksets[peakn]]-(PeakMLdata$massCorrection[[1]]),7)
		peakRT <- round(RTs[which.peaksets[peakn]],0)
		peakInt <- format(round(Intensities[which.peaksets[peakn]],0),scientific=TRUE)
		adduct <- PeakMLdata$GroupAnnotations$relation.ship[which.peaksets[peakn]]
		infmat <- matrix(c(peakmass,peakRT,peakInt,adduct),ncol=1,nrow=4)
		rownames(infmat) <- c("m/z","RT, s","Intensity","Derivative")
		textplot(infmat,halign="left",valign="top",mar=c(0,0,0,0),show.colnames=FALSE)
		idents <- id.resolved[[which.peaksets[peakn]]]

		if (!all(is.na(idents)))
		{
			idents[,6] <- round(as.numeric(idents[,6]),1)
			colnames (idents) <- c("id","Formula","Mass","Name","DB","ppm","Addcut")
			textplot(idents[,c(2,7,6,4,5)],halign="left",valign="top",mar=c(0,0,0,0.5),show.rownames=FALSE)
		} else
		{
			plot(1,1,xlab="",ylab="",pch="",axes=FALSE)
		}
		groupid <- which.peaksets[peakn]
		PeakML.Plot.Chromatograms(PeakMLdata=PeakMLdata,groupid=groupid,sampleClasses=sampleClasses,xaxis=FALSE)
	}


	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	id.resolved <- PeakML.Methods.DBidToCompoundName(DBS=DBS, PeakMLdata=PeakMLdata, collapse=FALSE)
	Masses <- apply(PeakTable[[2]],2,median,na.rm=TRUE)
	RTs <- apply(PeakTable[[3]],2,median,na.rm=TRUE)
	Intensities <- apply(PeakTable[[1]],2,max,na.rm=TRUE)

	## How many cluster there are

	peak.clusters <- unique(PeakMLdata$GroupAnnotations$relation.id)

	pdf (file=paste(outputfile,".pdf",sep=""),paper="A4",width=8,height=11.4)
	for (pclust in 1:length(peak.clusters))
	{
		which.peaksets <- which(PeakMLdata$GroupAnnotations$relation.id==peak.clusters[pclust])

		col.vector <- rep(NA,length(which.peaksets))
		col.vector[is.na(id.resolved[which.peaksets])] <- 1
		col.vector[is.na(col.vector)] <- 2
		layout (matrix(c(1,1,1,2:16),ncol=3,nrow=6,byrow=TRUE),widths=c(0.25,0.55,0.2))
		par (mar=c(4,4,1,0))
		ymax <- max(Intensities[which.peaksets])+(max(Intensities[which.peaksets])/10)
		plot (Masses[which.peaksets]-(PeakMLdata$massCorrection[[1]]),Intensities[which.peaksets],type="h",xlab="m/z",ylab="Intensity",lwd=1.5,col=col.vector,ylim=c(0,ymax),font.lab=2)
		legend ("topright", fill=c(2),c("putatively annotated"))
		text (Masses[which.peaksets]-(PeakMLdata$massCorrection[[1]]),Intensities[which.peaksets],round(Masses[which.peaksets]-(PeakMLdata$massCorrection[[1]]),5),pos=3,cex=0.9,font=2,col=col.vector)
		RTtime <- round(median(RTs[which.peaksets]),0)
		RTtimemin <- floor(RTtime/60)
		RTtimemin <- c(RTtimemin,RTtime-(RTtimemin*60))
		mtext (paste("RT: ",RTtimemin[1],"min ",RTtimemin[2],"s (",RTtime,"s)",sep=""),side=3,cex=0.9)
		if (length(which.peaksets)<=5)
		{
			maxnum <- length(which.peaksets)
		} else
		{
			maxnum <- 5
		}
		for (peakn in 1:maxnum)
		{
			plotPeak(peakn)
		}
		if (length(which.peaksets)>5)
		{
			# how many extra pages to create
			peak.dif <- length(which.peaksets)-5
			number.of.pages <- ceiling(peak.dif/6)
			for (pgn in 1:number.of.pages)
		{
				maxnum <- c((pgn*6):((pgn*6)+5))
				if (max(maxnum)>length(which.peaksets))
				{
					maxnum <- c(maxnum[1],length(which.peaksets))
					maxnum <- unique(maxnum)
				}

				layout (matrix(c(1:18),ncol=3,nrow=6,byrow=TRUE),widths=c(0.25,0.55,0.2))
				for (peakn in maxnum)
				{
					plotPeak(peakn)
				}
			}
		}
	}
	dev.off ()
}



