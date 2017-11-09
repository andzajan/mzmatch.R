PeakML.Plot.Chromatograms <- function(PeakMLdata,groupid,sampleClasses=NULL,xaxis=TRUE, Title=NULL, mar=c(2,1,0,1))
{
	MAXf <- function (x,i) {max(PeakMLdata$chromDataList[[x]][i,])}
	MINf <- function (x,i) {min(PeakMLdata$chromDataList[[x]][i,])}

	which.chromatograms <- which(PeakMLdata$peakDataMtx[,10]==groupid)
	plot.colors <- PeakMLdata$phenoData
	## plot only chromatograms which match selected samples
	if (!is.null(sampleClasses))
	{
		plot.colors[!plot.colors%in%sampleClasses] <- NA
	}
	plot.colors <- as.factor(plot.colors)
	plot.legend <- levels(plot.colors)
	plot.colors <- as.numeric(plot.colors)

	nas <- which(is.na(plot.colors[PeakMLdata$peakDataMtx[which.chromatograms,9]]))
	if (length(nas)>0)
	{
		which.chromatograms <- which.chromatograms[-c(nas)]
	}
	if (length(which.chromatograms)>0)
	{
		plot.colors <- plot.colors[PeakMLdata$peakDataMtx[which.chromatograms,9]]
		YMAX <- max(sapply(which.chromatograms,MAXf,i=2))
		XMIN <- min(sapply(which.chromatograms,MINf,i=3))
		XMAX <- max(sapply(which.chromatograms,MAXf,i=3))
		par (mar=mar)
		plot (1,1,pch="",xlab="",ylab="",ylim=c(0,YMAX),xlim=c(XMIN,XMAX),font.lab=2, axes=F)
		axis(2)
		if(xaxis==TRUE)
		{
			axis(1)
		}
		box()
		legend ("topright",fill=c(1:length(plot.legend)),plot.legend,cex=0.7)
		for (i in 1:length(which.chromatograms))
		{
			points(PeakMLdata$chromDataList[[which.chromatograms[i]]][3,],PeakMLdata$chromDataList[[which.chromatograms[i]]][2,],col=plot.colors[i],type="l")
		}
		if (!is.null(Title))
		{
			title(main=Title)
		}
	} else
	{
		plot(1,1,pch="",axes=F,xlab="",ylab="")
	}
}
