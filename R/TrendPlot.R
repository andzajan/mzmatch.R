TrendPlot <- function(peaktable,sampleclasses,plotnames=NULL,samplenames=NULL,outputfile,xaxisvalues=NULL,xaxislabels=NULL,log2scale=FALSE,collapse=FALSE,xaxislabel=NULL,percentiles=c(15,20,35),classaverage=FALSE,legendposition="bottomright",yminvalue=NULL,xmaxvalue=NULL,colors= c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"),fontfamily="Helvetica",pdflayout=c(3,1))
{
	# Check if classcolors are given as numeric vector, convert it to character definitions	
	if (is.numeric(colors))
	{
		colors <- colors()[colors]
	}
	
	if(length(which(colors()%in%colors))<length(colors))
	{
		cat ("One or more color names which you defined in parameter \"classcolors\" are note defined in R color palette. Please, check the color names again. \n")
	}

	if (is.null(samplenames))
	{
		sampnames <- rownames(peaktable) 
	} else
	{
		sampnames <- samplenames
	}
	
	## Measurementd ID's
	measurementids <- c(1:length(sampnames))-1	

	## How many groups (sets are there)	
	setsnumber <- ncol(peaktable)

	## Reorder sample names in peakML set according to ordert which should be plotted in trend plot
	sampleclasses <- as.factor(sampleclasses)	
	classlabels <- as.character(unique(sampleclasses))
	classcolors <- as.numeric(unique(sampleclasses))	
	classcolors <- colors[classcolors]

	reorderedsampleindexes <- NULL	
	reorderedsampleclasses <- NULL
	reorderedxaxisvalues <- NULL
	reorderedxaxislabels <- NULL
	
	reorderedsampleindexes <- measurementids	
	reorderedsampleclasses <- as.numeric(sampleclasses)
	reorderedxaxisvalues <- xaxisvalues
	reorderedxaxislabels <- xaxislabels

	## Replace classvalues with actual plotting color
	newsampleclasses <- reorderedsampleclasses
	for (i in 1:length(unique(reorderedsampleclasses)))
	{
		newsampleclasses[reorderedsampleclasses==unique(reorderedsampleclasses)[i]] <- classcolors[i]
	}
	oldgroups <- reorderedsampleclasses	
	reorderedsampleclasses <- newsampleclasses
	
	pdf (file=paste(outputfile,".pdf",sep=""),paper="A4",width=8,height=11.4,family=fontfamily)	
	par (mfrow=pdflayout)

	for (setnum in 1:setsnumber)
	{
		cat (setnum, ", ")		
		outtable <- peaktable[,setnum]
		plotcolors <- reorderedsampleclasses
		xvals <- reorderedxaxisvalues
		if (!is.null(plotnames))
		{
			plotname <- plotnames[setnum]
		} else
		{
			plotname <- colnames(peaktable)[setnum]
		}
		
		if (is.null(reorderedxaxisvalues))
		{
			xaxis <- 1:length(outtable)
		} else
		{
			xaxis <- xvals
		}
		if (log2scale==TRUE)
		{
			outtable <- log2(outtable)
		} 

		intavg <- mean(outtable,na.rm=TRUE)	
		if (all(is.na(outtable)))
		{
			intavg <- 1000
			outtable <- rep(0,length(reorderedsampleclasses))
		}
	
		if (!is.null(percentiles))		
		{
			ymax <- intavg+((intavg*max(percentiles))/100)
			ymin <- intavg-((intavg*max(percentiles))/100)
			yaxislimits <- c(min(ymin,min(outtable,na.rm=TRUE)),max(ymax,max(outtable,na.rm=TRUE)))
		} else
		{
			yaxislimits <- c(min(outtable,na.rm=TRUE),max(outtable,na.rm=TRUE))
		}

		if (!is.null(yminvalue))
		{
			yaxislimits <- c(yminvalue,yaxislimits[2])
		}
		
		if (classaverage==TRUE)
		{
						
			plottable <- NULL
			for (sclass in 1:length(unique(plotcolors)))	
			{
				VAL <- outtable[plotcolors==unique(plotcolors)[sclass]]
				plottable <- rbind(plottable,c(mean(VAL,na.rm=TRUE),sqrt(var(VAL,na.rm=TRUE))))
			}		
			
			PCH <- rep(21,nrow(plottable))
			PCH[is.na(plottable[,1])] <- 3	
			plottable[is.na(plottable[,1]),1] <- 0	
			plottable[is.na(plottable[,2]),2] <- 0
			
			yaxislimits <- c(min(plottable[,1]-plottable[,2],na.rm=TRUE),max(plottable[,1]+plottable[,2],na.rm=TRUE))

			if (!is.null(percentiles))		
			{
				intavg <- mean(plottable[,1],na.rm=TRUE)				
				ymax <- intavg+((intavg*max(percentiles))/100)
				ymin <- intavg-((intavg*max(percentiles))/100)
				yaxislimits <- c(min(ymin,yaxislimits[1]),max(ymax,yaxislimits[2]))
			}
			
			if (!is.null(yminvalue))
			{
				yaxislimits <- c(yminvalue,yaxislimits[2])
			}
				
			if (is.null(reorderedxaxisvalues))	
			{
				xaxis <- c(1:nrow(plottable))	
			} else
			{
				xaxis <- c(min(reorderedxaxisvalues),max(reorderedxaxisvalues))			
			}	

			plot (plottable[,1],type="b",pch=PCH,bg=unique(plotcolors),cex=1.2,ylim=yaxislimits,xlab=xaxislabel,ylab="intensity",xlim=c((min(xaxis)-0.5),(max(xaxis)+0.5)),col="gray50")	
			grid ()

			for (sclass in 1:nrow(plottable))
			{
				arrows (sclass,plottable[sclass,1]+plottable[sclass,2],sclass,plottable[sclass,1]-plottable[sclass,2],length=0.05,code=3,angle=90)			
			}
			points (plottable[,1],type="b",pch=PCH,bg=unique(plotcolors),cex=1.2,col="gray50")			
			title (plotname)
			if (!is.null(percentiles))	
			{		
				for (step in percentiles)
				{
					abline (h=mean(plottable[,1],na.rm=TRUE),col=2,lty=2)
					M <- ((mean(plottable[,1],na.rm=TRUE)*step)/100)
					yvalue <- mean(plottable[,1],na.rm=TRUE)
					abline (h=mean(plottable[,1],na.rm=TRUE)+M,col="gray50",lty=2)	
					abline (h=mean(plottable[,1],na.rm=TRUE)-M,col="gray50",lty=2)
					text (0.5,yvalue+M,paste(step,"%",sep=""),font=2)
					text (0.5,yvalue-M,paste(step,"%",sep=""),font=2)
				}
			}
		}

		if (collapse==FALSE & classaverage==FALSE)
		{		
			if (!is.null(xmaxvalue))
			{
				xlimmax <- xmaxvalue	
			} else
			{
				xlimmax <- (max(xaxis)+0.5)
			}	
			cat (xlimmax, "\n")	
			cat (outtable,"\n")
			cat (xaxis,"\n")		
			plot (xaxis,outtable,type="p",pch="",ylim=yaxislimits,xlab=xaxislabel,ylab="intensity",xlim=c((min(xaxis)-0.5),xlimmax),axes=F)	
			for (i in 1:length(unique(oldgroups)))	
			{
				if (!is.null(reorderedxaxisvalues))
				{
					xaxis2 <- reorderedxaxisvalues[which(oldgroups==unique(oldgroups)[i])]
				} else
				{
					xaxis2 <- xaxis[which(oldgroups==unique(oldgroups)[i])]
				}		
				points(xaxis2,outtable[which(oldgroups==unique(oldgroups)[i])],type="b",pch=21,cex=1.2,bg=classcolors[i])			
			}
			title (plotname)
			axis (2)
			if (!is.null(reorderedxaxislabels) & !is.null(reorderedxaxisvalues))
			{
				axis(1,at=reorderedxaxisvalues,labels=as.character(reorderedxaxislabels))
			} else
			{
				axis (1)
			}
			box ()
			
			if (!is.null(percentiles))	
			{		
				for (step in percentiles)
				{
					abline (h=mean(outtable,na.rm=TRUE),col=2,lty=2)
					M <- ((mean(outtable,na.rm=TRUE)*step)/100)
					yvalue <- mean(outtable,na.rm=TRUE)
					abline (h=mean(outtable,na.rm=TRUE)+M,col="gray50",lty=2)	
					abline (h=mean(outtable,na.rm=TRUE)-M,col="gray50",lty=2)
					text (min (xaxis)-0.5,yvalue+M,paste(step,"%",sep=""),font=2)
					text (min (xaxis)-0.5,yvalue-M,paste(step,"%",sep=""),font=2)
				}
			}
		}

		if (collapse==TRUE & classaverage==FALSE)
		{	
			if (!is.null(reorderedxaxisvalues))
			{	
				xaxis <- min(reorderedxaxisvalues):max(reorderedxaxisvalues)
			} else
			{		
				xaxis <- 1: max(table(plotcolors))
			}
	
			plot (1,1,type="p",pch="",xlab=xaxislabel,ylab="intensity",xlim=c((min(xaxis)-0.5),(max(xaxis)+0.5)),ylim=yaxislimits)
			if (!is.null(percentiles))	
			{			
				for (step in percentiles)
				{
					abline (h=mean(outtable,na.rm=TRUE),col=2,lty=2)					
					M <- ((mean(outtable,na.rm=TRUE)*step)/100)
					yvalue <- mean(outtable,na.rm=TRUE)
					abline (h=yvalue+M,col="gray50",lty=2)	
					abline (h=yvalue-M,col="gray50",lty=2)
					text (min(xaxis)-0.5,yvalue+M,paste(step,"%",sep=""),font=2)
					text (min (xaxis)-0.5,yvalue-M,paste(step,"%",sep=""),font=2)
				}
			}
			for (i in 1:length(unique(plotcolors)))	
			{
				if (!is.null(reorderedxaxisvalues))
				{
					xaxis <- reorderedxaxisvalues[plotcolors==unique(plotcolors)[i]]
				}			
				points(xaxis,outtable[plotcolors==unique(plotcolors)[i]],type="b",pch=21,cex=1.2,bg=unique(plotcolors)[i])			
			}
			title (plotname)
		}
		legend (legendposition,pch=21,pt.bg=classcolors,classlabels,cex=0.9)
	}
	dev.off ()
}

