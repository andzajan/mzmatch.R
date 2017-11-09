PeakML.trendPlot <- function(filename,sampleclasses,samplenames=NULL,outputfile,xaxisvalues=NULL,xaxislabels=NULL,log2scale=FALSE,collapse=FALSE,plotchromatograms=FALSE,xaxislabel=NULL,percentiles=c(15,20,35),classaverage=FALSE,legendposition="bottomright",yminvalue=NULL,xmaxvalue=NULL,colors= c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"),fontfamily="Helvetica",pdflayout=c(3,1))
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

	
	
	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")

	## Read samplenames from peakml file	
	sampnames <- .jcall(project,returnSig="[S", method="getMeasurementNames")
	## Measurementd ID's
	measurementids <- c(1:length(sampnames))-1	

	## Get masschromatograms
	cat ("Loading masschromatograms from file \n")		
	st <- system.time(chromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
	cat ("Done in:",st[3],"s \n")
	## Get setnums out of masschromatograms array
	peakdata <- NULL	
	for (chromnums in 1:length(chromatograms))
	{
		VEC <- .jevalArray(chromatograms[[chromnums]])[c(9,11,13,8)]
		peakdata <- rbind (peakdata,VEC)
	}

	## How many groups (sets are there)	
	setsnumber <- .jcall(project,returnSig="I",method="getNrPeaksets") 

	## Reorder sample names in peakML set according to ordert which should be plotted in trend plot
	sampleclasses <- as.factor(sampleclasses)	
	classlabels <- as.character(unique(sampleclasses))
	classcolors <- as.numeric(unique(sampleclasses))	
	classcolors <- colors[classcolors]

	reorderedsampleindexes <- NULL	
	reorderedsampleclasses <- NULL
	reorderedxaxisvalues <- NULL
	reorderedxaxislabels <- NULL
	
	if (is.null(samplenames))	
	{
		reorderedsampleindexes <- measurementids	
		reorderedsampleclasses <- as.numeric(sampleclasses)
		reorderedxaxisvalues <- xaxisvalues
		reorderedxaxislabels <- xaxislabels

	} else
	{
		for (sampleid in 1:length(samplenames))
		{
			NUM <- which (sampnames==samplenames[sampleid])	
			reorderedsampleindexes <- append(reorderedsampleindexes,measurementids[NUM])
			reorderedsampleclasses <- append(reorderedsampleclasses,sampleclasses[sampleid])	
			reorderedxaxisvalues <- append(reorderedxaxisvalues,xaxisvalues[sampleid])
			reorderedxaxislabels <- append(reorderedxaxislabels,xaxislabels[sampleid])
		}
	}

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
	plotnames <- rep(NA,setsnumber)

	for (setnum in 1:setsnumber)
	{
		cat (setnum, ", ")		
		settable <- peakdata[peakdata[,3]==(setnum-1),]
		settable <- rbind(settable,NULL)		
		outtable <- NULL
		plotcolors <- NULL	
		xvals <- NULL	
		for (index in 1:length(reorderedsampleindexes))
		{
			SEL <- settable[settable[,2]==reorderedsampleindexes[index],]
			#cat (length(SEL))			
			if (length(SEL)==0) 
			{
				SEL <- c(NA,reorderedsampleindexes[index],setnum-1,NA)			
				}
			SEL <- rbind(SEL,NULL)
			#cat (" ",nrow(SEL),index,"\n")
			plotcolors <- append(plotcolors,rep(reorderedsampleclasses[index],nrow(SEL)))
			xvals <- append(xvals,rep(reorderedxaxisvalues[index],nrow(SEL)))
			outtable <- rbind(outtable,SEL)		
		}

		plotname <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(setnum-1),"identification")
		if (is.null(plotname))
		{
			plotname <- round(mean(settable[,4],na.rm=TRUE),5)
		}
		plotnames[setnum] <- plotname
		if (is.null(reorderedxaxisvalues))
		{
			xaxis <- 1:nrow(outtable)
		} else
		{
			xaxis <- xvals
		}
		if (log2scale==TRUE)
		{
			outtable[,1] <- log2(outtable[,1])
		} 

		intavg <- mean(outtable[,1],na.rm=TRUE)	
		if (all(is.na(outtable[,1])))
		{
			intavg <- 1000
			outtable[,1] <- 0
		}
	
		if (!is.null(percentiles))		
		{
			ymax <- intavg+((intavg*max(percentiles))/100)
			ymin <- intavg-((intavg*max(percentiles))/100)
			yaxislimits <- c(min(ymin,min(outtable[,1],na.rm=TRUE)),max(ymax,max(outtable[,1],na.rm=TRUE)))
		} else
		{
			yaxislimits <- c(min(outtable[,1],na.rm=TRUE),max(outtable[,1],na.rm=TRUE))
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
				VAL <- outtable[plotcolors==unique(plotcolors)[sclass],1]
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
			plot (xaxis,outtable[,1],type="p",pch="",ylim=yaxislimits,xlab=xaxislabel,ylab="intensity",xlim=c((min(xaxis)-0.5),xlimmax),axes=F)	
			for (i in 1:length(unique(oldgroups)))	
			{
				if (!is.null(reorderedxaxisvalues))
				{
					xaxis2 <- reorderedxaxisvalues[which(oldgroups==unique(oldgroups)[i])]
				} else
				{
					xaxis2 <- xaxis[which(oldgroups==unique(oldgroups)[i])]
				}		
				points(xaxis2,outtable[which(oldgroups==unique(oldgroups)[i]),1],type="b",pch=21,cex=1.2,bg=classcolors[i])			
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
					abline (h=mean(outtable[,1],na.rm=TRUE),col=2,lty=2)
					M <- ((mean(outtable[,1],na.rm=TRUE)*step)/100)
					yvalue <- mean(outtable[,1],na.rm=TRUE)
					abline (h=mean(outtable[,1],na.rm=TRUE)+M,col="gray50",lty=2)	
					abline (h=mean(outtable[,1],na.rm=TRUE)-M,col="gray50",lty=2)
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
					abline (h=mean(outtable[,1],na.rm=TRUE),col=2,lty=2)					
					M <- ((mean(outtable[,1],na.rm=TRUE)*step)/100)
					yvalue <- mean(outtable[,1],na.rm=TRUE)
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
				points(xaxis,outtable[plotcolors==unique(plotcolors)[i],1],type="b",pch=21,cex=1.2,bg=unique(plotcolors)[i])			
			}
			title (plotname)
		}
		legend (legendposition,pch=21,pt.bg=classcolors,classlabels,cex=0.9)
	}
	dev.off ()

	if (plotchromatograms==TRUE)
	{
		pdf (file=paste(outputfile,"_masschromatograms.pdf",sep=""),paper="A4",width=8,height=11.4,family=fontfamily)	
		par (mfrow=pdflayout)	
		for (setnum in 1:setsnumber)
		{					
			#plotname <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(setnum-1),"identification")
			#if (is.null(plotname))
			#{
			#	plotname <- round(mean(settable[,4],na.rm=TRUE),5)
			#}
			chromnums <- which(peakdata[,3]==setnum-1)
			plotcolors <- peakdata[chromnums,2]
			plotcols <- NULL
			for (coli in 1:length(plotcolors))
			{
				plotcols <- append(plotcols,reorderedsampleclasses[reorderedsampleindexes==plotcolors[coli]])
		
			}		
			chromatogramslist <- vector ("list",length(chromnums))
			maxrt <- NULL
			minrt <- NULL
			maxint <- NULL
			minint <- NULL
			for (chrnum in 1:length(chromnums))	
			{
				retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(chromnums[chrnum]-1))
				intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(chromnums[chrnum]-1))
				chromatogramslist[[chrnum]] <- rbind(retentiontimes,intensities)
				maxrt <- append(maxrt,max(retentiontimes))
				minrt <- append(minrt,min(retentiontimes))
				maxint <- append(maxint,max(intensities))
				minint<- append(minint,min(intensities))			
			}
		plot (1,1,xlim=c(min(minrt),max(maxrt)),ylim=c(0,max(maxint)),main=plotnames[setnum],xlab="rt",ylab="intensity")		
		for  (i in 1:length(chromatogramslist))
		{
			points (chromatogramslist[[i]][1,],chromatogramslist[[i]][2,],type="l",col=plotcols[i])
		}
		legend ("topright",pch=21,pt.bg=classcolors,classlabels,cex=1)
		}	
		dev.off ()
	}
}

