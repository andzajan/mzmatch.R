PeakML.Plot.RankProducts <- function (pdffile, PeakMLtree, RankProducts, identifications, sortby="Class1")
{
	pdf (pdffile)
	if (is.null(ncol(RankProducts$RPs)))
	{
		ranks <- RankProducts$RPs
	} else
	{
		if (sortby=="Class1")
		{
			ranks <- RankProducts$RPs[,1]
		} else
		{
			ranks <- RankProducts$RPs[,2]
		}
	}

	if (sortby=="Class1")
	{
		pfp <- RankProducts$erp[,1]
	} else
	{
		pfp <- RankProducts$erp[,2]
	}

	orderedranks <- sort(pfp[1:length(pfp)])
	rankindex <- c(1:length(pfp))[order(pfp[1:length(pfp)])]

	outputtable <- data.frame(nrow=length(rankindex),ncol=8)

	erps <- RankProducts$erp[rankindex,]
	pfps <- RankProducts$pfp[rankindex,]

	intensities.Class1 <- rep(NA,length(rankindex))
	intensities.Class2 <- rep(NA,length(rankindex))
	keggids <- rep(NA,length(rankindex))
	idents <- rep(NA,length(rankindex))
	measmass <- rep(NA,length(rankindex))

	for (peakn in 1:length(rankindex))
	{
		num <- rankindex[peakn]
		rownums <- PeakMLtree$peakDataMtx[PeakMLtree$peakDataMtx[,10]==num,]
		rowhits <- which(PeakMLtree$peakDataMtx[,10]==num)

		## Select which chromatograms to plot
		selectedrows <- which(rownums[,9]%in%RankProducts$sampleindex)
		if (length(selectedrows)>0)
		{
			selectedrowhits <- rowhits[selectedrows]
			plotcolors <- as.factor(PeakMLtree$phenoData[rownums[selectedrows,9]])
			Ints <- vector("list")
			RTs <- vector("list")

			for (rwh in 1:length(selectedrowhits))
			{
				Ints[[rwh]] <- PeakMLtree$chromDataList[[selectedrowhits[rwh]]][2,]
				RTs[[rwh]] <- PeakMLtree$chromDataList[[selectedrowhits[rwh]]][3,]
			}

			YMAX <- max(unlist(Ints))
			XMIN <- min(unlist(RTs))
			XMAX <- max(unlist(RTs))

			layout (matrix (nrow=2,ncol=2,c(1,2,1,3)), widths=c(0.4,0.6))
			plot (1,1,pch="",xlab="RT",ylab="Intensity",xlim=c(XMIN,XMAX),ylim=c(0,YMAX))
			for (rwh in 1:length(selectedrowhits))
			{
				points (RTs[[rwh]],Ints[[rwh]],type="l",col=plotcolors[rwh])
			}
			if (all(is.na(identifications[[rankindex[peakn]]])))
			{
				subt <- NA
			} else
			{
				subt <- paste(identifications[[rankindex[peakn]]][,4],collapse=", ")
				idents[peakn] <- subt
				kegghits <- which(identifications[[rankindex[peakn]]][,5]=="kegg")
				if (length(kegghits)>0)
				{
					keggids[peakn] <- paste(identifications[[rankindex[peakn]]][kegghits,1],collapse=", ")
				}
			}
			title(main=paste("Rank=",round(ranks[rankindex[peakn]],0),", pfp_1=",round(RankProducts$pfp[rankindex[peakn],1],5),", pfp_2=",round(RankProducts$pfp[rankindex[peakn],2],5), ", Erp_1=",round(RankProducts$erp[rankindex[peakn],1],5), ", Erp_2=",round(RankProducts$erp[rankindex[peakn],2],5),", Mass:",round(median(rownums[selectedrows,1]-PeakMLtree$massCorrection[[1]]),5),sep=""), sub= subt, cex.sub=0.7)

			measmass[peakn] <- round(median(rownums[selectedrows,1]-PeakMLtree$massCorrection[[1]]),5)

			sampnums <- table(as.numeric(plotcolors))
			plot (1,1,pch="",xlab="",ylab="", axes=F)
			legend ("topleft", fill=c(1,2), c(paste(sampnums[1], "signals in Class 1"),paste(sampnums[2], "signals in Class 2")),cex=1)
			maxints <- rep(NA, length(Ints))
			for (setn in 1:length(Ints))
			{
				maxints[setn] <- max(Ints[[setn]], na.rm=TRUE)[1]
			}
			boxplot (maxints ~ plotcolors, col=c("gray50","red"))


			## Calculate intensities for peaks
			intensities.peaks <- PeakMLtree$peakDataMtx[selectedrowhits,8]
			hits <- which(as.numeric(plotcolors)==1)
			if (length(hits)>0)
			{
				intensities.Class1[peakn] <- median(intensities.peaks[hits], na.rm=TRUE)
			}
			hits <- which(as.numeric(plotcolors)==2)
			if (length(hits)>0)
			{
				intensities.Class2[peakn] <- median(intensities.peaks[hits], na.rm=TRUE)
			}

		} else
		{
			par(mfcol=c(1,1))
			plot (1,1,pch="")
			title(main=paste("Rank=",round(ranks[rankindex[peakn]],5),", pfp_1=",round(RankProducts$pfp[rankindex[peakn],1],5),", pfp_2=",round(RankProducts$pfp[rankindex[peakn],2],5),", Mass:",round(median(rownums[selectedrows,1]-PeakMLtree$massCorrection[[1]]),5),sep=""), sub= subt, cex.sub=0.7)
		}
	}
	dev.off ()

	outputtable <- data.frame(erps, pfps, intensities.Class1/intensities.Class2, intensities.Class2/intensities.Class1, measmass, idents, keggids)

	colnames(outputtable) <- c("Erp.Class1","Erp.Class2","pfp.Class1","pfp.Class2","ratio.Class1","ratio.Class2","measured.mass","putative.identifcation","KEGG.ids")
	write.table (row.names=FALSE, outputtable,file=paste(pdffile,".tsv",sep=""),sep="\t")

}
