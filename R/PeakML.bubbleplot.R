PeakML.bubbleplot <- function(xcmsset,peakidentifications,samples=NULL,samples2=NULL,DB1string=NULL,bubbletype="filled",plottype="identified",
outputfile="bubbleplot.pdf",cor.threshold=0.7,rtrange=NULL,mzrange=NULL,intensity.threshold.min=NULL,intensity.threshold.max=NULL, strict.relation.select=TRUE,scale="max",logbase=10,integral.val="maxint",bubblesize=3,smalllegend=FALSE)
{
	## Description of function parameters
	# xcmsset <- name of objet of class xcmsSet
	# peakidentifications <- matrix of identifications and peak relation ids, number of rows must be equal to number of groups in xcmsset@groups table
	# samples <- names or numbers for which bubble plots must be plotted, output is saved in pdf file, one plot per page. Default value 'NULL', all samples in set are plotted
	# samples2 <- name or numebr of sample with which to compare if plottype="compare"	
	# DB1string <- if this parameter is set, peaks which are identified with this DB tag are colored in red, string must contain unique name which is used for identifications
	# bubbletype <- "filled" (plot closed circles) or "open" for open circles on plot
	# plottype <- "identified" = bubble plot for identified compounds, "correlation" = correlation plot for dilution series, "compare" = compare two samples 
	# cor.threshold - corelation value treshold for no correlation
	# rtrange - retentiom range as vector c(rtmin,rtmax). Default value is set to 'NULL', full rtrange is plotted.
	# mazrange - m/z range as vector c(mzmin,mzax). Defaukt value 'NULL', full zmrange is plotted.
	# intensity.threshold.min - intensity level at which peaks are not plotted (min)
	# intensity.threshold.max - intensity level at which peaks are not plotted (max)
	# strict.relation.select - restrict selection of related peaks based on short list, isotopes, common adducts.
	##
	# scale - scale intensities across samples, peak intensity is represented as a bubble size
	#		"maxlog" - intensity values from all samples and all peak are scaled within interval 1..2, and log transformation is applied afterwards
	#					log(intensity/max(intensity)+1,logbase)*bubblesize
	#		"max" - intensities are scaled to constant "max" value of intensities vector from all samples, no log tranformation
	#		"log" - log is applied on intensities+0.001, adjust bubblesizs parameter for optimal bubble size
	#		"sumlog" - intensities(bubble size) are scaled to constant "total sum+1" value of intensities vector from all samples, log2 is taken afterwards.
	#		bysamplemaxlog - each sample is scaled by: log2((samp.int/max(samp.int))+1.1)
	#		bysamplelog - 
	##
	# logbase - "base" definition for logarithmic scaling
	# integral.val - if "maxint", maximum intensity value of the peak is used as circle diameter, if "integral" - area under peak (integral) are used
	# bubblesize - number by which bubles will magnified in the plots

	# Plot only these related peaks which are defined in list
	strict.relation.DB <- c("2M+ACN+H", "2M+K", "2M+Na", "2M+Na|C13 isotope #1", "2M+NH3|C13 isotope #1", "acetonitrile+sodium", "bp|C13 isotope #1", "bp|C13 isotope #2", "bp|O18 isotope #1", "methanol", "minus h2o", "minus h2o|C13 isotope #1", "minus nh3", "bp|N15 isotope #1", "bp|N15 isotope #2", "plus h2o", "plus h2o|C13 isotope #1", "plus h2o|N15 isotope #1", "salt (Na<>K)", "salt (Na<>NH4)", "sodium-formate","potassium","formic acid","acetonitrile","2*acetonitrile","2M+ACN+Na","ammonium","2M+NH3","amonium-formate","water","salt","salt (K<>NH4)","bp")

		if (integral.val=="maxint")
		{
			integral.val="maxo"		
		}
		if (integral.val=="integral")
		{
			integral.val="into"		
		}
	## Color definitions, using R full color palette
	## Hex codes, digtis 8 and 9 indicates color transparency level, these are removed for plottype="open" and colors plotted without trasnparency
		# Identified peaks from DB1, RED or positively correlated in dilution series	
		RED <- "#8B000090"
		# Related peaks of identifications from DB1, PINK
		PINK <- "#FF82AB90"
		# Identified peaks, GREEN or negatively correlated in dilution series
		GREEN <- "#00640090"
		# Related peaks of identified compound
		LIGHTGREEN <- "#C0FF3E90"
		## Unidentified peaks or no correlated
		BLACK <- "#00000070"
		## related peaks of unidentified
		GRAY <- "#8C8C8C90"
		## Dark blue, contaminants
		## Ligth blue, related peaks or contaminants
		LIGHTBLUE <- "#1874CD90"
		BLUE <- "#00008090"
	

	if (bubbletype=="filled")
	{
		PCH=21
		LWD=1
	} else
	{
		PCH=1
		LWD=2
	}

	##rtrange,mzrange
	if (is.null(rtrange))
	{
		rtrange <- 	c(min(xcmsset@groups[,"rtmed"]),max(xcmsset@groups[,"rtmed"]))
	}
	
	if (is.null(mzrange))
	{
		mzrange <- c(min(xcmsset@groups[,"mzmed"]),max(xcmsset@groups[,"mzmed"]))
	}
	
	## RT axis from scenods to minutes
	rtrange <- rtrange/60

	## Remove these peak groups and peaks whitch doesn't match intensity threshold criteria, as XCMS has a information about features in 'groups' table, we should extract it and, to get peak indexes in peak table, afterwards we can remove 'bad' peaks from bouth groups and peaks tables
		

	## Index of which peak in peaktable belongs to which peak group
	groupnums <- NULL
	for (i in 1:length(xcmsset@groupidx))
	{
		out <- rep(i,length(xcmsset@groupidx[[i]]))
		groupnums <- append(groupnums,out)
	}

	##Peaktable for all samples
	peaktable <- xcmsset@peaks
	if (is.null(intensity.threshold.max))
	{
		intensity.threshold.max <- max(peaktable[,integral.val])	
	}
	if (is.null(intensity.threshold.min))
	{
		intensity.threshold.min <- min(peaktable[,integral.val])	
	}
	
	## Indices of all features in peaktables which match intensity criteria
	goodfeatures <- which(peaktable[,integral.val] >= intensity.threshold.min & peaktable[,integral.val] <= intensity.threshold.max)
	
	## Which groupnumbers these features match, subselec a groupnums vector.
	goodgroups <- unique(groupnums[goodfeatures])
	cat (length(goodgroups),"peak groups are selected which matched intensity threshold criteria. \n") 

	## Group table of these peaks
	groupstable <- xcmsset@groups[goodgroups,]

	## identifications
	identifications <- as.factor(peakidentifications[goodgroups,1])

	## relation.id's for these groups
	relationids <- as.numeric(peakidentifications[goodgroups,2])

	## relationships for these groups
	relationships <- as.factor(peakidentifications[goodgroups,3])
	
	## empty two column matrix, in which label for each group will be added, second columnd will indicate the related peaks to this peak
	grouplabels <- matrix(data=NA,nrow=length(goodgroups),ncol=2)

	## Select only these peaks in peaktable which match intensity criteria
	peaktable <- peaktable[goodfeatures,]
	samplenumber <- peaktable[,"sample"]

	## New groupnums vector, which relays to new peak table
	groupnums <- groupnums[goodfeatures]
	
	## Filter out only these related peak identifications which are really common, rest or relationids are set to -1. Peaks will be marked as unidentified.
	if (strict.relation.select==TRUE)
	{
		relships <- NULL		
		for (i in 1:length(strict.relation.DB))
		{
			relship <- strict.relation.DB[i]
			relships <- append(relships,which(relationships==relship))
		}
		## Assign unique relationid for these peaks which are marked as unrelated		
		if (length(relships)!=0)
		{
			relationids[-c(relships)] <- max(relationids)+(1:length(relationids[-c(relships)]))
		} else
		{
			relationids <- c(1:length(relationids))
		}
	}
	

	## Sample names
	if (is.null(samples))
	{
		samples <- rownames(xcmsset@phenoData)	
	}
	
	if (is.character(samples))
	{
		sampnumvector <- which(rownames(xcmsset@phenoData)%in%samples)
	} else
	{
		sampnumvector <- samples	
	}	
	cat (nrow(peaktable)," ",sampnumvector,"\n")
	## Select only these peaks from peaktable which are requires by function call "samples" option
	samplenumbers <- NULL
	for (i in 1:length(sampnumvector))
	{
		samplenumbers <- append(samplenumbers,which(peaktable[,"sample"]==sampnumvector[i]))	
	}
	
	peaktable <- peaktable[samplenumbers,]
	groupnums <- groupnums[samplenumbers]
	samplenumber <- samplenumber[samplenumbers]
	cat (nrow(peaktable),"",sampnumvector,"\n")
	
	## Scaling of circle size
	if (scale=="maxlog")
	{
		intensities <- peaktable[,integral.val]
		intensities <- log((intensities/max(intensities))+1.1,logbase)*bubblesize
		peaktable[,integral.val] <- intensities
	}

	if (scale=="max")
	{
		intensities <- peaktable[,integral.val]
		intensities <- (intensities/max(intensities))*bubblesize
		peaktable[,integral.val] <- intensities
	}

	if (scale=="log")
	{
		intensities <- peaktable[,integral.val]
		peaktable[,integral.val] <- log(intensities+0.01,logbase)*bubblesize
	}
	
	if (scale=="sumlog")
	{
		intensities <- peaktable[,integral.val]
		intensities <- (intensities/sum(intensities))+1
		peaktable[,integral.val] <- log(intensities,logbase)*bubblesize
	}

	if (scale=="bysamplemaxlog")
	{
		for (i in 1:length(sampnumvector))		
		{
			intensities <- peaktable[peaktable[,"sample"]==sampnumvector[i],integral.val]
			intensities <- log((intensities/max(intensities))+1.1,logbase)*bubblesize
			peaktable[peaktable[,"sample"]==sampnumvector[i],integral.val] <- intensities
		}
	}	

	if (scale=="bysamplelog")
	{
		for (i in 1:length(sampnumvector))		
		{
			intensities <- peaktable[peaktable[,"sample"]==sampnumvector[i],integral.val]
			intensities <- log(intensities,logbase)*bubblesize
			peaktable[peaktable[,"sample"]==sampnumvector[i],integral.val] <- intensities
		}
	}
	## Let's try to make things much easier
	
	if (plottype=="identified")
	{
		## Identifications	
		grouplabels[which(is.na (identifications)==FALSE),1] <- "IDENT"		
		## Contaminants	
		grouplabels[grep("CONTAMINANTDB",identifications),1] <- "CONT"
		if (is.null(DB1string)==FALSE )
		{
		 grouplabels[grep(DB1string,identifications)] <- "DB1"
		}	
		## Related peaks	
		IND <- which(is.na (identifications)==FALSE)	
		for (i in 1:length(IND))
		{
			##assign peak group number
			peakgroupnumber <- IND[i]
			
			## does this peak group has a leationship marker ar basepek?
			relation <- relationships[peakgroupnumber]
			
			## relation ID for current peak			
			RID <- relationids[peakgroupnumber]

			## which peaks have this relation ID 
			
			rgroups <- which(relationids==RID)

			grouplabels[rgroups,2] <- rep(peakgroupnumber,length(rgroups))		
	
			if (length(rgroups)>1)
			{
				rgroups <- rgroups[-c(rgroups%in%peakgroupnumber)]
				## as in some cases related peaks are identified to be raleted to sveral peaks, only first relation will be plotted	
				for (a in 1:length(rgroups))
				{
					if (is.na(grouplabels[rgroups[a],1]))	
					{
						grouplabels[rgroups[a],1] <- paste(grouplabels[peakgroupnumber,1],"_re",sep="")					
					}			
				}		
			}
		}		
		## Unidenfied peaks
		IND <- which(is.na(grouplabels[,1]))
		grouplabels[IND[relationships[IND]=="bp" | relationships[IND]=="potential bp"],1] <- "UNIDENT"
		grouplabels[is.na(grouplabels[,1]),1] <- "UNIDENT_re"
	}

	## Only for dilution serries data. Features which are negatively correlted with decrease of concentration are colored.plottype="correlation"
	if (plottype=="correlation")
	{
		CORCOLORS <- rep(NA,length(goodgroups))
		for (i in 1: length(goodgroups))
		{
			peaktable2 <- peaktable[groupnums==i,]
			COR <- cor (1:nrow(peaktable2),peaktable2[,integral.val])			
			if (COR <= -cor.threshold)
			{
				CORCOLORS[i] <- RED
			}
			if (COR >= cor.threshold)
			{
				CORCOLORS[i] <- GREEN
			}
			if (abs(COR)<cor.threshold)
			{
				CORCOLORS[i] <- BLACK
			}
		}
	}

	## Code for comparison of two samples. plottype="compare"
	if (plottype=="compare")
	{
		COMPARECOLORS <- rep(NA,length(goodgroups))
		sample1 <- which(rownames(xcmsset@phenoData)%in%samples)
		sample2 <- which(rownames(xcmsset@phenoData)%in%samples2)
		sample1groups <- groupnums[which(samplenumber==sample1)]
		sample2groups <- groupnums[which(samplenumber==sample2)]
		BOUTH <- sample1groups[sample1groups%in%sample2groups]		
		COMPARECOLORS[BOUTH] <- BLACK
		COMPARECOLORS[sample1groups[sample1groups%in%BOUTH==FALSE]] <- RED
		COMPARECOLORS[sample2groups[sample2groups%in%BOUTH==FALSE]] <- GREEN
		peaktable2 <- peaktable[samplenumber==sample1,]
		peaktable2 <- rbind (peaktable2,peaktable[samplenumber==sample2,])

	}

	## Do ploting for each sample which was called
	## finaly ploting 
	palette (colors())
	pdf (file=outputfile)
	
	for (i in 1:length(sampnumvector))	
	{	
		peakrows <- which (samplenumber==sampnumvector[i])
		peaktable2 <- peaktable[peakrows,]
		groupnumber <- groupnums[peakrows]

		## defining colors
		if (plottype=="identified")
		{		
			identifieddb1 <- goodgroups[which(grouplabels[,1]=="DB1")]
			identifieddb2 <- goodgroups[which(grouplabels[,1]=="IDENT")]
			identifiedcont <- goodgroups[which(grouplabels[,1]=="CONT")]
			relateddb1 <- goodgroups[which(grouplabels[,1]=="DB1_re")]
			relateddb2 <- goodgroups[which(grouplabels[,1]=="IDENT_re")]
			relatedcont <- goodgroups[which(grouplabels[,1]=="CONT_re")]
			unidentifiedbp <- goodgroups[which(grouplabels[,1]=="UNIDENT")]
			unidentifiedre <- goodgroups[which(grouplabels[,1]=="UNIDENT_re")]
			
			groupcolors <- rep(NA,length(groupnumber))
			groupcolors[groupnumber%in%identifieddb1] <- RED
			groupcolors[groupnumber%in%identifieddb2] <- GREEN
			groupcolors[groupnumber%in%relateddb1] <- PINK
			groupcolors[groupnumber%in%relateddb2] <- LIGHTGREEN
			groupcolors[groupnumber%in%unidentifiedbp] <- BLACK
			groupcolors[groupnumber%in%unidentifiedre] <- GRAY
			groupcolors[groupnumber%in%identifiedcont] <- BLUE
			groupcolors[groupnumber%in%relatedcont] <- LIGHTBLUE

			RTs <- c(peaktable2[groupcolors==GRAY,4],peaktable2[groupcolors==BLACK,4],peaktable2[groupcolors==LIGHTBLUE,4],peaktable2[groupcolors==PINK,4],
			peaktable2[groupcolors==LIGHTGREEN,4],peaktable2[groupcolors==GREEN,4],peaktable2[groupcolors==BLUE,4],peaktable2[groupcolors==RED,4])
			RTs <- RTs/60
			MZs <- c(peaktable2[groupcolors==GRAY,1],peaktable2[groupcolors==BLACK,1],peaktable2[groupcolors==LIGHTBLUE,4],peaktable2[groupcolors==PINK,1],
			peaktable2[groupcolors==LIGHTGREEN,1],peaktable2[groupcolors==GREEN,1],peaktable2[groupcolors==BLUE,4],peaktable2[groupcolors==RED,1])
			size <- peaktable2[,integral.val]
			size <- c(size[groupcolors==GRAY],size[groupcolors==BLACK],size[groupcolors==LIGHTBLUE],size[groupcolors==PINK],
			size[groupcolors==LIGHTGREEN],size[groupcolors==GREEN],size[groupcolors==BLUE],size[groupcolors==RED])
			groupcolors <- c(groupcolors[groupcolors==GRAY],groupcolors[groupcolors==BLACK],groupcolors[groupcolors==LIGHTBLUE],groupcolors[groupcolors==PINK],
			groupcolors[groupcolors==LIGHTGREEN],groupcolors[groupcolors==GREEN],groupcolors[groupcolors==BLUE],groupcolors[groupcolors==RED])
		}
		if (plottype=="correlation")
		{
			groupcolors <- CORCOLORS[groupnumber]
			## Sort RT, m/z and size vectors, that BLACK dots are plotted first, then RED and GREEN are last.
			RTs <- c(peaktable2[groupcolors==BLACK,4],peaktable2[groupcolors==RED,4],peaktable2[groupcolors==GREEN,4])
			RTs <- RTs/60
			MZs <- c(peaktable2[groupcolors==BLACK,1],peaktable2[groupcolors==RED,1],peaktable2[groupcolors==GREEN,1])	
			size <- log2(peaktable2[,9]+0.001)
			size <- c(size[groupcolors==BLACK],size[groupcolors==RED],size[groupcolors==GREEN])	
			groupcolors <- c(groupcolors[groupcolors==BLACK],groupcolors[groupcolors==RED],groupcolors[groupcolors==GREEN])				
		}
		
		if (plottype=="compare")
		{
			groupcolors <- CORCOLORS[groupnumber]
		}
	
		if (bubbletype=="filled")		
		{
			plot (RTs,MZs,pch=PCH,bg=groupcolors,col=substr(groupcolors,1,7),cex=size,lwd=LWD,	
			xlab="Retention time, min.", ylab="m/z",xlim=rtrange,ylim=mzrange,axes=F)
			axis (1, at=seq(0,40,0.5))
			axis (2,at=seq(0,2000,25))
			box ()
			title (main=rownames(xcmsset@phenoData)[sampnumvector[i]])
			title (sub=paste("Int.thr.(min)=",round(intensity.threshold.min,0),", Int.thr (max)=",round(intensity.threshold.max,0),", RT range=",paste(round(rtrange,1),collapse=":"),"min, m/z range=",paste(round(mzrange,1),collapse=":"),", DB1=",DB1string,", scaling=",scale,", logbase=",logbase,sep=""),cex.sub=0.6)
		} else
		{
			plot (RTs,MZs,pch=PCH,col=substr(groupcolors,1,7),cex=size,lwd=LWD,	
			xlab="Retention time", ylab="m/z",xlim=rtrange,ylim=mzrange)
			title (main=rownames(xcmsset@phenoData)[sampnumvector[i]])
			title (sub=paste("Int.thr.(min)=",round(intensity.threshold.min,0),"Int.thr.(max)=",round(intensity.threshold.max,0),", RT range=",paste(round(rtrange,1),collapse=":"),"min, m/z range=",paste(round(mzrange,1),collapse=":"),", DB1=",DB1string,", scaling=",scale,", logbase=",logbase,sep=""),cex.sub=0.6)
		}

		if (smalllegend==TRUE)
		{		
			legend ("topright",pch=PCH,col=substr(c(GREEN,LIGHTGREEN,BLACK,GRAY),1,7),legend=c("Identified","Related","Unidentified","Related(unident.)"),pt.lwd=1,pt.bg=c(GREEN,LIGHTGREEN,BLACK,GRAY,BLUE,LIGHTBLUE),cex=0.8)
		}
		if (is.null(identifieddb1) & plottype=="identified" & smalllegend==FALSE | is.null(DB1string) & plottype=="identified" & smalllegend==FALSE)
		{		
			legend ("topright",pch=PCH,col=substr(c(GREEN,LIGHTGREEN,BLACK,GRAY,BLUE,LIGHTBLUE),1,7),legend=c("Identified","Related","Unidentified","Related(unident.)","Contaminants","Related(cont.)"),pt.lwd=1,pt.bg=c(GREEN,LIGHTGREEN,BLACK,GRAY,BLUE,LIGHTBLUE),cex=0.8)
		}
		if (!is.null(identifieddb1) & plottype=="identified" & !is.null(DB1string) & smalllegend==FALSE)
		{		
			legend ("topright",pch=PCH,col=substr(c(RED,PINK,GREEN,LIGHTGREEN,BLACK,GRAY,BLUE,LIGHTBLUE),1,7),pt.bg=c(RED,PINK,GREEN,LIGHTGREEN,BLACK,GRAY,BLUE,LIGHTBLUE),legend=c("Identified(DB1)","Related(DB1)","Identified","Related","Unidentified","Related(unident.)","Contaminants","Related(cont.)"),pt.lwd=1,cex=0.8)
		}
		if (plottype=="correlation")
		{		
			legend ("topright",pch=PCH,col=substr(c(RED,GREEN,BLACK),1,7),pt.bg=c(RED,GREEN,BLACK),legend=c("Neg.cor.","Pos.cor",paste("mod[cor.]<",cor.threshold,sep="")),pt.lwd=1)
		}
	}
	dev.off ()
}

