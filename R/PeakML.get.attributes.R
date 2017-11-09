PeakML.get.attributes <- function(filename,attribute="MeasurementNames",annotations=NULL)
{
	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")
		
	if (attribute=="MeasurementNames")
	{
		attr <- .jcall(project,returnSig="[S", method="getMeasurementNames")
	}
	
	if (attribute=="GroupNames")
	{
		setsnumber <- .jcall(project,returnSig="I",method="getNrPeaksets")
		attr <- rep(NA,setsnumber)
		for (setnum in 1:setsnumber)
		{	
			plotname <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(setnum-1),"identification")
			if (is.null(plotname))
			{
				plotname <- NA
			}
		attr[setnum] <- paste(plotname,sep=",")
		}
	}
	
	if (attribute=="getGroupAnnotation")
	{
		if(!is.null(annotations))
		{
			setsnumber <- .jcall(project,returnSig="I",method="getNrPeaksets")
			attr <- vector("list",length(annotations))
			for (ann in 1:length(annotations))
			{
				for (setnum in 1: setsnumber)
				{
					value <-.jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(setnum-1),as.character(annotations)[ann])
					if (is.null(value))
					{
						value <- NA
					}
					attr[[ann]][setnum] <- value
				}	
			}	
		} else
		{
			cat ("Provide annotation names with parameter 'annotations'. \n")
		}	
	}
	
	if (attribute=="getAnnotation")
	{
		if(!is.null(annotations))
		{
			setsnumber <- .jcall(project,returnSig="I",method="getNrPeaksets")
			attr <- vector("list",length(annotations))
			for (ann in 1:length(annotations))
			{
				for (setnum in 1: setsnumber)
				{
					value <-.jcall(project,returnSig="S", method="getAnnotation",as.integer(setnum-1),as.character(annotations)[ann])
					if (is.null(value))
					{
						value <- NA
					}
					attr[[ann]][setnum] <- value
				}
			}
		} else
		{
			cat ("Provide annotation names with parameter 'annotations'. \n")
		}
	}
	
	if (attribute=="IntensitiesTable")
	{
		cat ("Loading masschromatograms from file \n")
		st <- system.time(chromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
		cat ("Done in:",st[3],"s \n")
		## Table of intensities, sample ID, group ID, AvgMass,AvgRT
		peakdata <- NULL
		for (chromnums in 1:length(chromatograms))
		{
			VEC <- .jevalArray(chromatograms[[chromnums]])[c(9,11,13,8,3)]
			peakdata <- rbind (peakdata,VEC)
		}
		## How many groups (sets are there)
		samplenames <- .jcall(project,returnSig="[S", method="getMeasurementNames")
		## empty output data matrix
		attr <- matrix(ncol=max(peakdata[,3])+1,nrow=length(samplenames))
		for (chromnum in 1:nrow(peakdata))
		{
			intensity <- peakdata[chromnum,1]
			rownum <- peakdata[chromnum,2]+1
			colnum <- peakdata[chromnum,3]+1
			attr[rownum,colnum] <- intensity
		}
		rownames(attr) <- samplenames
		## Calculate average masses for groups
		colnames <- rep(NA,max(peakdata[,3])+1)
		for (groupid in 1:(max(peakdata[,3]+1)))
		{
			masses <- peakdata[peakdata[,3]==sort(unique(peakdata[,3]))[groupid],4]
			colnames[groupid] <- mean(masses,na.rm=TRUE)
		}
		colnames(attr) <- colnames
	}
	if (attribute=="CompleteTable")
	{
		cat ("Loading masschromatograms from file \n")
		st <- system.time(chromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
		cat ("Done in:",st[3],"s \n")
		## Table of: intensities, sample ID, group ID, AvgMass,AvgRT
		peakdata <- NULL	
		for (chromnums in 1:length(chromatograms))
		{
			VEC <- .jevalArray(chromatograms[[chromnums]])[c(9,11,13,8,3)]
			peakdata <- rbind (peakdata,VEC)
		}
		
	
		## Calculate RT at the maximum intensity of the peak to avoid shifted RT's for peaks with long tails
		peakRTs <- rep(NA,length(chromatograms))		
		for (chromnums in 1:length(chromatograms))
		{
			retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(chromnums-1))
			intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(chromnums-1))	
			peakRTs[chromnums] <- retentiontimes[which(intensities==max(intensities))[1]]
		} 

		samplenames <- .jcall(project,returnSig="[S", method="getMeasurementNames")
		## empty output data matrix
		attr <- vector("list",3)
		for (i in 1:3)
		{
			attr[[i]] <- matrix(ncol=max(peakdata[,3])+1,nrow=length(samplenames))
		}
		for (chromnum in 1:nrow(peakdata))
		{
			intensity <- peakdata[chromnum,1]
			mass <- peakdata[chromnum,4]
			RT <- peakRTs[chromnum]
			rownum <- peakdata[chromnum,2]+1
			colnum <- peakdata[chromnum,3]+1
			attr[[1]][rownum,colnum] <- intensity
			attr[[2]][rownum,colnum] <- mass
			attr[[3]][rownum,colnum] <- RT
		}
		for (i in 1:3)
		{
			rownames(attr[[i]]) <- samplenames
		}
		## Calculate average masses for groups
		colnames <- rep(NA,max(peakdata[,3])+1)
		for (groupid in 1:(max(peakdata[,3]+1)))
		{
			masses <- peakdata[peakdata[,3]==sort(unique(peakdata[,3]))[groupid],4]
			colnames[groupid] <- mean(masses,na.rm=TRUE)
		}
		for (i in 1:3)
		{
			colnames(attr[[i]]) <- colnames
		}
	}
	attr
}
