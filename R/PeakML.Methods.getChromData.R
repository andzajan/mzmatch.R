PeakML.Methods.getChromData <- function(filename, PeakMLtree, massCorrection)
{	
	# Generates a chromatorgram data list from peakDataMtx
	# PRE:	PeakMLtree: PeakML file XML tree object.
	# 	massCorrection: mass correction coefficient.
	# POST:	Chromatogram data list

	decode64 <- function (i,object,what,size,endian)
	{
		base64decode (object[[i]],what=what,size=size,endian=endian)
	}

	chromBuild.single <- function (i)
	{
		rts <- base64decode(retentiontimes[[i]],"double",size=4)
		ints <- base64decode(intensities[[i]],"double",size=4)
		massess <- base64decode(masses[[i]],"double",size=4)+massCorrection ## We have to correct masses for ionisation mode
		# Scan ID's are stored with Big Endian in peakml file
		scanidss <- base64decode(scanids[[i]],"integer",size=4,endian="swap")
		measurementidss <- base64decode(measurementids[[i]],"integer",size=4,endian="swap")
		chr <- rbind(massess,ints,rts,scanidss,measurementidss)
		chr
	}

	chromBuild <- function (i)
	{
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(i-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(i-1))
		## We have to correct masses for ionisation mode, that why there is a protonCoef defined
		masses <- .jcall(project, returnSig="[D", method="getMasses", as.integer(i-1))+massCorrection
		scanids <- .jcall(project, returnSig="[I", method="getScanIDs", as.integer(i-1))
		chrom <- rbind(masses,intensities,retentiontimes,scanids)
		chrom
	}

	# Check if peakML file contains data of single measurement or several
	
	numchroms <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peaks/peak/peakdata/retentiontimes"),xmlValue)

	if (length(numchroms)==0)
	{

		intensities <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peakdata/intensities"),xmlValue)
		masses <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peakdata/masses"),xmlValue)
		scanids <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peakdata/scanids"),xmlValue)
		measurementids <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peakdata/measurementids"),xmlValue)
		retentiontimes <- sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peakdata/retentiontimes"),xmlValue)
		
		system.time(chromDataList <- lapply (1:length(masses),chromBuild.single))

	#chromDataList
	
	# base64decode function is too slow, switching back to java there.
	} else
	{
	# create a new Project
		project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
		cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
		st <- system.time(.jcall(project, returnSig="V", method="load",filename))
		cat ("Done in:",st[3],"s \n")

		system.time( chromDataList <- lapply(1:length(numchroms),chromBuild))
	}
	chromDataList
}
