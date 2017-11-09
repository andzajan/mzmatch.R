PeakML.Methods.getPeakData <- function(PeakMLtree, chromDataList)
	{
	# Returns peak data from mass chromatograms generated from PeakML files
	# PRE:  javaProject: The java project
	#	ionisation: "detect" will detect the ionisation mode from PeakML files
	# POST: peakDataMtx: A matrix containing matrix of 11 columns, column names : "AVGMZ","MINMZ","MAXMZ","AVGRT (this is recalculated at maximum 
	#	intensity of peak, warning -this approach differs from that one which is used in XCMS)", "MINRT",  "MAXRT", "SUMINTENSITY",
	#	"MAXINTENSITY", "measurement id", "Peak group ID" , "SET ID".

	RTfunction <- function (i)
	{
		hit <- which(chromDataList[[i]][2,]==max(chromDataList[[i]][2,]))[1]
		chromDataList[[i]][3,hit]
	}

	GroupIDfunction <- function (i)
	{
		STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/peaks/peak[",i,"]/peaks/peak/measurementid\"),xmlValue)",sep="")
		REPS <- length(eval(parse(text=STR)))
		rep(i,REPS)
	}

	cat ("Extracting peak data from PeakMl file,\n")

	SetMeasurementids <- sapply(getNodeSet(PeakMLtree,"/peakml/header/sets/set/measurementids"),xmlValue)
	if (length(SetMeasurementids>0))
	{
		SetMeasurementids <- lapply (1:length(SetMeasurementids),function(i) {base64decode(SetMeasurementids[i], what="integer",endian="swap")+1})
	}
	
	nrPeakSets <- as.numeric(sapply(getNodeSet(PeakMLtree, "/peakml/header/nrpeaks"),xmlValue))

	st <- system.time({
	AVGMZ <- sapply (1:length(chromDataList),function(i) {weighted.mean(chromDataList[[i]][1,],chromDataList[[i]][2,])})
	MINMZ <- sapply (1:length(chromDataList),function(i) {min(chromDataList[[i]][1,])})
	MAXMZ <- sapply (1:length(chromDataList),function(i) {max(chromDataList[[i]][1,])})
	RT <- sapply (1:length(chromDataList),RTfunction)
	MINRT <- sapply (1:length(chromDataList),function(i) {min(chromDataList[[i]][3,])})
	MAXRT <- sapply (1:length(chromDataList),function(i) {max(chromDataList[[i]][3,])})
	SUMINTENSITY <- sapply (1:length(chromDataList),function(i) {sum(chromDataList[[i]][2,])})
	MAXINTENSITY <- sapply (1:length(chromDataList),function(i) {max(chromDataList[[i]][2,])})
	
	## If peakML file contains only 1 sample, there are no GROUPSETINDEX and GROUPID.
	if (length(SetMeasurementids)==0)
	{
		MEASUREMENTID <- 1
		GROUPSETSINDEX <- NA
		GROUPID <- NA
	} else
	{
		MEASUREMENTID <- as.numeric(sapply(getNodeSet(PeakMLtree,"/peakml/peaks/peak/peaks/peak/measurementid"),xmlValue))+1
		GROUPSETSINDEX <- rep(NA,length(chromDataList))
		for (i in 1:length(SetMeasurementids))
		{
			GROUPSETSINDEX[MEASUREMENTID%in%SetMeasurementids[[i]]] <- i
		}
		GROUPID <- unlist(lapply(1:nrPeakSets,GroupIDfunction))
	}
	peakDataMtx <- cbind (AVGMZ,MINMZ,MAXMZ,RT,MINRT,MAXRT,SUMINTENSITY,MAXINTENSITY,MEASUREMENTID,GROUPID,GROUPSETSINDEX)
	rownames(peakDataMtx) <- NULL
	})

	cat ("Peak data created in ",st[3],"s \n")
	peakDataMtx
}

