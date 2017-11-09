PeakML.GapFiller <- function(filename,ionisation="detect",Rawpath=NULL,outputfile,ppm=0,rtwin=0, nSlaves=1, fillAll=FALSE)
{
	require (rJava)
	version.1 <- get("version.1",envir=.GlobalEnv)
	#####
	# Methods
	#####

	# .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(0))
	# .jcall(project, returnSig="[D", method="getMasses", as.integer(0))
	# .jcall(project, returnSig="[D", method="getIntensities", as.integer(0))
	# .jcall(project, returnSig="[D", method="getScanIDs", as.integer(0))
	# .jcall(project, returnSig="I", method="getMeasurementID", 0) - for masschromatogram, by giving index of peak.
	# .jcall(project,returnSig="[[D", method="getMassChromatograms")	
	# .jcall(project,returnSig="[S", method="getSetNames")
	# .jcall(project,returnSig="[S", method="getMeasurementNames")
	# .jcall(project,returnSig="[S", method="getFileNames")	
	# .jcall(project,returnSig="S", method="getAnnotation",as.integer(0),"string")
	# .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(0),"string")
	# 				string = relation.id or identification
	# .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes","MeasurementName")
	# .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(MeasurementID))
	
	
	#####
	# values for GetMassChromatograms
	#####
	
	#final int MINRT   =  0;		1
	#final int MAXRT   =  1;		2
	#final int AVGRT   =  2;		3
	#final int MINSCAN  =  3;		4
	#final int MAXSCAN  =  4;		5
	#final int MINMZ   =  5;		6
	#final int MAXMZ   =  6;		7
	#final int AVGMZ   =  7;		8
	#final int INTENSITY  =  8;		9
	#final int SUMINTENSITY =  9; 	10
	#final int MEASUREMENTID = 10; 	11
	#final int SETID   = 11; 		12
	#final int GROUPID  = 12; 		13


	## Global header annotations
	## public void addHeaderAnnotation(String label, String value)	
	## public String getHeaderAnnotation(String label)
	
	## Create scan info, for each ot the measurements
	## public void addScanInfo(int measurementid, double retentiontime, String polarity)
	
	## Add extra info to scans, for example RAW RTs
	## public void addScanAnnotation(int measurementid, int scanid, String label, String value)
	## public String getScanAnnotation(int measurementid, int scanid, String label)

	## Groups or sets		
	## public String getGroupAnnotation(int groupid, String name)
	## public void addGroupAnnotation(int groupid, String label, String value)

	## Masschromatogramms level
	## public String getAnnotation(int index, String name)
	## public void addAnnotation(int index, String label, String value)
	
	## public int getNrScans(int measurementid)
	## public String getIonisation(int measurementid, int scanid)
	
	##.jcall (project,return="D",method="formulaToMass",as.character("[M1];[C4H4]n[+]"))


	## Functions
	FillinPeaks <- function(peaknum){
		whichpeakset <- numchromsexpected[fillinnums[peaknum],1]
		## get a peakset RT and mass window from PeakMLdata$peakDataMtx for all detected samples
		subtable <- PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[,10]==whichpeakset,]
		subtable <- rbind(subtable,NULL)
		#subtable <- apply(subtable,2,mean,na.rm=TRUE)

		## Extract chromatogram
		# now we can locate the retention time as they will be in the raw data
		rt_start <- min(subtable[,5])-rtwin
		rt_finis <- max(subtable[,6])+rtwin
		if (rt_finis > max(correctedRT)){
			rt_finis <- max(correctedRT)
		}
		if (rt_start > max(correctedRT)){
			rt_start <- max(correctedRT)
		}
		if (rt_finis < min(correctedRT)){
			rt_finis <- min(correctedRT)
		}
		if (rt_start < min(correctedRT)){
			rt_start <- min(correctedRT)
		}
		
		mz_start <- min(subtable[,2])
		mz_finis <- max(subtable[,3])
		mz_start<- mz_start-(mz_start*ppm/10^6)
		mz_finis<- mz_finis+(mz_finis*ppm/10^6)

		# Detect MS scan inteval for RT window, if RT correction was used use raw RT's.
		scan_start <- which(correctedRT >= rt_start)[1]-1
		if (scan_start==0) scan_start=1
		scan_finis <- which(correctedRT >= rt_finis)[1]

		######
		#  Extract data form raw data files
		######
		C <- try(PeakML.Methods.getRawMat (allRawPeaks, scan_start, scan_finis, mz_start, mz_finis,correctedRT,uncorrectedRT))
		if (class(C)=="try-error")
		{
			C <- c(1,1,1,1,1)
		}
		C <- rbind (C,NULL)

		## Artefacts where intensities in all scans is constant will be also removed.
		if (nrow(C)<=3 | length(unique(C[,5]))<=3){
			scanids <- c(-1,-1,-1)
			retentiontimes <- c(-1,-1,-1)
			masses <- c(-1,-1,-1)
			intensities <- c(-1,-1,-1)
#			zerocount <- zerocount+1
		} else {
			scanids <- C[,3]
			retentiontimes <- C[,2]
			masses <- C[,4]
			intensities <- C[,5]
#			nonzerocount <- nonzerocount+1
#			fillednums <- append(fillednums,fillinnums[peaknum])
		}
			#assign ("zerocount",zerocount,envir=.GlobalEnv)
#			assign ("nonzerocount",nonzerocount,envir=.GlobalEnv)
#			assign ("fillednums",fillednums,envir=.GlobalEnv)
			#cat (peaknum," ")
			OUT <- rbind(masses,intensities,retentiontimes,scanids-1)
	}

	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))

	ionisation <- PeakMLdata$massCorrection[[2]]
	massCorrection <- PeakMLdata$massCorrection[[1]]
	samplenames <- PeakMLdata$sampleNames
	rawdatafullpaths <- PeakMLdata$rawDataFullPaths
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	## Vector of all possible peaks for every sample in peakset
	numchromsexpected <- unlist(lapply(1:max(PeakMLdata$peakDataMtx[,10]),function (x) rep(x,length(samplenames))))

	## This table is used to generate information for chromatograms extraction/filling. 1st column - peakset number. 2nd-column presence/absence of the chromatogram in the peakml input file. 3rd column - data file number. 4th columns - index of chromatograms in PeakMLdata$chromDataList

	numchromsexpected <- cbind(numchromsexpected,NA,NA,NA)

	for (setnum in 1:max(PeakMLdata$peakDataMtx[,10])){
		inset <- c(1:length(samplenames))
		rownums <- which(PeakMLdata$peakDataMtx[,10]==setnum)
        hit <- PeakMLdata$peakDataMtx[rownums,9]
        oneEqualsSetnum = which(numchromsexpected[,1]==setnum)
		numchromsexpected[oneEqualsSetnum,2] <- as.numeric(inset%in%hit)
		missed <- which(inset%in%hit==FALSE)
		if (length(missed)>0)
		{
			detectedpeaks <- c(rep(1,length(hit)),rep(0,length(missed)))
			hit <- append(hit,missed)
			rownums <- append(rownums,rep(0,length(missed)))
		} else
		{
			detectedpeaks <- rep(1,length(hit))
		}
		numchromsexpected[oneEqualsSetnum,3] <- hit
		numchromsexpected[oneEqualsSetnum,2] <- detectedpeaks
		numchromsexpected[oneEqualsSetnum,4] <- rownums
	}
	colnames(numchromsexpected) <- NULL

	## List object with mass chromatograms which are already in file
	chromslist <- vector("list",nrow(numchromsexpected))
	
	## Insert chromatograms from peakml file in the list
	#detected <- which(numchromsexpected[,2]==1)
	
	#system.time(
	#for (chrnum in 1:nrow(PeakMLdata$peakDataMtx)){
	#	chromslist[[detected[chrnum]]] <- rbind(masses,intensities,retentiontimes,scanids)
	#})


	#if fillAll is set to TRUE, all peaks will be reintegrated with given RT and mass window.
	if (fillAll==TRUE)
	{
		detectedchromatograms <- numchromsexpected[,2]
		numchromsexpected[,2] <- 0
	}

	## Now fill in missing peaks
	notdetected <- which(numchromsexpected[,2]==0)
	whichfiles <- numchromsexpected[notdetected,3]
	samplenums <- unique(whichfiles)
#	zerocount <- 0
#	nonzerocount <- 0
#	fillednums <- NULL

	if (length(samplenums!=0)){
		for (filenum in 1:length(samplenums)){
			samplefile <- samplenums[filenum]

			# Retention times vector is used to detect scan numbers
			# rawfile <- xcmsRaw(rawdatafullpaths[samplefile])
            cat('Working on file: ', rawdatafullpaths[samplefile], '\n')
			rawfile <- openMSfile (rawdatafullpaths[samplefile],verbose=FALSE)
            allRawPeaks <- peaks(rawfile)
			correctedRT <- as.numeric(PeakMLdata$correctedRTList[[samplefile]])
			uncorrectedRT <- header(rawfile)$retentionTime
			if (all(correctedRT==uncorrectedRT))
			{
				rtCorrection <- FALSE
			} else
			{
				rtCorrection <- TRUE
			}

			## detect which peaks to fill in for current data file
			fillinnums <- notdetected[whichfiles==samplefile]
			
			isSnow <- FALSE
			if (nSlaves>1){
				HIT <- grep("snow",installed.packages()[,1])
				if (!is.null(HIT)) 
				{				
					isSnow==TRUE
				} else
				{
					cat("Pleae install package snow to use multiple processors. \n We will continue with a single processor for the time being.","\\n")
				}

			}
			if (isSnow==TRUE){
				#library (snow)
				if (filenum==1)
				{
					cat("Package snow loaded.","\n")
				}
				cl <- makeCluster (nSlaves)
				assign ("rtwin",rtwin,envir=.GlobalEnv)
				assign ("rawfile",rawfile,envir=.GlobalEnv)
				assign ("numchromsexpected",numchromsexpected,envir=.GlobalEnv)
				assign ("fillinnums",fillinnums,envir=.GlobalEnv)
				assign ("PeakMLdata$peakDataMtx",PeakMLdata$peakDataMtx,envir=.GlobalEnv)
				assign ("ppm",ppm,envir=.GlobalEnv)
				assign ("FillinPeaks",FillinPeaks,ppm,envir=.GlobalEnv)
				assign ("rawMat",rawMat,envir=.GlobalEnv)
				clusterExport(cl, list=c("rtwin","rawfile", "numchromsexpected", "fillinnums", "PeakMLdata$peakDataMtx", "ppm","FillinPeaks","rawMat"))
				filledlist <- parLapply(cl, c(1:length(fillinnums)), FillinPeaks)
				stopCluster (cl)
			}else{
				filledlist <- lapply(1:length(fillinnums),FillinPeaks)
			}
			
			## For debug purpose only
			#for (bd in 1: length(fillinnums))
			#{
			#	cat (bd,"\n")
			#	FillinPeaks(bd)
			#}

#			assign ("zerocount",zerocount,envir=.GlobalEnv)
#			assign ("nonzerocount",nonzerocount,envir=.GlobalEnv)
#			assign ("fillednums",fillednums,envir=.GlobalEnv)
			
			for (i in 1:length(filledlist)){
				chromslist[[fillinnums[i]]] <- filledlist[[i]]
			}
			mzR::close (rawfile)
			rm (filledlist)
		}
	}
	
	##SetNames, if peakml file has a several peaksets, they will be restored.
	#samplegroups=PeakMLdata$peakDataMtx[,11]
	#sampleclasses <- .jcall(project,returnSig="[S", method="getSetNames")

	#samplelookfunction <- function (x)	
	#{
	#	phenoData <- sampleclasses[samplegroups[PeakMLdata$peakDataMtx[,9]==x][1]]
	#	phenoData
	#}
	#phenoData <- PeakMLdata$phenoData

	## Wipe out some unused large objects before writing
	# rm (masschromatograms,project)
	# rm (PeakMLdata) # as masschromatograms is inside a function now
	## Write all of this out
	
	project <- .jnew("peakml/util/rjava/Project", samplenames, rawdatafullpaths, as.character(PeakMLdata$phenoData))
	# Insert peakpicker method name ir peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character("XCMS_Gapfilled"))
	
	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (measurementid in 1:length(samplenames)){
		for (scannum in 1:length(PeakMLdata$correctedRTList[[measurementid]])){
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(measurementid-1),as.numeric(PeakMLdata$correctedRTList[[measurementid]][scannum]),as.character(ionisation))
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measurementid-1),as.integer(scannum-1),as.character("RT_raw"),as.character(PeakMLdata$rawRTList[[measurementid]][scannum]))
		}
	}
	

	# finally we can store the mass chromatogram in our project
	# subtract 1 from the measurementid to get in line with java
	# Filled chromatograms are taken from chromslist object, but existing ones from PeakMLdata$chromDataList

	for (i in 1:length(chromslist)){
		if (numchromsexpected[i,2]==0)
		{
			chrom <- chromslist[[i]]
		} else
		{
			ind <- numchromsexpected[i,4]
			chrom <- PeakMLdata$chromDataList[[ind]]
		}
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(numchromsexpected[i,3]-1), as.integer(chrom[4,]), as.numeric(chrom[3,]),as.numeric(chrom[1,]), as.numeric(chrom[2,]), as.character(ionisation))
	}

	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	## Make a list with peak numbers which should be grouped together in sets
	setindexes <- vector("list",length(unique(numchromsexpected[,1])))
	
	for (indexnumber in 1:length(setindexes)){
		setindexes[[indexnumber]] <- which(numchromsexpected[,1]==indexnumber)
	}
	for (ind in 1:length(setindexes)){
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setindexes[[ind]]-1))
	}


	## Restore peak annotations data

	if (!is.null(PeakMLdata$GroupAnnotations))
	{
		PeakML.Methods.writeGroupAnnotations (project, PeakMLdata$GroupAnnotations)
	}

	## Add Annotations to the peaks which are filled in, and which not.
	
	
	## Now we can add extra atributes to masschromatogram sets (groups in XCMS)
	## public void addGroupAnnotation(int groupid, String label, String value)
#	fillindex <- rep (0,length(setindexes))
#	## Peak sets which were appended
#	fillindex[unique(numchromsexpected[fillednums,1])] <- 1
#	for (groupnumber in 1:length(setindexes)){
#		.jcall(project, returnSig="V",method="addGroupAnnotation",as.integer(groupnumber-1), as.character("gap.filled"),as.character(fillindex[groupnumber]))
#	}

	## Finally write file
	#cat(zerocount," of fillled in mass chromatograms are removed (length <3 scans)\n")
	#cat (nonzerocount,"mass chromatogramms filled in.\n")
	.jcall(project, returnSig="V", method="write", outputfile)
}
