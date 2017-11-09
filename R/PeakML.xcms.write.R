PeakML.xcms.write <- function(xset, filename,ionisation="detect")
{
	#### Java calls	
	## public int getNrPeaksets() - get the number of PeakSets (groups) in memory 
	## public int getNrMassChromatograms() - get the number of masstraces in memory
	####

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

	### ionisation: detect - default, from mzXML, negative, positive, unknown

	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rownames(xset@phenoData), xset@filepaths, as.character(xset@phenoData[,1]))
	## Setting up counter for mass chromatograms with length 0	
	zerocount <- 0	
	nonzerocount <- 0
	
	## Detect which peak pickging algorithm were used to generate XCMS object, for matchedFilter peaks table contain columns "intf", for centWave columns name ir "intb".
	intcolname <- colnames(xset@peaks)[8]
	if (intcolname=="intf")
	{
		peakpicker="matchedFilter"
	}

	if (intcolname=="intb")
	{
		peakpicker="centWave"
	}
	
	if (intcolname!="intb" & intcolname!="intf")
	{
		peakpicker="unknown"	
	}

	# Insert peakpicker method name ir peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character(peakpicker))

	# Insert column names of XCMS set peaktable in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peaktablecolnames"),as.character(paste(colnames(xset@peaks),collapse=";")))

	# Insert column names of XCMS set groups table in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("grouptablecolnames"),as.character(paste(colnames(xset@groups),collapse=";")))	

	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (measurementid in 1:nrow(xset@phenoData))
	{
		for (scannum in 1:length(xset@rt$corrected[[measurementid]]))
		{
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(measurementid-1),xset@rt$corrected[[measurementid]][scannum],as.character(ionisation))
		}		
	}
	
	## Appending RT raw values to PeakMl header
	for (measurementid in 1:nrow(xset@phenoData))
	{
		for (scannum in 1:length(xset@rt$corrected[[measurementid]]))
		{
		#public void addScanAnnotation(int measurementid, int scanid, String label, String value)
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measurementid-1),as.integer(scannum-1),as.character("RT_raw"),as.character(xset@rt$raw[[measurementid]][scannum]))
		}		
	}
	

	## Extract only these peaks from peaktable which are grouped in peaksets (after applying groups() command, xset@groupidx list ir created)
	## Which peak indices are not listed groups table or two features are selected in same group for one sample
	groupvals <- groupval (xset,method="maxint",value="index",intensity="maxo")
	ungroupedpeaks <- which(c(1:nrow(xset@peaks))%in%groupvals==FALSE)
	

	## For speading up, loop only trough these indexes in the peaktable which are related to the sample number	
	## Creating empty list of length equal to sample count in peak table
	sampleindexes <- vector ("list",length(levels(factor(xset@peaks[,"sample"]))))	
	
	cat ("Retrieving peak index from peak table. \n")	
	ST <- system.time(for (measurementid in 1:length(sampleindexes))
	{ 
		cat (rownames(xset@phenoData)[measurementid]," ")	
		sampleindex <- which(xset@peaks[,"sample"]==levels(factor(xset@peaks[,"sample"]))[measurementid])
		## Remove indices of ungrouped peaks		
		REM <- which(sampleindex%in%ungroupedpeaks)
		if (length(REM)!=0)
		{
			sampleindex <- sampleindex[-c(REM)]
		}
		sampleindexes[measurementid] <- list(sampleindex)
	})
	cat ("\n")
	
	## Set indices order to match new structure of PeakML file, the order of xmcs group indices will be changed as peaks are added by sample. We need to reorder xset@groupidx list
	
	setindexes <- vector("list",length(xset@groupidx))	
	for (indexnumber in 1:length(setindexes))
	{	
		gidx <- xset@groupidx[[indexnumber]]
		setindexes[[indexnumber]] <- which(unlist(sampleindexes)%in%gidx)
	}
	
	## Get indices of filled peaks and write to peakML file header as comma separated character
	filledpeaks <- which(unlist(sampleindexes)%in%xset@filled)
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("filledpeaks"),as.character(paste(filledpeaks,collapse=";")))
	
	# bollocks - now we need to load all the files again and retrieve the mass chromatogram data ... idiots; who thinks of these things
	cat("retrieving raw data\n")
	for (measurementid in 1:length(xset@filepaths))
	{
		#print (measurementid)		
		cat("-", xset@filepaths[measurementid], "\n")
		
		# load the raw data
		rawdata <- xcmsRaw(xset@filepaths[measurementid])

		# correction for ionisation mode		
		if (ionisation=="detect") 
		{
			ionisation <- as.character(rawdata@polarity[1])
		}
		
		nrpeaks <- sampleindexes[[measurementid]]
		for (peakid in 1:length(nrpeaks))
		{
			# cat (nrpeaks[peakid]," ")		
			# retrieve the current peak and check whether it is part of the current measurement 
			# this assumes that all the peaks are sorted on xset@peaks["sample"]
			peak <- xset@peaks[nrpeaks[peakid],]
						
			# first locate the index of the rtmin/rtmax in the corrected retention time table
			rettime <- as.numeric(xset@rt$corrected[measurementid][[1]])
			indx_start <- which(rettime == peak["rtmin"])[1]
			indx_finis <- which(rettime == peak["rtmax"])[1]
			
			# This is needed for this crappy FillPeaks routine in xcms, as the for reintegrated peaks rt is not matching actual rt we are choosing this one which is most closest.			
			if (is.na(indx_start))
			{
				
				minrt <- min(rettime[rettime > peak["rtmin"]])	
				if (peak["rtmin"]<min(rettime))
				{
					maxrt <- min(rettime)				
				} else
				{				
					maxrt <- max(rettime[rettime < peak["rtmin"]])
				}				
				if ( (minrt-peak["rtmin"]) <= (peak["rtmin"]-maxrt))
				{
					indx_start <- which(rettime==minrt)[1]
				} else
				{
					indx_start <- which(rettime==maxrt)[1]
				}
			}

			if (is.na(indx_finis))
			{
				if (peak["rtmax"]>max(rettime))
				{
					minrt <- max(rettime)				
				} else
				{				
					minrt <- min(rettime[rettime > peak["rtmax"]])	
				}				
				maxrt <- max(rettime[rettime < peak["rtmax"]])
				if ( (minrt-peak["rtmax"]) <= (peak["rtmax"]-maxrt))
				{
					indx_finis <- which(rettime==minrt)
				} else
				{
					indx_finis <- which(rettime==maxrt)
				}
			}

			# now we can locate the retention time as they will be in the raw data
			rt_start <- xset@rt$raw[measurementid][[1]][indx_start]
			rt_finis <- xset@rt$raw[measurementid][[1]][indx_finis]

			# now we can look up the index of the start and finis scan
			##indx_start <- which(rawdata@scantime == rt_start) - 10
			indx_start <- which(rawdata@scantime == rt_start)
			if (indx_start < 1)
				indx_start <- 1
				##indx_finis <- which(rawdata@scantime == rt_finis) + 10
			indx_finis <- which(rawdata@scantime == rt_finis)
			#print (peakid)
			if (indx_finis > length(rawdata@scantime))
				indx_finis <- length(rawdata@scantime)
			length <- indx_finis - indx_start

			######
			#  Extract data form raw data files
			######
				 
			C <- rawMat (rawdata,mzrange=cbind(peak["mzmin"], peak["mzmax"]),scanrange=c(indx_start,indx_finis))
			C <- rbind (C,NULL)			

			## At some cases, two or more identical RT's are extracted, getting rid of them			
			repeatingRTS <- as.numeric(names(which(table(C[,1])>=2)))
			
			## Removing values with repeating RT's
			## Row with largest intensity are selected				
			if (length(repeatingRTS!=0))
			{			
				for (z in 1:length(repeatingRTS))
				{	
					Csub <- which(round(C[,1],5)==round(repeatingRTS[z],5))
					Csub <- Csub[-c(which(C[Csub,3]==max(C[Csub,3]))[1])]
					C <- C[-c(Csub),]				
				}
			}
					
			## C[,1] : uncorrected RT's
			##  which(unlist(xset@rt$corrected[measurementid])%in%C[,1]): 
			##locate indices of scans which to extract from RT corrected data set								
			scans <- which(unlist(xset@rt$corrected[measurementid])%in%C[,1])
			## if RT correction was applied, scans should be extracted from raw RT's			
			if (length(scans)==0)
			{
				scans <- which(unlist(xset@rt$raw[measurementid])%in%C[,1])
			}			
			
			if (length(scans)<3) 
			{
				scans <- c(-1,-1,-1)
				retentiontimes <- c(-1,-1,-1)
				masses <- c(-1,-1,-1)
				intensities <- c(-1,-1,-1)
				zerocount <- zerocount+1			
			} else {
				retentiontimes <-  unlist(xset@rt$corrected[measurementid])[scans]
				masses <- C[,2]
				intensities <- C[,3]
				nonzerocount <- nonzerocount+1
			}
	
			# finally we can store the mass chromatogram in our project
			# subtract 1 from the measurementid to get in line with java
			.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(measurementid-1), as.integer(scans), retentiontimes,masses, intensities, as.character(ionisation))
			# cat(peakid,length(scans),length(masses),"\n")
				
		}
	}

	## public int getNrPeaksets() - get the number of PeakSets (groups) in memory 
	## public int getNrMassChromatograms() 
	
	## Now we can add extra atributes to masschromatogramm which are related to XCMS peak table
	## public void addAnnotation(int index, String label, String value)			
	peaks <- xset@peaks[unlist(sampleindexes),]
	peaks[is.na(peaks)] <- "NA"		
	for (peakrow in 1:nrow(peaks))
	{			
		for (peakcol in 1:ncol(peaks))
		{		
			.jcall(project, returnSig="V",method="addAnnotation",as.integer(peakrow-1),
			as.character(colnames(peaks)[peakcol]),as.character(peaks[peakrow,peakcol]))
		}
		#cat (peakrow,peakid,"\n")
	}	

	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	for (i in 1:length(setindexes))
	{
		# subtract 1 from the indices to get in line with java
		# print (i)
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setindexes[i][[1]]-1))
	}
	
	## Now we can add extra atributes to masschromatogram sets (groups in XCMS)
	## public void addGroupAnnotation(int groupid, String label, String value)		
	for (groupnumber in 1:length(setindexes))
	{			
		for (i in 1:ncol(xset@groups))
		{		
			.jcall(project, returnSig="V",method="addGroupAnnotation",as.integer(groupnumber-1),
			as.character(colnames(xset@groups)[i]),as.character(xset@groups[groupnumber,i]))
		}
	}

	# and finally store the resulting data
	cat(zerocount,"mass chromatograms of length <3 are removed \n")	
	cat (nonzerocount,"mass chromatograms extracted \n")
	.jcall(project, returnSig="V", method="write", filename)
}

