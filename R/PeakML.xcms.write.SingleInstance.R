PeakML.xcms.write.SingleInstance <- function(xset, outputfile, ionisation=ionisation, addscans=addscans, ppm=ppm, writeRejected=writeRejected, ApodisationFilter=ApodisationFilter)
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

	project <- .jnew("peakml/util/rjava/ProjectSingleMeasurement", as.character(rownames(xset@phenoData)), as.character(xset@filepaths))
	## Setting up counter for mass chromatograms with length 0	
	rejected <- NULL
	accepted <- NULL
	
	# load the raw data
	rawdata <- xcmsRaw(xset@filepaths[1])
	#rawdata <- mzR::openMSfile (xset@filepaths[1])

	# correction for ionisation mode		
	if (ionisation=="detect") 
	{
		#ionisation <- as.character(rawdata@polarity[1])
		ionisation <- PeakML.Methods.getESIpolarity(xset@filepaths[1])
	}

	if (is.na(ionisation)) ionisation <- "neutral"

	cat(ionisation,"\n")
	
	#if (ionisation=="positive" | ionisation=="negative")
	#{
	#	if (length(rawdata@polarity)!=0)
	#	{	
	#		rawdata <- split(rawdata,rawdata@polarity)
	#		rawdata <- rawdata[[which(names(rawdata)==ionisation)]]
	#	}	
	#}

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

	# Insert peakpicker method name in peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character(peakpicker))

	# Insert column names of XCMS set peaktable in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peaktablecolnames"),as.character(paste(colnames(xset@peaks),collapse=";")))
	

	# Insert column names of XCMS set groups table in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("grouptablecolnames"),as.character(paste(colnames(xset@groups),collapse=";")))

	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (scannum in 1:length(xset@rt$corrected[[1]]))
	{
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(0),xset@rt$corrected[[1]][scannum],as.character(ionisation))
	}
	
	## Appending RT raw values to PeakMl header
	for (scannum in 1:length(xset@rt$corrected[[1]]))
	{
		#public void addScanAnnotation(int measurementid, int scanid, String label, String value)
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(0),as.integer(scannum-1),as.character("RT_raw"),as.character(xset@rt$raw[[1]][scannum]))
	}
	
	# bollocks - now we need to load all the files again and retrieve the mass chromatogram data ... idiots; who thinks of these things
	cat("retrieving raw data\n")
	
	#print (measurementid)		
	cat("-", xset@filepaths[1], "\n")
	
	## Make a list of masschromatograms for every peak, so we can filter out peaks with crappy peakshape.
	

	chromatograms <- vector("list",nrow(xset@peaks))
	for (peakid in 1:nrow(xset@peaks))
	{
		#cat (peakid," ")
		##cat (nrpeaks[peakid]," ")		
		# retrieve the current peak and check whether it is part of the current measurement 
		# this assumes that all the peaks are sorted on xset@peaks["sample"]
		peak <- xset@peaks[peakid,]
		mz_start <- peak["mzmin"]-(peak["mzmin"]*ppm*10^-6)
		mz_finis <- peak["mzmax"]+(peak["mzmax"]*ppm*10^-6)

		# first locate the index of the rtmin/rtmax in the corrected retention time table
		rettime <- as.numeric(xset@rt$corrected[1][[1]])
		indx_start <- which(rettime == peak["rtmin"])[1]-addscans
		indx_finis <- which(rettime == peak["rtmax"])[1]+addscans
		
		# if RT's are corrected and does not match original measured ones, detect scan index from the closest RT value
		if (is.na(indx_start))
		{
			indx_start <- which(peak["rtmin"] <= rettime)[1]-1
			if (is.na(indx_start))
			{
				indx_start <- 1
			}
		}
		if (is.na(indx_finis))
		{
			indx_finis <- which(peak["rtmax"] <= rettime)[1]
			if (is.na(indx_finis))
			{
				indx_finis <- length(rettime)
			}
		}

		if (indx_start < 1)
			indx_start <- 1
		if (indx_finis > length(rawdata@scantime))
			indx_finis <- length(rawdata@scantime)

		######
		#  Extract data form raw data files
		######
		C <- rawMat (rawdata,mzrange=cbind(mz_start, mz_finis),scanrange=cbind(indx_start,indx_finis))
		C <- rbind(C,NULL)
		
	
		## At some cases, two or more identical RT's are extracted, getting rid of them			
		## Removing values with repeating RT's
		## Row with largest intensity are selected

		repeatingRTS <- as.numeric(names(which(table(C[,1])>=2)))
		if (length(repeatingRTS!=0)){
			for (z in 1:length(repeatingRTS)){
				Csub <- which(round(C[,1],5)==round(repeatingRTS[z],5))
				Csub <- Csub[-c(which(C[Csub,3]==max(C[Csub,3]))[1])]
				C <- C[-c(Csub),]
				C <- rbind(C,NULL)
			}
		}
		

		##plot (C[,3],type="l",axes=F)
		
		scans <- which(xset@rt$raw[[1]]%in%C[,1])

		if (length(scans)<3 | length(unique(C[,3]))<=3)
		{
			scans <- c(-1,-1,-1)
			retentiontimes <- c(-1,-1,-1)
			masses <- c(-1,-1,-1)
			intensities <- c(-1,-1,-1)
			rejected <- append(rejected,peakid)
		} else {
			retentiontimes <- rettime[scans]
			intensities <- C[,3]
			masses <- C[,2]
			accepted <- append(accepted,peakid)
		}
		# In java scan index should start from 0, so I deduct 1 from all scans
		chromatograms[[peakid]] <- rbind(scans-1,retentiontimes,masses,intensities)
	}

	## Now we can filter on mass chromatogram as we want and how much we want. 
	## Filter out these chromatograms for which intensities at the end or beginnings is max of all intensity. 
	## So removing half detected peaks
	for (chromnum in 1:length(accepted))
	{
		chrom <- chromatograms[[accepted[chromnum]]]
		if (chrom[4,1] >= max(chrom[4,]) | chrom[4,ncol(chrom)] >= max(chrom[4,]))
		{
			rejected <- append (rejected,accepted[[chromnum]])
		}
	}

	if (!is.null(rejected))
	accepted <- accepted[-c(which(accepted%in%rejected))]

	## Filter out centroiding artefacts. It checks a mass withing 0.9 mass units and RT shift of 5s. Remove peaks with intensity less than 2% of maximum intensity.

	
	if (ApodisationFilter==TRUE)
	{
		## mass, intensities for all peaks in accepted list
		## Columns: chromID, RT, mass, intensity
		intensitiesMatrix <- matrix(ncol=4,nrow=length(accepted))
		for (chromnum in 1:length(accepted))
		{
			chrom <- chromatograms[[accepted[chromnum]]]
			intensitiesMatrix[chromnum,]<- c(accepted[chromnum],chrom[2,which(chrom[4,]==max(chrom[4,]))[1]],weighted.mean(chrom[3,],chrom[4,]),max(chrom[4,]))
		}

		## Reorder matrix by mass
		intensitiesMatrix <- intensitiesMatrix[order(intensitiesMatrix[,3]),]

		## Split intensitiesMatrix in the list, by mass window of 0.9 Da
		StartRow <- 1
		massWinList <- vector("list")
		listnum=1
		while (!StartRow>nrow(intensitiesMatrix))
		{
			initMass <- intensitiesMatrix[StartRow,3]
			finMass <- initMass + 0.9
			hits <- which(intensitiesMatrix[,3]>=initMass & intensitiesMatrix[,3]<=finMass)
			massWinList[[listnum]] <- intensitiesMatrix[hits,]
			listnum <- listnum +1
			StartRow <- max(hits)+1
		}

		## Now within each list where masses were stored, order these masses by RT, and split new list, based on RT window.

		RTmassWinList <- vector ("list")
		listnum=1
		for (listindex in 1:length(massWinList))
		{
			datMat <- massWinList[[listindex]]
			datMat <- rbind(datMat,NULL)
			if (nrow(datMat)==1)
			{
				RTmassWinList[[listnum]] <- datMat
				listnum=listnum+1
			} else
			{
				datMat <- datMat[order(datMat[,2]),]
				StartRow <- 1
				while (!StartRow>nrow(datMat))
				{
					initRT <- datMat[StartRow,2]
					finRT <- initRT + 5
					hits <- which(datMat[,2]>=initRT & datMat[,2]<=finRT)
					RTmassWinList[[listnum]] <- datMat[hits,]
					listnum <- listnum+1
					StartRow <- max(hits)+1
				}
			}
		}

		## Now check each list in object RTmassWinList which has more than 1 row, and remove these peaks with intensity less than 2% of max int

		REM <- NULL
		for (listindex in 1:length(RTmassWinList))
		{
			datMat <- RTmassWinList[[listindex]]
			datMat <- rbind(datMat,NULL)
			if (nrow(datMat)!=1)
			{
				maxInt <- max(datMat[,4])
				# Intensity threshold 2%
				intTresh <- maxInt*0.02
				hits <- which(datMat[,4]<intTresh)
				if (length(hits)!=0)
				{
					REM <- append(REM,datMat[hits,1])
				}
			}
		}

		REM <- sort(unique(REM))
		if (length(REM)!=0)
		{
			rejected <- sort(append(rejected,REM))
			accepted <- accepted[-c(which(accepted%in%REM))]
		}
	}

	
	# finally we can store the mass chromatogram in our project
		# subtract 1 from the measurementid to get in line with java
	TESting <- matrix(ncol=4,nrow=length(accepted))
	for (chromnum in 1:length(accepted))
	{
		chrom <- chromatograms[[accepted[chromnum]]]
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(0), as.integer(chrom[1,]), chrom[2,],chrom[3,], chrom[4,], as.character(ionisation))
		TESting[chromnum,]<- c(accepted[chromnum], ncol(chrom),mean(chrom[3,]),max(chrom[4,]))
	}
	
	# and finally store the resulting data
	.jcall (project,returnSig="V",method="writeMeasurements",outputfile)
	
	## Write rejected chromatograms in separate file
	if (writeRejected==TRUE)
	{
		project <- .jnew("peakml/util/rjava/ProjectSingleMeasurement", as.character(rownames(xset@phenoData)), as.character(xset@filepaths))
		for (scannum in 1:length(xset@rt$corrected[[1]]))
		{
			.jcall(project, returnSig="V", method="addScanInfo", as.integer(0),xset@rt$corrected[[1]][scannum],as.character(ionisation))
		}
		for (chromnum in 1:length(rejected))
		{
			chrom <- chromatograms[[rejected[chromnum]]]
			.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(0), as.integer(chrom[1,]), chrom[2,],chrom[3,], chrom[4,], as.character(ionisation))
		}
		for (i in 1:length(rejected))
		{
			.jcall(project, returnSig="V", method="addPeakSet", as.integer(i-1))
		}
		filename2 <- sub (".peakml","",outputfile)
		filename2 <- paste(filename2,"_rejected.peakml",sep="")
		.jcall(project, returnSig="V", method="write", filename2)
	}
	cat (length(accepted)," peaks exported.\n")
	#mzR::close(rawdata)
}

