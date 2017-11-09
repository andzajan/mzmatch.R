PeakML.xcms.read <- function(filename,ionisation="detect",Rawpath=NULL,annotations=FALSE, version.1)
{
	#####
	# Methods
	#####

	# .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(0))
	# .jcall(project, returnSig="[D", method="getMasses", as.integer(0))
	# .jcall(project, returnSig="[D", method="getIntensities", as.integer(0))
	# .jcall(project, returnSig="I", method="getMeasurementID", 0) - for masschromatogram, by giving index of peak.
	# .jcall(project,returnSig="[[D", method="getMassChromatograms")	
	# .jcall(project,returnSig="[S", method="getSetNames")
	# .jcall(project,returnSig="[S", method="getMeasurementNames")
	# .jcall(project,returnSig="[S", method="getFileNames")	
	# .jcall(project,returnSig="S", method="getAnnotation",as.integer(0),"string")
	# .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(0),"string")
	# 								string = relation.id or identification
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

	#Correction factor by ionisation mode
	hydrogen <- 1.00782503214
	electron <- 0.00054857990924
	protonmass <- hydrogen-electron
		
	
	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")
	
	#Detect ionisation mode, only first scan from first sample is used, no ionisation switching is supported yet
	
	
	if (ionisation=="detect")
	{
		ionisation<- .jcall (project,returnSig="S",method="getIonisation",as.integer(0),as.integer(0))		
	}
	
	if (version.1==FALSE)
	{
		protonCoef <- 0
	} 
	if (ionisation=="positive" & version.1==TRUE)
	{
		protonCoef <- 1
	}
	if (ionisation=="negative" & version.1==TRUE)
	{
		protonCoef <- -1
	} 
	
	## Read samplenames from peakml file	
	samplenames <- .jcall(project,returnSig="[S", method="getMeasurementNames")	
	
	# Check if data directory for raw files is readable, and if files are present there, 
	# if RaeDataPath is not set, check for raw data in default folder which is defined in PeakML file	
	if (!is.null(Rawpath))
	{
		directorycontent <- dir(Rawpath,recursive=TRUE)
		matchedfiles <- which(directorycontent%in%samplenames)
	  	if(length(directorycontent[directorycontent%in%samplenames])==length(samplenames))
			{
			fileslocated="Y"	
			rawdatafullpaths <- dir(Rawpath,full.names=TRUE,recursive=TRUE)[matchedfiles]	
			cat(paste ("Raw data file located at: ",rawdatafullpaths,"\n",sep=""),"\n")
			} else
			{
			fileslocated="N"
			errormessage = paste("Raw data file: ",samplenames,"cant't be found in folder: ",Rawpath,"\n",sep="")
			}
	} else
	{
		rawdatafullpaths <- .jcall(project,returnSig="[S", method="getFileNames")
		## Check if files are readable at location defined in peakml file		
		locatedfiles <- which(file.exists(rawdatafullpaths)==TRUE)	
		if (length(locatedfiles)!=length(rawdatafullpaths))	
		{
		fileslocated="N"
		if (length(locatedfiles)!=0)
		{
			rawdatafullpaths <- rawdatafullpaths[-c(locatedfiles)]			
		}	
		errormessage = paste("Raw data file can't be read from location': ",rawdatafullpaths,"\n",sep="")
		} else
		{
			fileslocated="Y"
			cat(paste ("Raw data file located at: ",rawdatafullpaths,"\n",sep=""),"\n")
		}
	}

	if (fileslocated=="N")
	{
	cat (errormessage)
	cat ("Warning: some of the raw data files are not accessible, you will not be able to extract ion chromatograms from XCMS object. \n")
	}

	# Extract mass chromatograms
	cat ("Extracting mass chromatograms from peakml data file")	
	st <- system.time(masschromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
	cat (", done in:",st[3],"s \n")
		
	## Creating epmty xcmsSet objects, and filling in xcms structure
	xset <- new ("xcmsSet")
	
	## profinfo, by default these arbitary values are inserted for centWave, should check it
	xset@profinfo <- data.frame(method="bin",step=0.1,stringsAsFactors=FALSE)

	## Check if PeakML file ir generated by xcms to PeakML(PeakML.xcms.write function) writer.
	## Extra attributes: peakproc,peaktablecolnames, grouptablecolnames and filledpeaks are defined in PeakML file header 
	## public String getHeaderAnnotation(String label)
	peakproc <- .jcall (project,returnSig="S",method="getHeaderAnnotation",as.character("peakproc"))
	if (!is.null(peakproc))
	{
		peaktablecolnames <- .jcall(project,returnSig="S",method="getHeaderAnnotation",as.character("peaktablecolnames"))
		grouptablecolnames <- .jcall(project,returnSig="S",method="getHeaderAnnotation",as.character("grouptablecolnames"))
		filledpeaks <- .jcall(project,returnSig="S",method="getHeaderAnnotation",as.character("filledpeaks"))
		peaktablecolnames <- unlist(strsplit (peaktablecolnames,";"))
		grouptablecolnames <- unlist(strsplit (grouptablecolnames,";"))
		filledpeaks <- unlist(strsplit (filledpeaks,";"))
	}		

	## This tables reads generic masschrograms information from PeakML file, which should be present in all PeakML files
	## Output: matrix of 11 columns, column names : "AVGMZ","MINMZ","MAXMZ","AVGRT (this is recalculated at maximum intensity of peak, warning - this approach differs from that one which is used in XCMS)", "MINRT","MAXRT","SUMINTENSITY","MAXINTENSITY","measurement id","group id","peakcount"
	reconstructpeaktable <- function (i)
	{
		peakmldata <- .jevalArray(masschromatograms[[i]])
		peakdata <- matrix(ncol=11,nrow=1)		
		# AVGMZ		
		peakdata[1] <-  peakmldata[8]+(protonmass*protonCoef)
		# MINMZ		
		peakdata[2] <-  peakmldata[6]+(protonmass*protonCoef)
		# MAXMZ
		peakdata[3] <-  peakmldata[7]+(protonmass*protonCoef)
		# AVGRT				
		#peakdata[4] <-  peakmldata[3]
		## Calculate RT at the maximum intensity of the peak to avoid shifted RT's for peaks with long tails
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(i-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(i-1))	
		peakdata[4] <- retentiontimes[which(intensities==max(intensities))[1]]
		# MINRT		
		peakdata[5] <-  peakmldata[1]
		# MAXRT		
		peakdata[6] <-  peakmldata[2]
		# SUMINTENSITY		
		peakdata[7] <-  peakmldata[10]
		# MAXINTENSITY
		peakdata[8] <-  peakmldata[9]
		# measurement id, +1 because java index starts at 0		
		peakdata[9] <- peakmldata[11]+1
		# group indexes, in next loop converted to list as required by xcms structure		
		peakdata[10] <- peakmldata[13]+1
		## Fort this stupid group table we need to count how many peaks from each sample class are grouped together.
		## So this vector will be used to genarate proper groups table structure
		peakdata[11] <- peakmldata[12]+1	
		peakdata	
	}	
		
	cat ("Extracting peak data from PeakMl file,")
	ST <- system.time(peakdata <- do.call (rbind,lapply(1:length(masschromatograms),reconstructpeaktable)))
	cat ("done in",ST[3],"s \n")
		
	groupindex <- peakdata[,10]
	samplegroups <- peakdata[,11]
		
	## Building up peaks table, if peakml file was exported by XCMS original colnames will be restored, in other case default table with column names that match "matchedFilter" algorithm will be used.
	## public String getAnnotation(int index, String name)
	if (!is.null(peakproc))
	{
		cat ("Reading content of original XCMS peak table, ")		
		getPeakAnnotations <- function (rownumber,annotations)
		{
			VAL <- rep (NA,length(peaktablecolnames))
			for (annotations in 1:length(peaktablecolnames))
			{
				VAL[annotations] <- .jcall (project, returnSig="S",method="getAnnotation",as.integer(rownumber-1),as.character(peaktablecolnames[annotations]))
			}
			VAL
		}				
		ST <- system.time (peakstable <- lapply(1:nrow(peakdata),getPeakAnnotations))
		peakstable <- do.call(rbind,peakstable)	
		colnames(peakstable) <- peaktablecolnames
		cat ("done in",ST[3],"s \n")
	}
	
	## If there is no peakproc set in PeakML file header, make a peaks matrix with dummy values
	if (exists("peakstable"))
	{
		restoredpeakstable <- matrix(nrow=nrow(peakstable),ncol=length(peaktablecolnames)) 
		colnames (restoredpeakstable) <- c(peaktablecolnames)
		restoredpeakstable[,"mz"] <- peakdata[,1]
		restoredpeakstable[,"mzmin"] <- peakdata[,2]
		restoredpeakstable[,"mzmax"] <- peakdata[,3]
		restoredpeakstable[,"rt"] <- peakdata[,4]
		restoredpeakstable[,"rtmin"] <- peakdata[,5]
		restoredpeakstable[,"rtmax"] <- peakdata[,6]
		restoredpeakstable[,"into"] <- peakdata[,7]
		restoredpeakstable[,"maxo"] <- peakdata[,8]
		restoredpeakstable[,"sample"] <- peakdata[,9]
		## Detect wihich columns are not filled yet		
		emptycolumns <- peaktablecolnames[!peaktablecolnames%in%c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","maxo","sample")]
		for (colnum in 1:length(emptycolumns))
		{
			vals <- peakstable[,emptycolumns[colnum]]
			vals[vals=="NA"] <- NA
			restoredpeakstable[,emptycolumns[colnum]] <- as.numeric(vals)
		}
	} else
	{
		cat ("XCMS package was not used to process this data set, generating \"generic\" peak table(mimics find.Peaks.centWave output). \n")
		restoredpeakstable <- matrix(nrow=nrow(peakdata),ncol=12) 
		colnames (restoredpeakstable) <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn","egauss","sample")
		restoredpeakstable[,"mz"] <- peakdata[,1]
		restoredpeakstable[,"mzmin"] <- peakdata[,2]
		restoredpeakstable[,"mzmax"] <- peakdata[,3]
		restoredpeakstable[,"rt"] <- peakdata[,4]
		restoredpeakstable[,"rtmin"] <- peakdata[,5]
		restoredpeakstable[,"rtmax"] <- peakdata[,6]
		restoredpeakstable[,"into"] <- peakdata[,7]
		restoredpeakstable[,"intb"] <- peakdata[,7]
		restoredpeakstable[,"maxo"] <- peakdata[,8]
		restoredpeakstable[,"sample"] <- peakdata[,9]
		restoredpeakstable[,"sn"] <- 1
		restoredpeakstable[,"egauss"] <- 0.1
	}

	xset@peaks <- restoredpeakstable	
	
	##PhenoData	
	sampleclasses <- .jcall(project,returnSig="[S", method="getSetNames")
	
	## extract and sort unique sample numbers from peak table	
	samplecount <- sort(unique (xset@peaks[,"sample"]))
	## for each unique sample number, look for class labels and generate phenoData table
	## First row from peaks table which much sample number is taken
	samplelookfunction <- function (x)	
	{
		phenoData <- data.frame(class=sampleclasses[samplegroups[xset@peaks[,"sample"]==x][1]])
		phenoData
	}
	phenoData <- do.call(rbind,lapply(1:length(samplenames),samplelookfunction)) 
	rownames(phenoData) <- samplenames
	
	## Remove samples which are not present in all peaksets
	#REM <- which(is.na(xset@phenoData[,1]))
	#if (length(REM)!=0)
	#{
	#	phenoData <- data.frame(as.matrix(phenoData)[-c(REM),])
	#	rawdatafullpaths <- rawdatafullpaths[-c(REM)]
	#	colnames(phenoData) <- "class"
	#}
	
	xset@phenoData <- phenoData
	xset@filepaths <- rawdatafullpaths
	
	
	
	
	##Restoring this stupid list of group indexes	
	groupidx <- vector("list",length(unique(groupindex)))
	
	cat ("Restoring group index, ")
	ST <- system.time(for (i in 1: length(unique(groupindex)))
	{
		groupidx[[i]] <- which(unique(groupindex)[i]==groupindex)
	})
		
	xset@groupidx <- groupidx
	cat ("done in",ST[3],"s. \n")	
	
	## Ionisation mode
	xset@polarity <- ionisation
	
	##Filling in groups table ???
	
	groupeaksfunction <- function (i)
	{
		groups <- matrix(ncol=7+length(levels(xset@phenoData[,1])),nrow=1)
		grouppeaks <- xset@peaks[groupidx[[i]],]		
		## avoid R to interpret matrix wiht 1 row as vector
		grouppeaks <- rbind(grouppeaks,NULL)		
		if (nrow(grouppeaks)>1)
		{		
			groups[c(1:6)] <- apply (grouppeaks[,c(1:6)],2,median)
		} else
		{
			groups[1:6] <- grouppeaks[,1:6]			
		}
		groups[7] <- nrow(grouppeaks)
		## fill in for which sample classes peaks ar detected, I am not sure that we need this at all
		detections <- rep(0,length(levels(xset@phenoData[,1])))
		## vector of length=nsamp and filled with zeros		
		isdetectedinsample <- rep(0,nrow(xset@phenoData))
		## If peak is detected for current sample, replace it with "1"	
		isdetectedinsample[c(1:nrow(xset@phenoData))%in%grouppeaks[,"sample"]] <- 1	
		## for each sample class, count amount of samples in which peak is deteced		
		for (a in 1:length(levels(xset@phenoData[,1])))
		{
			detections[a] <- sum(isdetectedinsample[!is.na(xset@phenoData[,1])][xset@phenoData[!is.na(xset@phenoData[,1]),1]==levels(xset@phenoData[,1])[a]])		
		}
		groups[8:(7+length(detections))] <- detections	
		groups
	}	
		
	ST <- system.time(groups <- do.call(rbind,lapply(1:length(groupidx),groupeaksfunction)))	
	groups <- rbind(groups,NULL)
	colnames(groups) <- c("mzmed","mzmin","mzmax","rtmed","rtmin","rtmax","npeaks",as.character(levels(xset@phenoData[,1])))
	xset@groups <- groups
	cat ("Groups table extracted in",ST[3],"s \n")

	## Filling in RT's, no corrected RT's are availiable
	rts <- vector("list",nrow(xset@phenoData))
	rts_raw <- vector("list",nrow(xset@phenoData))
	ST <- system.time(for (measurementid in 1:nrow(xset@phenoData))
	{
		rts[[measurementid]] <- .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(measurementid-1))					
	})	
	cat ("RT's extracted in",ST[3],"s. \n")
	
	## Check if raw RT's are stored, for first sample, first scan 
	rawrt <- as.numeric(.jcall(project,returnSig="S",method="getScanAnnotation",as.integer(0), as.integer(0), as.character("RT_raw")))	
	if (length(rawrt)!=0)
	{
		cat ("Values for uncorrected RT's are found, ")	
		ST <- system.time(for (measurementid in 1:nrow(xset@phenoData))
		{
			rtsvec <- rep(NA,length(rts[[measurementid]]))
			for (scannum in 1:length(rtsvec))
			{
				rtsvec[scannum] <- as.numeric(.jcall(project,returnSig="S",method="getScanAnnotation",as.integer(measurementid-1), as.integer(scannum-1), as.character("RT_raw")))	
			}	
			rts_raw[[measurementid]] <- rtsvec	
		})
		cat ("extracted in",ST[3],"s \n")
	} else
	{
		rts_raw <- rts	
	}

	xset@rt <- list(raw=rts_raw,corrected=rts)		
		
	## Adding peakml extra information (identification, relation.id - related peak indexes, relationship - relation identificator) to XCMS object (this is not used by xcms package), output is a list of xcms object and extra peakml information.		
	if (annotations==TRUE)
	{
		cat ("Extracting peakml extra information, identifications and relation id's")
		identifications <- rep (NA,nrow(xset@groups))
		relationids <- rep (NA,nrow(xset@groups))
		relationships <- rep(NA,nrow(xset@groups))	
		dummies <- rep(NA,nrow(xset@groups))		
		ST <- system.time(for (i in 1:nrow(xset@groups))
		{			
			identification <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(i-1),"identification")
			relationid <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(i-1),"relation.id")
			relationship <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(i-1),"relation.ship")	
			dummy <- .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(i-1),"relation.ship")
			if (is.null(identification)) {identification <- NA}
			if (is.null(relationid)) {relationid <- NA}
			if (is.null(relationship)) {relationship <- NA}
			if (is.null(dummy)) {dummy <- NA}
			identifications[i] <- identification
			relationids[i] <- relationid
			relationships[i] <- relationship
			dummies[i] <- dummy 		
		})
		peakmldata <- cbind (identifications,relationids,relationships)		
		cat (ls (),"\n")				
		cat ("peakml extra information extracted in",ST[3],"s. \n")	
		set <- list("vector",2)
		set[[1]] <- xset
		cat("what","\n")
		set[[2]] <- peakmldata
	} else
	{
		set <- xset
	}
	set
}

