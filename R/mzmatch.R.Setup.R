mzmatch.R.Setup <- function (projectFolder=NULL, samplelist=NULL, outputfolder="peakml", sampleClassFromFolderNames=FALSE, extractMeasurentDate=FALSE)
{
	if (is.null(projectFolder))
	{
		projectFolder <- tclvalue (tkchooseDirectory())
	}
	setwd (projectFolder)

	## supported MS data formats
	extensions <- c("\\.mzXML$", "\\.mzML$", "\\.mzData$", "\\.cdf$", "\\.CDF$")

	if (is.null(samplelist))
	{		

		mzXMLfiles.fullnames <- NULL
		outputfilenames <- NULL
		for (ext in extensions)
		{
			HIT <- dir(full.names=TRUE,pattern=ext,recursive=TRUE)		
			mzXMLfiles.fullnames <- append(mzXMLfiles.fullnames, HIT)
			outputfilenames <- append(outputfilenames,sub(ext,".peakml", HIT))
		}
		
		for (i in 1:length(mzXMLfiles.fullnames))
		{
			outputfilenames[i] <- PeakML.Methods.extractFileName(outputfilenames[i])
		}
		outputfilenames <- paste(outputfolder,"/",outputfilenames, sep="")
		
		duplicate <- table (outputfilenames)
		duplicate <- duplicate[duplicate>1]
		if (length(duplicate)!=0)
		{	
			cat ("Oh no! Output file names:", names(duplicate)," have duplicates. Please check that input files doesn't have a duplicate file names. \n")
			stop()
		}

		sampleList <- data.frame (filenames=mzXMLfiles.fullnames, sampleClass=rep("",length(mzXMLfiles.fullnames)), globalClass=rep("",length(mzXMLfiles.fullnames)), measurementTimeStamp=rep("", length(mzXMLfiles.fullnames)))

		if (sampleClassFromFolderNames==FALSE)
		{
			fix (sampleList)
			## fix writes edited object to the .GlobalEnv, we have to read it back into fucntion environmnet.
			sampleList <- get ("sampleList", envir=.GlobalEnv)
		} else
		{
			#create sample classes, supports subfolders.
			foldernames <- do.call(rbind,strsplit (mzXMLfiles.fullnames,"/"))
			cnames <- apply(foldernames,2, function(x) (length(unique(x))))
			hits <- which(cnames>1)
			if (length(hits)>0)
			{	
				sampleList$sampleClass <- foldernames[,hits[1]]
			}
			
		}

		if (extractMeasurentDate==TRUE)
		{
			sampleList$measurementTimeStamp <- sapply (mzXMLfiles.fullnames,PeakML.Methods.mzML.startTimeStamp)
		}
		
		sampleList[,1] <- as.character(sampleList[,1])
		
		write.table (sampleList,file="sample_setup.tsv",sep="\t",row.names=FALSE)
	}
	else
	{
		sampleList <- read.table (samplelist, sep="\t", header=TRUE)
		mzXMLfiles.shortnames <- rep(NA,nrow(sampleList))
		sampleList[,1] <- as.character(sampleList[,1])
		for (i in 1:nrow(sampleList))
		{
			mzXMLfiles.shortnames[i] <- PeakML.Methods.extractFileName(sampleList$filenames[i])
		}
		outputfilenames <- paste(outputfolder,"/", mzXMLfiles.shortnames, sep="")

		for (ext in extensions)
		{
			outputfilenames <- sub(ext, ".peakml", outputfilenames)
		}

		if (extractMeasurentDate==TRUE)
		{
			sampleList$measurementTimeStamp <- sapply (mzXMLfiles.fullnames,PeakML.Methods.mzML.startTimeStamp)
		}
	}
	if (!file.exists(outputfolder))
	{
		dir.create (outputfolder)
	}
	sampleList$outputfilenames <- outputfilenames
	assign("sampleList", sampleList, envir=.GlobalEnv)
}
