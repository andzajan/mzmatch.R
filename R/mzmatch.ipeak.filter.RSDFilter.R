mzmatch.ipeak.filter.RSDFilter <- function(JHeapSize=1425,i=NULL, o=NULL, rejected=NULL, rsd=NULL, h=NULL, v=NULL, sampleList=NULL, nSlaves=1, inputfolder="combined")
{
	version.1 <- get("version.1",envir=.GlobalEnv)

	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")

	if (version.1==TRUE)
	{
		mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	
		## setup the tool
		tool <- paste(mzmatch, "mzmatch.ipeak.filter.RSDFilter")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
		if (!is.null(o))
		tool <- paste(tool, "-o", o)
		if (!is.null(rejected))
		tool <- paste(tool, "-rejected", rejected)
		if (!is.null(rsd))
		tool <- paste(tool, "-rsd", rsd)
		if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
		if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

		system(tool)
	} else
	{
		MainClasses <- levels(sampleList$sampleClass)
		outputfiltered <- paste(inputfolder,"_RSD_filtered", sep="")
		outputrejected<- paste(inputfolder,"_RSD_rejected", sep="")
		if (!file.exists(outputfiltered))
		{
			dir.create (outputfiltered)
		}
		if (!file.exists(outputrejected))
		{
			dir.create (outputrejected)
		}

		RSDfunction <- function (fnum)
		{
			FILESf <- sampleList$outputfilenames[which(sampleList$sampleClass==MainClasses[fnum])]
			if (length(FILESf)>1)
			{
				i <- paste(inputfolder,"/",MainClasses[fnum],".peakml",sep="")
				o <- paste(outputfiltered,"/",MainClasses[fnum],".peakml",sep="")
				rejected <- paste(outputrejected,"/",MainClasses[fnum],".peakml",sep="")

				mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch_2.0.jar", sep="")
	
				## setup the tool
				tool <- paste(mzmatch, "mzmatch.ipeak.filter.RSDFilter")
				if (!is.null(i))
				tool <- paste(tool, "-i", i)
				if (!is.null(o))
				tool <- paste(tool, "-o", o)
				if (!is.null(rejected))
				tool <- paste(tool, "-rejected", rejected)
				if (!is.null(rsd))
				tool <- paste(tool, "-rsd", rsd)
				if (!is.null(h) && h==T)
				tool <- paste(tool, "-h")
				if (!is.null(v) && v==T)
				tool <- paste(tool, "-v")

				system(tool)
			} else
			{
				i <- paste(inputfolder,"/",MainClasses[fnum],".peakml",sep="")
				o <- paste(outputfiltered,"/",MainClasses[fnum],".peakml",sep="")
				file.copy (i,o)
			}
		}

		if (nSlaves>1)
		{
			cl <- makeCluster (nSlaves)
			envname <- environment()
			clusterExport (cl, list=c("RSDfunction","MainClasses","sampleList","rsd","java","h","v"), envir=envname)
			system.time(clusterApply(cl,1:length(MainClasses),RSDfunction))
			stopCluster(cl)
		} else

		{
			for (fnum in 1:length(MainClasses))
			{
				RSDfunction(fnum)
			}
		}

	}
}
