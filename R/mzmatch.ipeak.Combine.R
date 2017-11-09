mzmatch.ipeak.Combine <- function(JHeapSize=1425,i=NULL, o=NULL, label=NULL, labels=NULL, ppm=NULL, rtwindow=NULL, combination=NULL, h=NULL, v=NULL, sampleList=NULL, nSlaves = 1, outputfolder="combined")
{
	version.1 <- get("version.1",envir=.GlobalEnv)

	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")

	if (!is.null(i))
	{
		if (version.1==TRUE)
		{
			mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
		} else
		{
			mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch_2.0.jar", sep="")
		}
		## setup the tool
		tool <- paste(mzmatch, "mzmatch.ipeak.Combine")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
		if (!is.null(o))
		tool <- paste(tool, "-o", o)
		if (!is.null(label))
		tool <- paste(tool, "-label", label)
		if (!is.null(labels))
		tool <- paste(tool, "-labels", labels)
		if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
		if (!is.null(rtwindow))
		tool <- paste(tool, "-rtwindow", rtwindow)
		if (!is.null(combination))
		tool <- paste(tool, "-combination", combination)
		if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
		if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

		system(tool)
	} else
	{
		MainClasses <- levels(as.factor(sampleList$sampleClass))
		if (!file.exists(outputfolder))
		{
			dir.create (outputfolder)
		}

		combineFunction <- function (fnum)
		{
			FILESf <- sampleList$outputfilenames[which(sampleList$sampleClass==MainClasses[fnum])]
			OUTPUTf <- paste(outputfolder,"/",MainClasses[fnum],".peakml",sep="")
			mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch_2.0.jar", sep="")
			## setup the tool
			tool <- paste(mzmatch, "mzmatch.ipeak.Combine")
			i=paste(FILESf, collapse=",")
			tool <- paste(tool, "-i", i)
			o=OUTPUTf
			tool <- paste(tool, "-o", o)
			label=MainClasses[fnum]
			tool <- paste(tool, "-label", label)
			if (!is.null(labels))
			tool <- paste(tool, "-labels", labels)
			if (!is.null(ppm))
			tool <- paste(tool, "-ppm", ppm)
			if (!is.null(rtwindow))
			tool <- paste(tool, "-rtwindow", rtwindow)
			if (!is.null(combination))
			tool <- paste(tool, "-combination", combination)
			if (!is.null(h) && h==T)
			tool <- paste(tool, "-h")
			if (!is.null(v) && v==T)
			tool <- paste(tool, "-v")
			cat(label,"\n")
			system(tool)
		}

		if (nSlaves>1)
		{
			cl <- makeCluster (nSlaves)
			envname <- environment()
			clusterExport (cl, list=c("combineFunction","MainClasses","sampleList","labels","ppm","rtwindow","combination","java","h","v"), envir=envname)
			system.time(clusterApply(cl,1:length(MainClasses),combineFunction))
			stopCluster(cl)
		} else

		{
			for (fnum in 1:length(MainClasses))
			{
				combineFunction(fnum)
			}
		}
	}
}
