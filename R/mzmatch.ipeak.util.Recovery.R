mzmatch.ipeak.util.Recovery <- function(JHeapSize=1425,i=NULL, o=NULL, recovered=NULL, bins=NULL, ppm=NULL, rtwindow=NULL, h=NULL, v=NULL)
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
	} else
	{
		mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch_2.0.jar", sep="")
	}
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.util.Recovery")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(recovered))
		tool <- paste(tool, "-recovered", recovered)
	if (!is.null(bins))
		tool <- paste(tool, "-bins", bins)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(rtwindow))
		tool <- paste(tool, "-rtwindow", rtwindow)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

	system(tool)
}
