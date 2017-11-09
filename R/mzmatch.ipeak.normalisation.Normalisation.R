mzmatch.ipeak.normalisation.Normalisation <- function(JHeapSize=1425,i=NULL, o=NULL, type=NULL, v=NULL, h=NULL)
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
	tool <- paste(mzmatch, "mzmatch.ipeak.normalisation.Normalisation")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(type))
		tool <- paste(tool, "-type", type)
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")

	system(tool)
}
