mzmatch.ipeak.align.CowCoda <- function(JHeapSize=1425,i=NULL, o=NULL, ppm=NULL, codadw=NULL, order=NULL, maxrt=NULL, image=NULL, selection=NULL, h=NULL, v=NULL)
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
	tool <- paste(mzmatch, "mzmatch.ipeak.align.CowCoda")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(codadw))
		tool <- paste(tool, "-codadw", codadw)
	if (!is.null(order))
		tool <- paste(tool, "-order", order)
	if (!is.null(maxrt))
		tool <- paste(tool, "-maxrt", maxrt)
	if (!is.null(image))
		tool <- paste(tool, "-image", image)
	if (!is.null(selection))
		tool <- paste(tool, "-selection", selection)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

	system(tool)
}
