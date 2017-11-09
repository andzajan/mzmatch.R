mzmatch.ipeak.normalisation.VanDeSompele <- function(JHeapSize=1425,i=NULL, o=NULL, selection=NULL, selection_normalized=NULL, img=NULL, factors=NULL, ppm=NULL, database=NULL, v=NULL, h=NULL)
{
	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")
	mzmatch <- paste(java, " ", find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.normalisation.VanDeSompele")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(selection))
		tool <- paste(tool, "-selection", selection)
	if (!is.null(selection_normalized))
		tool <- paste(tool, "-selection_normalized", selection_normalized)
	if (!is.null(img))
		tool <- paste(tool, "-img", img)
	if (!is.null(factors))
		tool <- paste(tool, "-factors", factors)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(database))
		tool <- paste(tool, "-database", database)
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")

	system(tool)
}
