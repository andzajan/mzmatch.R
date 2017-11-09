mzmatch.ipeak.filter.SimpleFilter <- function(JHeapSize=1425,i=NULL, o=NULL, rejected=NULL, databases=NULL, ppm=NULL, n=NULL, offset=NULL, mindetections=NULL, minscanid=NULL, maxscanid=NULL, minretentiontime=NULL, maxretentiontime=NULL, minmass=NULL, maxmass=NULL, minintensity=NULL, maxintensity=NULL, annotations=NULL, h=NULL, v=NULL)
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
	tool <- paste(mzmatch, "mzmatch.ipeak.filter.SimpleFilter")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(rejected))
		tool <- paste(tool, "-rejected", rejected)
	if (!is.null(databases))
		tool <- paste(tool, "-databases", databases)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(n))
		tool <- paste(tool, "-n", n)
	if (!is.null(offset))
		tool <- paste(tool, "-offset", offset)
	if (!is.null(mindetections))
		tool <- paste(tool, "-mindetections", mindetections)
	if (!is.null(minscanid))
		tool <- paste(tool, "-minscanid", minscanid)
	if (!is.null(maxscanid))
		tool <- paste(tool, "-maxscanid", maxscanid)
	if (!is.null(minretentiontime))
		tool <- paste(tool, "-minretentiontime", minretentiontime)
	if (!is.null(maxretentiontime))
		tool <- paste(tool, "-maxretentiontime", maxretentiontime)
	if (!is.null(minmass))
		tool <- paste(tool, "-minmass", minmass)
	if (!is.null(maxmass))
		tool <- paste(tool, "-maxmass", maxmass)
	if (!is.null(minintensity))
		tool <- paste(tool, "-minintensity", minintensity)
	if (!is.null(maxintensity))
		tool <- paste(tool, "-maxintensity", maxintensity)
	if (!is.null(annotations))
		tool <- paste(tool, "-annotations", annotations)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

	system(tool)
}
