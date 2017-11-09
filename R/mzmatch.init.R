mzmatch.init <- function (memorysize=1024, version.1=TRUE)
{	
	DBS.path <- paste(find.package("mzmatch.R"),"/dbs",sep="")
	mzmatch.path <- paste(find.package("mzmatch.R"),"/java",sep="")

	# Install required JAR files and XML databeses on the first run
	if (!file.exists(DBS.path))
	{
		current.dir <- getwd ()
		setwd (find.package("mzmatch.R"))
		download.file("http://sourceforge.net/projects/mzmatch/files/mzmatch.R/dbs.zip/download","dbs.zip",cacheOK=FALSE,mode="wb")
		unzip ("dbs.zip")
		file.remove ("dbs.zip")
		setwd (current.dir)
	}
	
	if (version.1==TRUE)
	{
		if (!file.exists(paste(mzmatch.path,"/mzmatch.jar",sep="")))
		{
			current.dir <- getwd ()
			if (!file.exists(mzmatch.path))
			{
				dir.create (mzmatch.path)
			}
			setwd (mzmatch.path)
			download.file("http://sourceforge.net/projects/mzmatch/files/mzmatch.R/mzmatch.jar/download","mzmatch.jar",cacheOK=FALSE,mode="wb")
			setwd (current.dir)
		}
		lib <- paste(mzmatch.path,"/mzmatch.jar",sep="")
	}

	else
	{
		if(!file.exists(paste(mzmatch.path,"/mzmatch_2.0.jar",sep="")))
		{
			current.dir <- getwd ()
			if (!file.exists(mzmatch.path))
			{
				dir.create (mzmatch.path)
			}
			setwd (mzmatch.path)
			download.file("http://sourceforge.net/projects/mzmatch/files/mzmatch.R/mzmatch_2.0.jar/download","mzmatch_2.0.jar",cacheOK=FALSE,mode="wb")
			setwd (current.dir)
		}
		lib <- paste(mzmatch.path,"/mzmatch_2.0.jar",sep="")
	}
	params <- paste("-Xmx",memorysize,"m",sep="")
	.jinit(classpath=lib, parameters=params, force.init=TRUE)
	cat(lib,"\n")
	assign("version.1", version.1, envir=.GlobalEnv)
}
