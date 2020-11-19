PeakML.Viewer <- function(arch="detect",install.path=NULL, JHeapSize=1024, uninstall=FALSE)
{
	if (is.null(install.path))
	{
		install.path <- paste(find.package("mzmatch.R"),"/PeakMLViewer",sep="")
	}

	dir.create (install.path, showWarnings=FALSE)
	
	if (file.access(install.path)!=0)
	{
		cat ("Folder:", install.path,"can't be created. Please specify writable location as the function argument \"install.path\".","\n")
		stop ()
	}

	if (arch=="detect")
	{
		OS <- R.Version()$platform
		cat ("OS type:",OS,"\n")
		platform <- R.Version()$arch
		if (length(grep("linux",OS))!=0)
		{
			if (platform=="x86_64")
			{
				arch="Linux_64"
				platform=64
			}
		}
		if (length(grep("darwin",OS))!=0)
		{
			if (platform=="x86_64")
			{
				arch="OSX_64"
				platform=64
			} 
		}
		if (length(grep("mingw",OS))!=0)
		{
			if (platform=="x86_64")
			{
				arch="Windows_64"
				platform=64
			} 
		}
	}
	
	viewerfile <- paste(install.path,"/PeakML_Viewer_",arch,"-bit.jar",sep="")
	viewerfile_short <- paste("PeakML_Viewer_",arch,"-bit.jar",sep="")

	if (!file.exists(viewerfile))
	{
		filename <- paste("PeakML_Viewer_",arch,"-bit.jar",sep="")
		download.file(paste("http://sourceforge.net/projects/mzmatch/files/PeakML%20Viewer/",filename,"/download",sep=""),destfile=viewerfile,cacheOK=FALSE,mode="wb")
	}

	## Create settings.xml file
	DBS <- dir(paste(find.package("mzmatch.R"),"/dbs",sep=""),full.names=TRUE)
	settingsfile <- paste(install.path,"/settings.xml",sep="")
	
	databases <- xmlNode("databases")

	nodeGen <- function (nodenum)
	{
		xmlNode ("database",DBS[nodenum])
	}
	dbs <- lapply (1:length(DBS),nodeGen)
	databases <- addChildren(databases,kids=dbs)

	settings <- xmlNode("settings")
	smooth <- xmlNode("smooth","false")
	settings <- xmlNode ("settings",smooth,databases)
	pkmlviewer <- xmlNode ("peakmlviewer", settings)

	saveXML(doc=pkmlviewer,file=settingsfile,prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n \n")

	## Generate command to star PeakML Viewer
	## java -d32 -XstartOnFirstThread -jar PeakMLViewerOsX.jar
	if (arch == "Windows_64")
	{
		java <- "java"
	}
	
	if (arch=="OSX_64")
	{
		java <- paste(java, " -XstartOnFirstThread",sep="")
	}
	
	start.command <- paste(java," -Xms",JHeapSize,"m -Xmx",JHeapSize,"m -jar ",viewerfile_short,sep="")
	
	## Temporary change working folder to PeakMl viewer path, to load also proper settings.xml file. Satr viewer and change folder back.
	if (uninstall==FALSE)
	{	
		currentfolder <- getwd ()
		setwd (install.path)
		if (arch=="Windows_64")
		{
			shell (start.command,wait=FALSE)
		} else
		{
			system (start.command,wait=FALSE,show.output.on.console = FALSE)
		}
		setwd (currentfolder)
	}
	if (uninstall==TRUE)
	{
		unlink(install.path, recursive=TRUE)
		cat ("PeakML Viewer files are removed from the given folder.","\n")
	}
}
