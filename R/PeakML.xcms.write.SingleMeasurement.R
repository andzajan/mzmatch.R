PeakML.xcms.write.SingleMeasurement <- function(xset, filename, ionisation="detect", addscans=5, ppm=5, writeRejected=FALSE, ApodisationFilter=TRUE, nSlaves=1)
{
	version.1 <- get("version.1",envir=.GlobalEnv)
	if (version.1==TRUE)
	{
		PeakML.xcms.write.SingleInstance (xset=xset, outputfile=filename, ionisation=ionisation, addscans=addscans, ppm=ppm, writeRejected=writeRejected, ApodisationFilter=ApodisationFilter)
	} else
	{
		#Dirty hack to get filenames right if samples does not follow alphabetic order in the sampleList.
		

		xseto <- split (xset,xset@filepaths)
		#filename <- sort (filename)
		xsetFunction <- function (i)
		{
			require (mzmatch.R)
			mzmatch.init (version.1=version.1)
			PeakML.xcms.write.SingleInstance(xset=xseto[[i]], outputfile=filename[i],ionisation=ionisation,addscans=addscans,ppm=ppm,writeRejected=writeRejected,ApodisationFilter=ApodisationFilter)
		}

		if (nSlaves>1)
		{
			cl <- makeCluster (nSlaves)
			envname <- environment()
			clusterExport (cl, list=c("xsetFunction","xseto","filename", "ionisation","addscans","ppm","writeRejected","ApodisationFilter","version.1"), envir=envname)
			system.time(clusterApply(cl,1:length(filename),xsetFunction))
			stopCluster(cl)
		} else

		{
			for (i in 1:length(filename))
			{
				PeakML.xcms.write.SingleInstance(xset=xseto[[i]], outputfile=filename[i],ionisation=ionisation,addscans=addscans,ppm=ppm,writeRejected=writeRejected,ApodisationFilter=ApodisationFilter)
			}
		}
	}
cat ("Use version 1: ",version.1,"\n")
}

