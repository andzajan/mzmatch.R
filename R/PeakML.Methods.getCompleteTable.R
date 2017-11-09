PeakML.Methods.getCompleteTable <- function(PeakMLData,sumintensity=FALSE)
{
	nrpeakgroups <- max(PeakMLData$peakDataMtx[,10])

	nrsamples <- length(PeakMLData$sampleNames)

	intensities_matrix <- matrix(ncol=nrpeakgroups,nrow=nrsamples)
	mass_matrix <- matrix(ncol=nrpeakgroups,nrow=nrsamples)
	rt_matrix <- matrix(ncol=nrpeakgroups,nrow=nrsamples)

	system.time(
	{
	for (i in 1:nrsamples)
	{
		hit <- which(PeakMLData$peakDataMtx[,9]==i)
		whichsets <- PeakMLData$peakDataMtx[hit,10]
		if (sumintensity==FALSE)
		{
			intensities_matrix[i,whichsets] <- PeakMLData$peakDataMtx[hit,8]
		} else
		{
			intensities_matrix[i,whichsets] <- PeakMLData$peakDataMtx[hit,7]
		}
		mass_matrix[i,whichsets] <- PeakMLData$peakDataMtx[hit,1]
		rt_matrix[i,whichsets] <- PeakMLData$peakDataMtx[hit,4]
	}
	})

	out <- list(intensities_matrix,mass_matrix,rt_matrix)
	names(out) <- c("Intensities","Masses","Retentiontimes")
	out
}
