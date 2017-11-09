PeakML.Methods.combinePeakMLfiles <- function(PeakMLfiles, outputfile, ionisation="detect",Rawpath=NULL)
{
	PeakMLtree <- PeakML.Read (PeakMLfiles[1])

	system.time({
	for (filen in 2:length(PeakMLfiles))
	{
		cat (filen, " -- ")
		PeakMLtree2 <- PeakML.Read (PeakMLfiles[filen])
		PeakMLtree2$peakDataMtx[,10] <- PeakMLtree2$peakDataMtx[,10]+max(PeakMLtree$peakDataMtx[,10])
		PeakMLtree$peakDataMtx <- rbind(PeakMLtree$peakDataMtx,PeakMLtree2$peakDataMtx)
		PeakMLtree$chromDataList <- c(PeakMLtree$chromDataList,PeakMLtree2$chromDataList)
		PeakMLtree$sampleGroups <- c(PeakMLtree$sampleGroups, PeakMLtree2$sampleGroups)
		for (annot in 1:length(PeakMLtree$GroupAnnotations))
		{
			PeakMLtree$GroupAnnotations[[annot]] <- c(PeakMLtree$GroupAnnotations[[annot]],PeakMLtree2$GroupAnnotations[[annot]])
		}
	}
	})
	
	system.time(PeakML.Write (peakMLdata=PeakMLtree, outFileName=outputfile))
}
