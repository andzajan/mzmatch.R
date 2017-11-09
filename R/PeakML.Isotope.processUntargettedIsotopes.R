PeakML.Isotope.processUntargettedIsotopes <- function(peakMLFile, analyse, databases, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, numSlaves, label, exclude_from_plots) {

	molFrames <- list()

	if (analyse == "databases"){
		mFrame <- PeakML.Isotope.ReadDB(databases)

		if (!numSlaves==1){
			# Split the matrix into small chunks and save as a list of dataframes

			startPoints <- seq(1, length(row.names(mFrame)), round(length(row.names(mFrame))/numSlaves))
			for (sp in 1:length(startPoints)){
				start <- startPoints[sp]
				if (!is.na(startPoints[sp+1])){
					end <- startPoints[sp+1]-1
				} else {
					end <- length(row.names(mFrame))
				}
				dFrame <- mFrame[start:end,]
				dFrame$rt <-NA
				dFrame$follow <- NA
				molFrames [[sp]] <- dFrame
			}
		} else {
			dFrame <- mFrame
			dFrame$rt <- NA
			dFrame$follow <- NA
			molFrames [[1]] <- dFrame
		}
	}

	if (analyse == "identified"){
		mFrame <- PeakML.Isotope.DB2Text(peakMLFile, databases)
		if (!numSlaves ==1 ){

			startPoints <- seq(1, length(row.names(mFrame)), round(length(row.names(mFrame))/numSlaves))
			for (sp in 1:length(startPoints)){
				start <- startPoints[sp]
				if (!is.na(startPoints[sp+1])){
					end <- startPoints[sp+1]-1
				} else {
					end <- length(row.names(mFrame))
				}
				dFrame <- mFrame[start:end,]
				dFrame$rt <- NA
				dFrame$follow <- NA
				molFrames [[sp]] <- dFrame
			}
		} else {
			molFrames[[1]] <- mFrame
		}
	}

	if (!numSlaves == 1){
		plotUtgt <- function(index){
			require(mzmatch.R)
			mzmatch.init()
			PeakML.Isotope.plotUntargettedIsotopes (peakMLFile, molFrames[[index]], outDirectory, paste(outFileName, index, sep="_"), layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, label, exclude_from_plots)
		}

		cl <- makeCluster (numSlaves)
		assign ("peakMLFile", peakMLFile, envir=.GlobalEnv)
		assign ("molFrames", molFrames, envir=.GlobalEnv)
		assign ("outDirectory", outDirectory, envir=.GlobalEnv)
		assign ("outFileName", outFileName, envir=.GlobalEnv)
		assign ("layoutMtx", layoutMtx, envir=.GlobalEnv)
		assign ("ppm", ppm, envir=.GlobalEnv)
		assign ("stdRTWindow", stdRTWindow, envir=.GlobalEnv)
		assign ("sampleNames",sampleNames, envir=.GlobalEnv)
		assign ("peakDataMtx", peakDataMtx, envir=.GlobalEnv)
		assign ("chromDataList", chromDataList, envir=.GlobalEnv)
		assign ("phenoData", phenoData, envir=.GlobalEnv)
		assign ("sampleGroups", sampleGroups, envir=.GlobalEnv)
		assign ("plotOrder",plotOrder, envir=.GlobalEnv)
		assign ("mzXMLSrc", mzXMLSrc, envir=.GlobalEnv)
		assign ("fillGaps",fillGaps, envir=.GlobalEnv)
		assign ("massCorrection",massCorrection , envir=.GlobalEnv)
		assign ("useArea",useArea, envir=.GlobalEnv)
		assign ("filterStringency",filterStringency, envir=.GlobalEnv)
		assign ("baseCorrection",baseCorrection, envir=.GlobalEnv)
		assign ("label",label, envir=.GlobalEnv)
                assign ("exclude_from_plots", label, envir=.GlobalEnv)
		vars <- list("peakMLFile", "molFrames", "outDirectory", "outFileName", "layoutMtx", "ppm", "stdRTWindow", "sampleNames", "peakDataMtx", "chromDataList", "phenoData", "sampleGroups", "plotOrder", "mzXMLSrc", "fillGaps", "massCorrection", "useArea", "filterStringency", "baseCorrection", "label", "exclude_from_plots")
		clusterExport (cl, list=vars)
		system.time(clusterApply(cl, 1:length(molFrames), plotUtgt))
		stopCluster(cl)

	} else {
		PeakML.Isotope.plotUntargettedIsotopes (peakMLFile, molFrames[[1]], outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, label,  exclude_from_plots)

	}

}
