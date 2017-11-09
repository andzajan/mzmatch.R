\name{PeakML.Isotope.plotUntargettedIsotopes}
\Rdversion{1.1}
\alias{PeakML.Isotope.plotUntargettedIsotopes}
\title{PeakML.Isotope.plotUntargettedIsotopes}
\description{
	Retreives the list of all identified compounds from the PeakML file and generate isotope plots for all of them, if present. The output is stored under the given filename at the given location.
}
\usage{
	PeakML.Isotope.plotUntargettedIsotopes(peakMLFile, molFrame, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, filterStringency, baseCorrection, label, exclude_from_plots)
}

\arguments{
	\item{peakMLFile}{
		The name of the PeakML file which contains the combined peaks of all samples that are filtered, gap filled and identified. See the advanced metabolomics data processing tutorial at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.html for more information.
	}
	\item{molFrame}{
		A list containing containing data.frames of compound id, name and formula.
	}
	\item{outDirectory}{
		The output directory.
	}
	\item{outFileName}{
		The output filename
	}
	\item{layoutMtx}{
		The layout matrix that specifies the layout of the plot to be generated. Leave this as NULL if the number of sample to be plotted is less than or equal to 10. See layout in R more details.
	}
	\item{ppm}{
		The ppm to specify the mass window for searching the isotopes.
	}
	\item{stdRTWindow}{
		The standard retention time  window. Leave it to NULL if no standard retention time is available for the compound.
	}
	\item{sampleNames}{
		The names of samples.
	}
	\item{peakDataMtx}{
		The peak data matrix. See PeakML.Methods.getPeakData.
	}
	\item{chromDataList}{
		The list of chromatograms. See PeakML.Methods.getChromData.
	}
	\item{phenoData}{
		The phenotype data of the samples. See PeakML.Methods.getPhenoData.
	}
	\item{sampleGroups}{
		The sample groups.
	}
	\item{plotOrder}{
		The order in which samples has to be plotted.
	}
	\item{mzXMLSrc}{
		The path to the raw data mzXML files.
	}
	\item{fillGaps}{
		"ALLPEAKS" if all the peaks are to gap filled, requied to fill missing signals in some part of the curve. "MISSINGPEAKS", if only the missing peaks have to be filled Or "NONE" if no peaks have to be gap filled. Defaults to "ALLPEAKS".
	}
	\item{massCorrection}{
		Requied if gap fill is TRUE. This is to correct the natural mass with the ionised mass. See PeakML.Methods.getMassCorrection(). see PeakML.Methods.getMassCorrection. 
	}
	\item{useArea}{
		TRUE if sum of intensities has to be used instead of the maximum of intensities for peak intensity calculations.
	}
	\item{filterStringency}{
		This is the percentage of samples that needs to be labelled. It can be left to "NULL" if the filter has to be turned off.
	}
	\item{baseCorrection}{
		TRUE or FALSE based on the requirement of base correction for noisy peaks. Enabling this will make the processing slower.
	}
	\item{label}{
		This is the label to be identified. 1 = C13; 2= N15. More will be added as and when implemented.
	}
	\item{exclude_from_plots}{
                Exclude sample groups in this list from the plots. Can be used for removing blanks and standards from the result plots.
        }
}


\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.plotUntargettedIsotopes generates a pdf document containing the plots of isotopes of each identified compound in the PeakML file specified.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getIsotopes}, \link{PeakML.Isotope.getChromData} \link{PeakML.Isotope.plotSamples}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
