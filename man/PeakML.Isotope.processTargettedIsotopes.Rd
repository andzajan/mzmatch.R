\name{PeakML.Isotope.processTargettedIsotopes}
\Rdversion{1.1}
\alias{PeakML.Isotope.processTargettedIsotopes}
\title{PeakML.Isotope.processTargettedIsotopes}
\description{
	Read the compound id, name and formula from a file and process the plot for all isotopes for that compound. The resulting pdf is stored in the given filename at the given location.
}
\usage{
	PeakML.Isotope.processTargettedIsotopes(molFormulaFile, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, fillGaps, massCorrection, useArea, baseCorrection, label, exclude_from_plots)
}

\arguments{
	\item{molFormulaFile}{
		A tab separated file containing compound id, compound name, compound formula, followCarbon (the labelling of the carbon to be tracked) along the rows.
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
	\item{baseCorrection}{
		Specifies if baseCorrection has to be applied for all peaks. This is especially helpfull while analysing noisy peaks. Leave it to TRUE in most cases.
	}
	\item{label}{
		Specify the isotope used for labelling as a string. The available options are C13, N15, etc. More will be added once tested with with some real dataset.
	}
    \item{exclude_from_plots}{
        Exclude the groups in this vector from the plots in the output.
    }

}


\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.processTargettedIsotopes generated a pdf document containing the plots of isotopes of each compound listed in the file (molFormulaFile) specifed.
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
