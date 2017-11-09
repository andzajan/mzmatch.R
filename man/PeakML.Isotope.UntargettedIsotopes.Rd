\name{PeakML.Isotope.UntargettedIsotopes}
\Rdversion{1.1}
\alias{PeakML.Isotope.UntargettedIsotopes}
\title{PeakML.Isotope.UntargettedIsotopes}
\description{
	Plot isotopes of all identified compounds in the PeakML file.
}
\usage{
	PeakML.Isotope.UntargettedIsotopes(baseDir, outFileName, mzXMLSrc=NULL, outDirectory = "untargettedIsotops", peakMLFile="final_combined_related_identified.peakml", analyse = "databases", databases = c("kegg"), sampleGroups = NULL, layoutMtx = NULL, ppm = 3, trendPlots = NULL, fillGaps = "ALLPEAKS", useArea = FALSE, stdRTWindow = NULL, filterStringency=30, baseCorrection=FALSE, numSlaves=1, label=1, exclude_from_plots = NULL)
}

\arguments{
	\item{baseDir}{
		The parent directory where PeakML files are located.
	}
	\item{outFileName}{
		The output filename
	}
	\item{mzXMLSrc}{
		The path to the raw data mzXML files.
	}
	\item{outDirectory}{
		The output directory.
	}
	\item{peakMLFile}{
		The name of the PeakML file which contains the combined peaks of all samples that are filtered, gap filled and identified. See the advanced metabolomics data processing tutorial at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.html for more information.
	}
	\item{analyse}{
		Indicates what type of analysis you want to run. The tree possible values at the moment are "identified", "databases" and "allmasses". "identified" and "allmasses" can be used to run an untargetted isotope identification procedure on the given PeakML file looking for all the identified compounds or all masses detected. "databases" can be used to run entire databases such as KEGG and HMDB.
	}
	\item{databases}{
		A vector containing databases as numbers or strings based on which initial annotations were made or the databases that has to be used for identification. For more information, see 'Identify peaks from databases' section at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.php
	}
	\item{sampleGroups}{
		The sample groups.
	}
	\item{layoutMtx}{
		The layout matrix that specifies the layout of the plot to be generated. Leave this as NULL if the number of sample to be plotted is less than or equal to 10. See layout in R more details.
	}
	\item{ppm}{
		The ppm to specify the mass window for searching the isotopes.
	}
	\item{trendPlots}{
		Leave it to NULL if the default trend plots have to be plots. Else specify the name of the custom trend plots as a vector. Corresponding plot function has to be incorporated into the PeakML.Isotope.plotSamples function.
	}
	\item{fillGaps}{
		"ALLPEAKS" if all the peaks are to gap filled, requied to fill missing signals in some part of the curve. "MISSINGPEAKS", if only the missing peaks have to be filled Or "NONE" if no peaks have to be gap filled. Defaults to "ALLPEAKS".
	}
	\item{useArea}{
		TRUE if sum of intensities has to be used instead of the maximum of intensities for peak intensity calculations.
	}
	\item{stdRTWindow}{
		The standard retention time  window. Leave it to NULL if no standard retention time is available for the compound.
	}
	\item{filterStringency}{
		This is the percentage of samples that needs to be labelled. It can be left to "NULL" if the filter has to be turned off.
	}
	\item{baseCorrection}{
		TRUE or FALSE based on the requirement of base correction for noisy peaks. Enabling this will make the processing slower.
	}
	\item{numSlaves}{
		This parameter can be set to indicate the number of CPU cores that can be used for isotope identification and quantification.
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
	PeakML.Isotope.UntargettedIsotopes generates a pdf document containing the plots of isotopes of each identified compound in the PeakML file specified.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getIsotopes} \link{PeakML.Isotope.getChromData} \link{PeakML.Isotope.plotSamples} \link{PeakML.Isotope.processTargettedIsotopes}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
