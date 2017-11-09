\name{PeakML.Isotope.getIsotopes}
\Rdversion{1.1}
\alias{PeakML.Isotope.getIsotopes}
\title{PeakML.Isotope.getIsotopes}
\description{
	This function will return location of isotope peaks in the peak data matrix that fall in the given mass window and retention time window (if provided).
}
\usage{
	PeakML.Isotope.getIsotopes(peakDataMtx, mzXMLSrc, sampleNames, label, numElements, metMass, ppm, massCorrection, baseCorrection, stdRT=NULL, stdRTWindow=NULL, fillGaps="ALLPEAKS")
}

\arguments{
	\item{peakDataMtx}{
		The peak data matrix.
	}
	\item{mzXMLSrc}{
		Path to where mzXML files are stored as string.
	}
	\item{sampleNames}{
		Sample names as a vector.
	}
	\item{label}{
		The element whose isotope we are interested in.
	}
	\item{numElements}{
		The number of specified element in the compound plus one (this is to handle the unlabelled as well).
	}
	\item{metMass}{
		The unlabelled mass of the compound in question.
	}
	\item{ppm}{
		The ppm to calculate the mass window.
	}
	\item{massCorrection}{
		Requied if gap fill is TRUE. This is to correct the natural mass with the ionised mass. See PeakML.Methods.getMassCorrection().
	}
	\item{baseCorrection}{
		Specifies if baseCorrection has to be applied for all peaks. This is especially helpfull while analysing noisy peaks. Leave it to TRUE in most cases.
	}
	\item{stdRT}{
		The stand RT if available else provide NULL to this function.
	}
	\item{stdRTWindow}{
		The time in seconds to get the retention time window. Else leave as NULL.
	}
	\item{fillGaps}{
		"ALLPEAKS" if all the peaks are to gap filled, requied to fill missing signals in some part of the curve. "MISSINGPEAKS", if only the missing peaks have to be filled Or "NONE" if no peaks have to be gap filled. Defaults to "ALLPEAKS".
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getIsotopes returns the location of isotope peaks in the peak data matrix that fall in the given mass window and retention time window as a list of the form isotope[[peak_group]][[isotope]][[sample]].
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getMassWindow} \link{PeakML.Isotope.filterPeaks} \link{PeakML.Methods.getRawSignals} \link{PeakML.Methods.getPPMWindow}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
