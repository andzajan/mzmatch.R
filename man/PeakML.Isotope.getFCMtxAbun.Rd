\name{PeakML.Isotope.getFCMtxAbun}
\Rdversion{1.1}
\alias{PeakML.Isotope.getFCMtxAbun}
\title{PeakML.Isotope.getFCMtxAbun}
\description{
	Returns the trend of labelled or unlabelled carbons in the compound as specified by the followCarbon parameter.
}
\usage{
	PeakML.Isotope.getFCMtxAbun(trendList, sampleGroups, followCarbon, element)
}

\arguments{
	\item{trendList}{
		The trend of labelled and unlabelled carbons in the compound. See PeakML.Isotope.getTrendList.
	}
	\item{sampleGroups}{
		The sample groups.
	}
	\item{followCarbon}{
		The number of the carbon in a compound whose trend is to be followed.
	}
	\item{element}{
		The element such as C, N that is labelled in the data being analysed.
	}

}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getFCMtxAbun returns the trend of labelled or unlabelled carbons in the compound as specified by the followCarbon parameter as a matrix.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getTrendList} \link{PeakML.Isotope.getChromData} \link{PeakML.Isotope.getIsotopes}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
