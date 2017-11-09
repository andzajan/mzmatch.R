\name{PeakML.Isotope.getTrendMtx}
\Rdversion{1.1}
\alias{PeakML.Isotope.getTrendMtx}
\title{PeakML.Isotope.getTrendMtx}
\description{
	Returns the trend of labelled and unlabelled carbons in the compound. Convert output of PeakML.Isotope.getTrendList as a matrix
}
\usage{
	PeakML.Isotope.getTrendMtx(trendList, sampleGroups)
}

\arguments{
	\item{trendList}{
		The trend of labelled and unlabelled carbons in the compound. See PeakML.Isotope.getTrendList
	}
	\item{sampleGroups}{
		The sample groups.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getTrendMtx returns the trend of labelled and unlabelled carbons in the compound as a matrix of mean of the sum or max of internsities of all replicates.
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
