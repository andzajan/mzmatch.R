\name{PeakML.Methods.getRawSignals}
\Rdversion{1.1}
\alias{PeakML.Methods.getRawSignals}
\title{PeakML.Methods.getRawSignals}
\description{
	Returns the raw signals for the given sample in the rt and mass window specified as a table with masses, intensities, retentiontimes and scanids along the rows.
}
\usage{
	PeakML.Methods.getRawSignals(mzXMLSrc, sampleName, rtWindow, massWindow, massCorrection)
}

\arguments{
	\item{mzXMLSrc}{
		mzXML raw data source path
	}
	\item{sampleName}{
		Name of the sample
	}
	\item{rtWindow}{
		The retention time window as list [min, max]
	}
	\item{massWindow}{
		The mass window as list [min, max]
	}
	\item{massCorrection}{
		See PeakML.Methods.getMassCorrection()
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Methods.getRawSignals returns a table with masses, intensities, retentiontimes and scanids along the rows.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.GapFiller}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
