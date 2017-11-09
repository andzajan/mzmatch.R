\name{PeakML.Methods.getStdErr}
\Rdversion{1.1}
\alias{PeakML.Methods.getStdErr}
\title{PeakML.Methods.getStdErr}
\description{
	Returns the standard error based on the standard deviation.
}
\usage{
	PeakML.Methods.getStdErr(x)
}

\arguments{
	\item{x}{
		A numeric vector whose standard error has to be computed.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Methods.getStdErr returns the standard error based on the standard deviation as a number.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getTrendList}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
