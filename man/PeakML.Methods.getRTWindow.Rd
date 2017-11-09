\name{PeakML.Methods.getRTWindow}
\Rdversion{1.1}
\alias{PeakML.Methods.getRTWindow}
\title{PeakML.Methods.getRTWindow}
\description{
	Returns the retention time window wrt the time (in seconds) specified.
}
\usage{
	PeakML.Methods.getRTWindow(retTime, time)
}

\arguments{
	\item{retTime}{
		Original retention time in seconds.
	}
	\item{time}{
		Required time window in seconds.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	getRTWindow returns the retention time window wrt to the time specified in seconds, as a list.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Methods.getPPMMass}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
