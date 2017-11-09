\name{PeakML.Methods.baseCorrection}
\Rdversion{2}
\alias{PeakML.Methods.baseCorrection}
\title{PeakML.Methods.baseCorrection}
\description{
	Correct the baseline of noisy signals by estimating the trend based on asymmetric least squares
}
\usage{
	PeakML.Methods.baseCorrection(signals, lambda=100)
}

\arguments{
	\item{signals}{
		Intensities that are to be corrected
	}
	\item{lambda}{
		The smoothing parameter, the larger it is, the smoother the trendline will be. See documentation for asysm method from ptw package for more details.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Methods.baseCorrection returns the baseline corrected signals.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Methods.formula2mass}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
