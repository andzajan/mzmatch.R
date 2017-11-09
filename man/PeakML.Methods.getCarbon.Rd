\name{PeakML.Methods.getCarbon}
\Rdversion{1.1}
\alias{PeakML.Methods.getCarbon}
\title{PeakML.Methods.getCarbon}
\description{
	This function will return the number of carbons in a given formula.
}
\usage{
	PeakML.Methods.getCarbon(formula)
}

\arguments{
	\item{formula}{
		The formula.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Methods.getCarbon returns the number of carbons in a given formula, as a number.
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
