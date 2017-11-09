\name{PeakML.Isotope.getMassWindow}
\Rdversion{2}
\alias{PeakML.Isotope.getMassWindow}
\title{PeakML.Isotope.getMassWindow}
\description{
	This function will return mass window of the isotope of the given mass.
}
\usage{
	PeakML.Isotope.getMassWindow(mass, nIsotope, ppm, element)
}

\arguments{
	\item{mass}{
		Unlabelled mass of the compound.
	}
	\item{nIsotope}{
		The isotope of interest.
	}
	\item{ppm}{
		The ppm window required.
	}
	\item{element}{
		The element whose isotopes has to be calculated.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getMassWindow, returns mass window of the isotope of the given mass, as a list.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getMassWindow}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
