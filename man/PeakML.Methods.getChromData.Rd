\name{PeakML.Methods.getChromData}
\Rdversion{1.1}
\alias{PeakML.Methods.getChromData}
\title{PeakML.Methods.getChromData}
\description{

	This function returns a list containing chromatograms for each peak in the peak data matrix.
}
\usage{
	PeakML.Methods.getChromData(filename, PeakMLtree, massCorrection)
}

\arguments{
	\item{filename}{
		Name of the peakML data fall including full path.
	}
	\item{PeakMLtree}{
		XML tree of PeakML file content.
	}
	\item{massCorrection}{
		To correct for masses based on the ionisation detected. 
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
	Chromatograms for peaks in the the peak data matrix, in R list.
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getMassCorrection}, \link{PeakML.Methods.getPeakData}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
