\name{PeakML.Methods.getRawDataPaths}
\Rdversion{1.1}
\alias{PeakML.Methods.getRawDataPaths}
\title{PeakML.Methods.getRawDataPaths}
\description{
	This function returns a list containing samples names and the full path of mzXML files that represent these sample names.
}
\usage{
	PeakML.Methods.getRawDataPaths(PeakMLtree, Rawpath=NULL)
}

\arguments{
	\item{PeakMLtree}{
		XML tree of PeakML file content.
	}
	\item{Rawpath}{
		Defaulted to NULL. A new path can be specified if mzXML files (from which PeakML files were created) are relocated.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
	R list containing the samples names and the full path of mzXML files that represent these sample names.
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getMassCorrection}, \link{PeakML.Methods.getPeakData}, \link{PeakML.Methods.getChromData}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
