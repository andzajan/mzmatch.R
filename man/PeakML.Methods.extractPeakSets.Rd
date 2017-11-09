\name{PeakML.Methods.extractPeakGroups}
\Rdversion{1.1}
\alias{PeakML.Methods.extractPeakGroups}
\title{PeakML.Methods.extractPeakGroups}
\description{
 This function creates a peakML file of subset of given peak groups.
}
\usage{
  PeakML.Methods.extractPeakGroups (PeakMLData,outputfile,sets)
}
\arguments{
\item{PeakMLData}{
 A R object of the PeakML.Read function output.
}
\item{outputfile}{
	Filename of the output peakML file to be generates
}
\item{sets}{
	Numeric vector of the peak group indiced to be extracted
}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.Methods.getGroupAnnotations function allows to access attributes defined in PeakML file.
}
\value{
	List of 3 data tables.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.Read}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
