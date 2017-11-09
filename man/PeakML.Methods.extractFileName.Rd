\name{PeakML.Methods.extractFileName}
\Rdversion{1.1}
\alias{PeakML.Methods.extractFileName}
\title{PeakML.Methods.extractFileName}
\description{
	Simple method to extract eaxct file name from the longer file path. Used to auto generate output file names.
}
\usage{
PeakML.Methods.extractFileName (path)
}
\arguments{
\item{path}{
	Full path to then input file.
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
	Stripped file name which does not contains path to folders and subfolders.
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
