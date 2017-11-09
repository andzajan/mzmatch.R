\name{PeakML.Methods.getGroupAnnotations}
\Rdversion{1.1}
\alias{PeakML.Methods.getGroupAnnotations}
\title{PeakML.Methods.getGroupAnnotations}
\description{
  This function allows to read group annotations (identificaton, relation ship etc) which are defined in PeakML file.
}
\usage{
  PeakML.Methods.getGroupAnnotations (PeakMLtree=NULL,filename=NULL)
}
\arguments{
\item{PeakMLtree}{
  For internal use of the function PeakMLtree containing XML tree of PeakMl file contents can be provided. Parameter 'filename' should be set to NULL.
}
\item{filename}{
	The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail. If filenames is provided PeakMLtree argument will be igonred, and XML tree created from given peakml file.
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
	List of the existing group annotations in given PeakML file of XML tree structute.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.xcms.read}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
