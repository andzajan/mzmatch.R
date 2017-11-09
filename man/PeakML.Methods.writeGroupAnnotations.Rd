\name{PeakML.Methods.writeGroupAnnotations}
\Rdversion{1.1}
\alias{PeakML.Methods.writeGroupAnnotations}
\title{PeakML.Methods.writeGroupAnnotations}
\description{
  This function allows to write group annotations (identificaton, relation ship etc) which are defined in PeakML file.
}
\usage{
  PeakML.Methods.writeGroupAnnotations (project, GroupAnnotations)
}
\arguments{
\item{project}{
	The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail. If filenames is provided PeakMLtree argument will be igonred, and XML tree created from given peakml file.
}
\item{GroupAnnotations}{
  List of annotations which should be written to PeakML file. List names should contain annotation names. And each list should contain a vector of the length equal to number of peak groups. NA values should be used if annotation is not present in given group.
}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
}
\value{
	No value is returned.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getGroupAnnotations}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
