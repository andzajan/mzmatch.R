\name{PeakML.Methods.combinePeakMLfiles}
\Rdversion{1.1}
\alias{PeakML.Methods.combinePeakMLfiles}
\title{PeakML.Methods.combinePeakMLfiles}
\description{
 This function combines several peakML files in one peakML file.
}
\usage{
  PeakML.Methods.combinePeakMLfiles (PeakMLfiles,outputfile,ionisation="detect",Rawpath=NULL)
}
\arguments{
\item{PeakMLfiles}{
	Caharacter vecoter of the file names of file which should be combined.
}
\item{outputfile}{
	Name of the resulting peakML file.
}
\item{ionisation}{
	Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
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
