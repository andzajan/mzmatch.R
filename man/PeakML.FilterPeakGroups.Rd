\name{PeakML.FilterPeakGroups}
\Rdversion{1.1}
\alias{PeakML.FilterPeakGroups}
\title{PeakML.FilterPeakGroups}
\description{
 This function will filter inconsistien (if several peak groups ar close to each other within ppm and RT range) peak groups for peakML file. It suggested to use this tool before running gap fillling tool.
}
\usage{
PeakML.FilterPeakGroups (filename, ionisation="detect", Rawpath=NULL, ppm=5, rtwin=60, outputfile)
}
\arguments{
\item{filename}{
	Name ot the input file in peakML format.
}
\item{ionisation}{
	Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
}
\item{ppm}{
	Mass window to much peak grops in ppm.
}
\item{rtwin}{
	Retention time window to match peak groups in seconds.
}
\item{outputfile}{
	Name of the resulting peakML file.
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
