\name{PeakML.BlankFilter}
\Rdversion{1.1}
\alias{PeakML.BlankFilter}
\title{PeakML.BlankFilter}
\description{
	This function compares a signal intensities of "real" samples with signal present in blanks (samples usuallu aquired at the beginnning or between analytical samples to equalibrate a LC system and flush out any residues.)
}
\usage{
PeakML.BlankFilter (filename, ionisation="detect", Rawpath=NULL, outputfile, BlankSample=NULL, IgnoreIntensity=FALSE, detectedInNumOfBlanks=NULL, rtwindow=NULL)
}
\arguments{
\item{filename}{
	File in peakML format to be used as input file.
}
\item{ionisation}{
	Ionisation polarity at which data set is acquired. If set to "detect" wil be automatically detected from the data set.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
}
\item{outputfile}{
	File name of the peakml output file  to be created.
}
\item{BlankSample}{
	Character string fefined a sample groups containing blank samples.
}
\item{IgnoreIntensity}{
	If set to TRUE all peak sets what are present in blank will be removed from a real sample ignoring signal intensity.
}
\item{detectedInNumOfBlanks}{
	Only remove peak set if siganl is detecet in n number of blank samples.
}
\item{rtwindow}{
	Numeric vector of two values, defining retention time window in which to apply filter. Fore examples rtwindow=c(5,20).
}
}

\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
}
\value{
	Output file in peakml file format with signal matching blank filter criteria filtered out.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Plot.RankProducts}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
