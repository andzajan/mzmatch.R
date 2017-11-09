\name{PeakML.Plot.RelatedPeaks}
\Rdversion{1.1}
\alias{PeakML.Plot.RelatedPeaks}
\title{PeakML.Plot.RelatedPeaks}
\description{
	Function to generate an overview pdf document containing summary information shown fragment ad other related signals. Peaks shapes and groups annotations. Including putative identifications.
}
\usage{
PeakML.Plot.RelatedPeaks (filename, ionisation="detect", Rawpath=NULL, DBS=dir(paste(find.package("mzmatch.R"), "/dbs", sep=""),full.names=TRUE),outputfile,sampleClasses=NULL)
}
\arguments{
\item{filename}{
	File in peakML format to be used as input file.
}
\item{ionisation}{
	Ionisation mode of the data set.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
}
\item{DBS}{
	Character vector containing names of the compound data bases files in XML format.
}
\item{outputfile}{
	File name of the pdf output to be created.
}
\item{sampleClasses}{
	Character vector indicating whihc sample groups should be shown on the plot. If set to NULL all sample groups from the peakML file will be plotted.
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
	Plot to the current graphical device.
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
