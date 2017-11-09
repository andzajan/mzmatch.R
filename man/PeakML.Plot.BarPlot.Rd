\name{PeakML.Plot.BarPlot}
\Rdversion{1.1}
\alias{PeakML.Plot.BarPlot}
\title{PeakML.Plot.BarPlot}
\description{
	Simple method to create a barplot with peak set intensities for different sample groups.
}
\usage{
PeakML.Plot.BarPlot (PeakMLdata,groupid,sampleClasses=NULL)
}
\arguments{
\item{PeakMLdata}{
	PeakML.Read function output object of the peakml file.
}
\item{groupid}{
	Number of the peak set to plot
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
