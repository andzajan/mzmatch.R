\name{PeakML.Plot.Chromatograms}
\Rdversion{1.1}
\alias{PeakML.Plot.Chromatograms}
\title{PeakML.Plot.Chromatograms}
\description{
	Simple method to create a EIC plot for different sample groups.
}
\usage{
PeakML.Plot.Chromatograms (PeakMLdata,groupid,sampleClasses=NULL, xaxis=TRUE, Title=NULL, mar=c(2,1,0,1))
}
\arguments{
\item{PeakMLdata}{
	PeakML.Read function output object of the peakml file.
}
\item{groupid}{
	Number of the peak set to plot
}
\item{sampleClasses}{
	Character vector indicating which sample groups should be shown on the plot. If set to NULL all sample groups from the peakML file will be plotted.
}
\item{xaxis}{
	If set to "TRUE" x-axis (RT) will be shown on the plot.
}
\item{Title}{
	Character string to show as a title of the plot. Default is "NULL" (empty).
}
\item{mar}{
	Parameter to define a plot margins, see a help(par) for more details.
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
