\name{PeakML.Methods.PeakShapeCorrelation}
\Rdversion{1.1}
\alias{PeakML.Methods.PeakShapeCorrelation}
\title{PeakML.Methods.PeakShapeCorrelation}
\description{
 This function will calculate peak correlation between two chhromatograms.
}
\usage{
PeakML.Methods.PeakShapeCorrelation (chromatogram1,chromatogram2)
}
\arguments{
\item{chromatogram1}{
	Chromatogram from PeakML.Read$chromDataList
}
\item{chromatogram2}{
	Chromatogram from PeakML.Read$chromDataList
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
	Vector containing correlation value and p.valu from correlation test.
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
