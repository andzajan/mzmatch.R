\name{PeakML.Methods.getCompleteTable}
\Rdversion{1.1}
\alias{PeakML.Methods.getCompleteTable}
\title{PeakML.Methods.getCompleteTable}
\description{
 This function creates a list of intensities, mass and retention times tables. Data in this format are suitable for PCA analysis or similar statistical routines.
}
\usage{
  PeakML.Methods.getCompleteTable (PeakMLData,sumintensity=FALSE)
}
\arguments{
\item{PeakMLData}{
 A R object of the PeakML.Read function output.
}
\item{sumintensity}{
	If set to "TRUE" sum of all intensities is used as a matrix value for intensities table instead of intensity the max intensity of thr peak.
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
