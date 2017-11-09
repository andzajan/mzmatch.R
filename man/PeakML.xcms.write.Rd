\name{PeakML.xcms.write}
\Rdversion{1.1}
\alias{PeakML.xcms.write}
\title{PeakML.xcms.write}
\description{
  This function writes the given XCMS set to the given filename in the PeakML format. This file format
  supports the mass traces (features) and the sets-structure utilized by XCMS and as such provides a
  fully reproducible representation of the contents of an XCMS set structure.
}
\usage{
  PeakML.xcms.write(xset, filename, ionisation="detect")
}
\arguments{
\item{xset}{
  The XCMS set structure containing the data to be written to the PeakML file.
}
\item{filename}{
  The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail.
}
\item{ionisation}{
  Character value to be set on to "positive" if positive ionisation mode data are processed. And "negative" if negative ionisaton mode was used. Default value is seto to "detect", ionisation mode will be detected from the RAW data file.
}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.xcms.write function extends the XCMS package with write functionality for exporting data to other applications.
  For example, the results of an XCMS step can then be loaded in the user interface environment PeakML Viewer and manipulated
  there.
}
\value{
  This function returns no value.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.xcms.read}}
}

\keyword{mzmatch.R}
\keyword{ XCMS }
