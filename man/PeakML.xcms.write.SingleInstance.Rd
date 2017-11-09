\name{PeakML.xcms.write.SingleInstance}
\Rdversion{1.1}
\alias{PeakML.xcms.write.SingleInstance}
\title{PeakML.xcms.write.SingleInstance}
\description{
  This is a single core function to write the given XCMS set (for a single measurement) to the given filename in the PeakML format. This file format
  supports the mass traces (features) and the sets-structure utilized by XCMS and as such provides a
  fully reproducible representation of the contents of an XCMS set structure.
}
\usage{
  PeakML.xcms.write.SingleInstance(xset, outputfile,ionisation=ionisation,addscans=addscans,ppm=ppm,writeRejected=writeRejected,ApodisationFilter=ApodisationFilter)
}
\arguments{
\item{xset}{
  The XCMS set structure containing the data to be written to the PeakML file.
}
\item{outputfile}{
  The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail.
}
\item{ionisation}{
  Character value to be set on to "positive" if positive ionisation mode data are processed. And "negative" if negative ionisaton mode was used. Default value is seto to "detect", ionisation mode will be detected from the RAW data file.
}
\item{addscans}{
  Number of scans which to add to the tails of extracted mass chromatograms on the each side.
}
\item{ppm}{
  Mass window in ppm which to use when mass chromatogram is extracted from RAW data file, before exporting to peakml file.
}
\item{writeRejected}{
  If values is set to TRUE, then rejected chromatograms will be writen in separate file with prefix: *_rejected.peakml. Criteria for adding mass chromatogram to the list of rejected are scan count less than 3 scans or if all intensities in RT region detected by XCMS peak picking algorithm are 0.
}

\item{ApodisationFilter}{
  If values is set to TRUE, peak set will be filtered on centroid artefacts which are common in Orbitrap data sets. Every mass within 0.9 mass units and 5s retention time window is selected. Peak with the highest signal intensity is retained.
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
  \code{\link{PeakML.xcms.write.SingleMeasurement}}
}

\keyword{mzmatch.R}
\keyword{ XCMS }
