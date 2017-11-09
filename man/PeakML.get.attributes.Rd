\name{PeakML.get.attributes}
\Rdversion{1.1}
\alias{PeakML.get.attributes}
\title{PeakML.get.attributes}
\description{
  This function allows to read attributes which are defined in PeakML file.
}
\usage{
  PeakML.get.attributes(filename,attribute="MeasurementNames",annotations=NULL)
}
\arguments{
\item{filename}{
	The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail.
}
\item{attribute}{
  Name of the attribute from PeakML file header which should be read. Currently supports 'MeasurementNames': will return character vector of sample names present in PeakML file. "IntensitiesTable":  will return matrix of peak intensities - samples as rows, peaks as columns. 'NA' values is set if peak is not present in this sample. "CompleteTable": will return list which contains three matrices, first contain peak intensity values, second measured mass and third contains retention time. "getGroupAnnotation" and "getAnnotation": returns character vector of attribute values in peakml groups and masschromatograms levels. Values is set to 'NA' if annoation is not set for current peak or group.
}
\item{annotations}{
	Character vector of annotation names which sould should be extracted, if function argument attribute is set to 'getGroupAnnotation' or 'getAnnotation'. For example: c("identification","relation.id","ppm","relation.ship") for 'getGroupAnnotations'.
}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.get.attributes function allows to access attributes defined in PeakML file.
}
\value{
  Character vector for "MeasurementNames",matrix  for "IntensitiesTable" or list of attribute values for "getGroupAnnotation" and "getAnnotation" are returned.
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

\keyword{ mzmatch.R }
\keyword{ peakML }
