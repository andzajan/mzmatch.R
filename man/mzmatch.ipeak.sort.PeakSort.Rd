\name{mzmatch.ipeak.sort.PeakSort}
\alias{mzmatch.ipeak.sort.PeakSort}
\title{mzmatch.ipeak.sort.PeakSort}

\description{
This tools sorts the contents of a PeakML and writes the result to a new file.
}

\usage{
mzmatch.ipeak.sort.PeakSort(JHeapSize=1425, i, o, sort, annotation, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). The only allowed format is PeakML and when it is not set the input is read from standard in.
}
\item{o}{
Option for the ouput file. The file is written in the PeakML file format and peaks that passed the defined filter are saved here. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{sort}{
- scanascending|scandescending
  can be used to sort the list of peaks on retention time
- massascending|massdescending
  can be used to sort the list of peaks on mass
- intensityascending|intensitydescending
  can be used to sort the list of peaks on intensity
- annotationascending|annotationdescending  can be used to sort the list of peaks on an annotation
}
\item{annotation}{
When either annotationascending or annotationdescending has been chosen for the sort option, this option needs to be filled in for the label of the annotation to sort on.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
This tools sorts the contents of a PeakML and writes the result to a new file. A number of sorting mechanisms is supported: intensity, mass, scan and annotation. When sorting on annotation the -annotation options needs to be used to indicate on which label to sort the list.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.


}

\author{
RA Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
\code{\link{PeakML.xcms.read}}
}

\keyword{mzMatch}
\keyword{peakML}
\keyword{java tools}
