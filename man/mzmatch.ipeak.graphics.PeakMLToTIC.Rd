\name{mzmatch.ipeak.graphics.PeakMLToTIC}
\alias{mzmatch.ipeak.graphics.PeakMLToTIC}
\title{mzmatch.ipeak.graphics.PeakMLToTIC}

\description{
Creates Total Ion Current graphs from the given PeakML file.
}

\usage{
mzmatch.ipeak.graphics.PeakMLToTIC(JHeapSize=1425, i, o, width, height, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
The input file in the PeakML format. When this option is not set, the input is read from the standard input.
}
\item{o}{
The output file for the graph(s). This option is required to be set.
}
\item{width}{
The width of a graph.
}
\item{height}{
The height of a graph.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Creates Total Ion Current graphs from the given PeakML file. The contents is analyzed and correctly summed to get to the real TIC. When the file contains multiple measurements, a seperate TIC is made for each measurement.
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
