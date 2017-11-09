\name{mzmatch.ipeak.util.BreakupRelated}
\alias{mzmatch.ipeak.util.BreakupRelated}
\title{mzmatch.ipeak.util.BreakupRelated}

\description{
Splits the derivative clusters found by mzmatch.
}

\usage{
mzmatch.ipeak.util.BreakupRelated(JHeapSize=1425, i, o, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
The input file (in PeakML format), resulting from the tool mzmatch.ipeak.sort.RelatedPeaks. When this option is not set the input is read from the standard input.
}
\item{o}{
Specifies the output directory where all the PeakML files are written with the clusters of derivatives.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Splits the derivative clusters found by mzmatch.ipeak.sort.RelatedPeaks into separate files. Each file is named with the relation.id found in the annotations.
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
