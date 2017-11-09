\name{mzmatch.ipeak.util.ExtractSets}
\alias{mzmatch.ipeak.util.ExtractSets}
\title{mzmatch.ipeak.util.ExtractSets}

\description{
This tool extracts the entries of the given setnames (stored in the header) from the given PeakML file and stores the result in the output file.
}

\usage{
mzmatch.ipeak.util.ExtractSets(JHeapSize=1425, i, o, setnames, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
The input file (in PeakML format). When this option is not set the input is read from the standard input.
}
\item{o}{
The output file (in PeakML format). When this option is not set the output is written to the standard output.
}
\item{setnames}{
The list of set names to extract.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
This tool extracts the entries of the given setnames (stored in the header) from the given PeakML file and stores the result in the output file.
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
