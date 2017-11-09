\name{mzmatch.ipeak.db.DBToText}
\alias{mzmatch.ipeak.db.DBToText}
\title{mzmatch.ipeak.db.DBToText}

\description{
Exports the given moleculedb database to an excel readable text file.
}

\usage{
mzmatch.ipeak.db.DBToText(JHeapSize=1425, i, o, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Input file of the database
}
\item{o}{
Output file for the tab-separated version of the database
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Exports the given moleculedb database to an excel readable text file. Only the general name is exported (ie the synonyms are not exported)
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
