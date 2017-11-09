\name{mzmatch.ipeak.convert.ConvertToText}
\alias{mzmatch.ipeak.convert.ConvertToText}
\title{mzmatch.ipeak.convert.ConvertToText}

\description{
Converts the contents of a PeakML file to a tab-separated text file.
}

\usage{
mzmatch.ipeak.convert.ConvertToText(JHeapSize=1425, i, o, databases, annotations, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
The input file in the PeakML file format. When this option is not set, the input is read from the standard input.
}
\item{o}{
The output file in text format. When this option is not set, the ouput is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{databases}{
Optional list of database-files in the moleculedb format. These databases are used to translate the id's in the annotation 'identifications' to full molecule names.
}
\item{annotations}{
Optional list of annotations to be exported.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Converts the contents of a PeakML file to a tab-separated text file. In order to make the file useable in a spreadsheet environment, only the top-level data for each measurement is given in the file. This application has been extended with support for outputting annotations.
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
