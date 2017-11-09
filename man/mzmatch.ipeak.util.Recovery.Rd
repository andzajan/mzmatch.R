\name{mzmatch.ipeak.util.Recovery}
\alias{mzmatch.ipeak.util.Recovery}
\title{mzmatch.ipeak.util.Recovery}

\description{
Complements missing signals to a combination file produced by mzmatch.
}

\usage{
mzmatch.ipeak.util.Recovery(JHeapSize=1425, i, o, recovered, bins, ppm, rtwindow, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files. The only allowed file format is PeakML containing either mass chromatograms or backgroundions at the lowest level (ie the result of another Combine can be used).
}
\item{o}{
Option for the ouput file. The resulting matches are written to this file in the PeakML file format.
When this option has not been set the output is written to the standard output. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{recovered}{
Optional output file where the all the sets that have been complemented with missing signals are written.
}
\item{bins}{
Option for the files containing the missing signals. The only allowed file format is PeakML containing either mass chromatograms or backgroundions at the lowest level (ie the result of another Combine can be used).
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for the matching of mass chromatogram (collections) and needs to bereasonable for the equipment used to make the measurement (the LTQ-Orbitrap manages approximately 3 ppm).
}
\item{rtwindow}{
The retention time window in seconds, defining the range where to look for matches.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Complements missing signals to a combination file produced by mzmatch.ipeak.Combine (or another combine tool for that matter), containing sets of signals (mass chromatograms or backgroundions) produced by the same analyte in different measurements. It is possible that the combiner incorrectly missed signals, or was offered files lacking the missing signals (eg signals from a timepoint removed by the tool mzmatch.filter.RSDFilter, which can be linked to a stable signal from another timepoint).

For this application to function correctly it is required that the files containing the missing signals are generated from the same measurement, as the id's for the measurement needs to correspond.
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
