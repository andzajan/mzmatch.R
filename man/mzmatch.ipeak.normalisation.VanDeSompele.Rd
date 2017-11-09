\name{mzmatch.ipeak.normalisation.VanDeSompele}
\alias{mzmatch.ipeak.normalisation.VanDeSompele}
\title{mzmatch.ipeak.normalisation.VanDeSompele}

\description{
Applies a basic normalisation scheme to the PeakML file resulting from the combine process.
}

\usage{
mzmatch.ipeak.normalisation.VanDeSompele(JHeapSize=1425, i, o, selection, selection_normalized, img, factors, ppm, database, v, h)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file. The only allowed format is PeakML and when it is not set the input is read from standard in. The tool expects a combination of peaks from different sets and will exit when this is not encountered.
}
\item{o}{
Option for the ouput file. The file is written in the PeakML file format and contains all the peaks with normalized intensities. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{selection}{
Option for the file where the un-normalized selection of peaks is written.
}
\item{selection_normalized}{
Option for the file where the normalizaed selection of peaks is written.
}
\item{img}{
Option for the file where a graph of the normalization factors is written. This file is in PDF format.
}
\item{factors}{
Option for the file where the normalization factors are written.
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for matching the masses to those found in the supplied databases. This value is obligitory.
}
\item{database}{
Option for the molecule databases to match the contents of the input file to. This file should adhere to the compound-xml format.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}
\item{h}{
When this is set, the help is shown.
}

}

\details{
Applies a basic normalisation scheme to the PeakML file resulting from the combine process. The method looks for house-hold metabolites not expected to change. These metabolites are defined as the metabolites showing the least changes compared to the detections in all of the other measurements. This is reflected by a stability score (stored as an annotation for each of the entries labeled 'stability factor'). These stability scores are calculated from all the metabolites that could be identified with the given database, within the given ppm-range. The identified metabolites need to have been identified in all of the used measurements.

From the top 10% (at least 10) most stable, identifed metabolites the normalisation factors are calculated and applied to the rest of the data.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
1. Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol. 2002.

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
