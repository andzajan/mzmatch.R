\name{mzmatch.ipeak.align.CowCoda}
\alias{mzmatch.ipeak.align.CowCoda}
\title{mzmatch.ipeak.align.CowCoda}

\description{
Performs retention time alignment on the given set of peaks.
}

\usage{
mzmatch.ipeak.align.CowCoda(JHeapSize=1425, i, o, ppm, codadw, order, maxrt, image, selection, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files, fow which the only allowed format is PeakML. Either a list of files (comma-separated) or a wild-card can be entered here.
}
\item{o}{
Option for the output directory, where the retention time aligned files are written.
}
\item{ppm}{
The precision of the data in parts-per-million. This accuracy value is used for matching the peaks from all of the measurement.
}
\item{codadw}{
The minimum mass chromatogram quality value a mass chromatogram needs have in order to qualify for the alignment process. This value is [0..1], where a higher value is better. The standard value is set to 0.8, which appears to give good results.
}
\item{order}{
The order of the polynomial fitted on the found deviations of the retention times. The higher this order is, the more flexible the polynomial will be, but less likely to give good results in low population areas.
}
\item{maxrt}{
The maximum retention time shift allowed in seconds.
}
\item{image}{
Option for a graph with all the used peaks and their deviation in seconds. The fitted polynomials are also displayed.
}
\item{selection}{
Option for the selection of peaks used to calculate the shifts in retention times.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Performs retention time alignment on the given set of peaks. The tool uses the CowCoDA algorithm, which looks for strong peaks in the set based on the CoDA algorithm. This particular approach uses the Durbin-Watson statistic for finding the best peaks. The best scoring peaks are then used to match over the different measurements and calculate the retention time shift. This data is consequently used to fit a polynomial function with the given order, which is used to align the data.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
Christin C, Smilde AK, Hoefsloot HC, Suits F, Bischoff R, Horvatovich PL. Optimized time alignment algorithm for LC-MS data: correlation optimized warping using component detection algorithm-selected mass chromatograms. 2008.

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
