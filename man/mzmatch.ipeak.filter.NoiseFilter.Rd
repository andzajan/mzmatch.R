\name{mzmatch.ipeak.filter.NoiseFilter}
\alias{mzmatch.ipeak.filter.NoiseFilter}
\title{mzmatch.ipeak.filter.NoiseFilter}

\description{
Filters noise from PeakML files, containing mass chromatograms at the lowest level.
}

\usage{
mzmatch.ipeak.filter.NoiseFilter(JHeapSize=1425, i, o, rejected, codadw, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). The only allowed format is PeakML and when it is not set the input is read from standard in.
}
\item{o}{
Option for the ouput file. The file is written in the PeakML file format and peaks that passed the noise filter are saved here. When this option is not set the output is written to the standard out.
Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{rejected}{
Option for the file where to write all the rejected peaksets. When this option has not been set the rejected peaksets are not written.
}
\item{codadw}{
With this option the CoDA-DW filter is activated. The CoDA algorithm calculates a quality value for mass chromatograms (MCQ value [0..1] where higher is better) with the Durbin-Watson statistic, which is used to remove invalid peaks. The double value defines the hard threshold on which the mass chromatograms are removed from the set.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Filters noise from PeakML files, containing mass chromatograms at the lowest level. When the file contains a list of sets of mass chromatograms, the maximum score for the used method is calculated and compared to the given threshold. This is the best approach, as we expect that high quality can match up to low quality signals at the same mass and retention time. Only those entries scoring above the given threshold are stored in the output file. The rejected can be stored in a separate file (option 'rejected') for inspection or recovery.

The option 'codadw' can be used to set the threshold for the CoDA Durbin-Watson noise filtering approach. Normally the Durbin-Watson criterion results in a value between 0 and 4, where higher means a large amount of periodicity in the signal and lower vica versa. For mass chromatograms we expect little periodicity in the signal, thus a lower value is preferable. However, in order to preserve unity in our quality scores the CoDA-DW score is scaled between 0..1, where higher is better (less periodicity in the signal). As a general rule-of-thumb, for high quality mass chromatograms a score >0.8 is expected.

Remarks:
- CoDA-DW is scaled between 0..1, where higher is better mass chromatogram quality.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
1. Windig W: The use of the Durbin-Watson criterion for noise and background reduction of complex liquid chromatography/mass spectrometry data and a new algorithm to determine sample differences. Chemometrics and Intelligent Laboratory Systems. 2005, 77:206-214.

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
