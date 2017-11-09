\name{mzmatch.ipeak.normalisation.Normalisation}
\alias{mzmatch.ipeak.normalisation.Normalisation}
\title{mzmatch.ipeak.normalisation.Normalisation}

\description{
General tool to perform basic normalisation on a PeakML file resulting from the Combine process (ie each entry is the combination of peaks from multiple measurements).
}

\usage{
mzmatch.ipeak.normalisation.Normalisation(JHeapSize=1425, i, o, type, v, h)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file. The only allowed format is PeakML and when it is not set the input is read from standard in. The tool expects a combination of peaks from different sets and will exit when this is not encountered.
}
\item{o}{
Option for the ouput file. The file is written in the PeakML file format and peaks that passed the defined filter are saved here. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{type}{
- znormalization
  Z-score normalization, scales the intensity values to 0 mean and 1 stddev.- quantile  Quantile normalisation, essentially rank-based normalisation.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}
\item{h}{
When this is set, the help is shown.
}

}

\details{
General tool to perform basic normalisation on a PeakML file resulting from the Combine process (ie each entry is the combination of peaks from multiple measurements). The commonly used z-score and quantile normalisation procedures have been implemented, which adjust the intensity-values according to the factors they find.
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
