\name{mzmatch.ipeak.sort.RelatedPeaks}
\alias{mzmatch.ipeak.sort.RelatedPeaks}
\title{mzmatch.ipeak.sort.RelatedPeaks}

\description{
LC-MS experiments yield large amounts of peaks, many of which correspond to derivatives of peaks of interest (eg, isotope peaks, adducts, fragments, multiply charged molecules), termed here as related peaks.
}

\usage{
mzmatch.ipeak.sort.RelatedPeaks(JHeapSize=1425, i, o, basepeaks, ppm, rtwindow, minrt, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file. The only allowed format is PeakML and when it is not set the input is read from standard in. The contents of the file is enforced to be peaksets (result of Combine) as this tool utilizes the full information of peaksets in order to identify related peaks.
}
\item{o}{
Optional filename where the output is written. If this is not set the output is written to the standard output.
}
\item{basepeaks}{
Optional filename where the output is written. If this is not set the file with the basepeaks is not written.
}
\item{ppm}{
The mass tolerance in ppm. This value is used for the identification of the relations in the found sets of related peaks.
}
\item{rtwindow}{
The retention time window in seconds, defining the range where to look for matches.
}
\item{minrt}{
Denotes the minimum retention time in seconds peaks should occur before being taken into account for the relation process. When this value is not set all peaks are taken into account.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
LC-MS experiments yield large amounts of peaks, many of which correspond to derivatives of peaks of interest (eg, isotope peaks, adducts, fragments, multiply charged molecules), termed here as related peaks. This tool clusters groups of related peaks together and attempts to identify their relationship. Each cluster is given a unique id, which is stored as an annotation with the name 'relation.id'. The relationship is stored as an annotation with the name 'relation.ship'.

The possibility is offered only to store the most intense peak of each cluster (option 'basepeak'). This is useful for cleaning up the data before attempting statistical approaches to mine the data. However, as the peak of interest is not per definition the most intense peak one cannot rely on this file for identification purposes.

It is advised to only run this application on a combined file containing a complex experiment (e.g. time-series, comparison, etc), as the performance of the approach is expected to improve with increasingly complex setups. The approach iterates through the peak list sorted on descending intensity and works threefold for identifying all of the derivatives. (1) the most intense, not yet processed peak is selected and all not yet processed peaks within the specified retention time window are collected. (2) First all of the collected peaks are correlated to the most intense peak on signal shape. For this only the peaks from the same measurement are compared, as otherwise distorted peak shapes would unfairly reduce the score. (3) The surviving peaksets are then correlated on intensity trend (hence the need for complex setups). This is the most informative feature, as it is not likely that that two co-eluting peaks have the same intensity trend for complex setups. All signals to correlate well enough to the selected, most intense signal are then clustered with this most intense signal.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
1. Tautenhahn R, Bottcher C, Neumann S: Annotation of LC/ESI-MS Mass Signals. In: Bioinformatics Research and Development. Hochreiter S, Wagner R (Ed.), Springer-Verlag, Heidelberg, Germany, 371-380 (2007).
2. Scheltema RA, Decuypere S, Dujardin JC, Watson DG, Jansen RC, and Breitling R: A simple data reduction method for high resolution LC/MS data in metabolomics. Bioanalysis - in press.

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
