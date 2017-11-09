\name{mzmatch.ipeak.Combine}
\alias{mzmatch.ipeak.Combine}
\title{mzmatch.ipeak.Combine}

\description{
Combines the contents of a set of PeakML files, containing either mass chromatograms or backgroundions at the lowest level (the signal of which is needed in order to make a correct assesment of similarness).
}

\usage{
mzmatch.ipeak.Combine(JHeapSize=1425, i, o, label, labels, ppm, rtwindow, combination, h, v, sampleList=NULL, nSlaves=1, outputfolder="combined")
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files. Multiple files can be passed by separating them with a comma (ie ,) or the use of a name with a wildcard (eg samples_*hrs.xml). The only allowed file format is PeakML containing either mass chromatograms or backgroundions at the lowest level (ie the result of another Combine can be used).
}
\item{o}{
Option for the ouput file. The resulting matches are written to this file in the PeakML file format.
When this option has not been set the output is written to the standard output. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{label}{
Optional label for the set being made. The labels are stored in the header of the resulting file and used for display purposes.
}
\item{labels}{
Optional labels for the input files. When these are used make sure to give as many labels as there are input files. The labels are stored in the header of the resulting file and used for display purposes.
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for the matching of mass chromatogram (collections) and needs to bereasonable for the equipment used to make the measurement (the LTQ-Orbitrap manages approximately 3 ppm).
}
\item{rtwindow}{
The retention time window in seconds, defining the range where to look for matches.
}
\item{combination}{
- set
  The files are to be combined as a true set.
- technical
  The files are to be combined as technical replicates.
- biological
  The files are to be combined as biological replicates.

}
\item{h}{
When this is set, the help is shown.
}

\item{v}{
When this is set, the progress is shown on the standard output.
}

\item{sampleList}{
Output objetc created by "mzmatch.R.Setup" function, it is used to assign sample names classes.
}

\item{nSlaves}{
number of slaves/cores to be used for parallel peak detection. Snow package is employed for multicore support.
}
\item{outputfolder}{
  Character string containin a name of the folder where output files in peakML format will be created.
}

}

\details{
Combines the contents of a set of PeakML files, containing either mass chromatograms or backgroundions at the lowest level (the signal of which is needed in order to make a correct assesment of similarness). The approach starts from the most intense, unprocessed peak in the complete set signals, covering all the measurements, and attempts to find all those signals from the other measurements falling within the mass window (option 'ppm') and the retention time window (option 'rtwindow'). The correct match from each measurement to the currently most intense signal is consequently identified by optimizing on the difference in area under the curve. Signals caused by the same analyte are expected to roughly have a similar shape and retention time. All matched signals are then clustered for the output and marked as processed, after which a new iteration is begun.

This tool can be used to iteratively compound files. For example, when analyzing a time series experiment, firstly the biological replicates of each timepoint can first be combined into a single file, the set of which can be labeled as biological replicates with the option 'combination'. Additional filtering operations can be applied in order to the timepoint combinations before proceeding to combine all the timepoints in a final set.
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
