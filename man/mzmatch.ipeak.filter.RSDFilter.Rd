\name{mzmatch.ipeak.filter.RSDFilter}
\alias{mzmatch.ipeak.filter.RSDFilter}
\title{mzmatch.ipeak.filter.RSDFilter}

\description{
Filters all the non-reproducibly detected signals from the given PeakML file.
}

\usage{
mzmatch.ipeak.filter.RSDFilter(JHeapSize=1425, i, o, rejected, rsd, h, v, sampleList, nSlaves, inputfolder="combined")
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). The only allowed format is PeakML and when it is not set the input is read from standard in.
}
\item{o}{
Option for the ouput file(s). The file is written in the PeakML file format and peaks that passed the defined filter are saved here. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs. When multiple input files are defined this option needs to point to a directory.
}
\item{rejected}{
Option for the file where to write all the rejected peaksets. When this option has not been set the rejected peaksets are not written.
}
\item{rsd}{
The RSD value as a percentage between [0..1]. This percentage indicates the maximum an individual measurement can deviate from the mean intensity of all the peaks in the set.
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
\item{inputfolder}{
  Character string containin a name of the folder where output files in peakML format will be created.
}

}

\details{
Filters all the non-reproducibly detected signals from the given PeakML file. The approach loads the PeakML file and determines based on the set-information, which of the sets of peaks contain signals that were not reproducibly detected in the replicates. The reproducability measure is defined as the Reproducibility Standard Deviation (RSD; stddev/mean), which is calculated from all the intensity values in a set. When the PeakML file contains multiple sets (in the header) the RSD is evaluated for each of the sets, and an entry is discarded when all of the sets fail the set threshold.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
1. Shah V, Midha K, Findlay J, et al: Bioanalytical Method Validation - A Revisit with a Decade of Progress. Pharmaceutical Research. 2000, 17:1551-1557

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
