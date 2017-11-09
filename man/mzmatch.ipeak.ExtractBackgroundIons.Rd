\name{mzmatch.ipeak.ExtractBackgroundIons}
\alias{mzmatch.ipeak.ExtractBackgroundIons}
\title{mzmatch.ipeak.ExtractBackgroundIons}

\description{
Extracts background ions (x-axis: RT; y-axis: Intensity) from 2D mass spectrometry data (LC/MS or GC/MS).
}

\usage{
mzmatch.ipeak.ExtractBackgroundIons(JHeapSize=1425, i, o, ppm, threshold, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files, which should be in one of the open standard file formats (mzML, mzXML or mzData) and contain data from a 2D mass spectrometry setup (LC/MS or GC/MS).
When this option has not been set, the input is read from the stdin (allowing for pipeline building). When a single input file is defined, the output '-o' should contain the output filename. When multiple input files are defined, the output '-o' should define an output directory.
For now only centroid input data is supported.
}
\item{o}{
Option for the ouput file(s); refer to the option input '-i' for a description of behaviours with regards to multiple input files. The extracted backgroundions are written here in the PeakML format.
When this option has not been set the output is written to the standard output (works only when there is a single input file).Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for the collection of the data-points belonging to a background ion and needs to be reasonable for the equipment used to make the measurement (the LTQ-Orbitrap manages approximatetly 3 ppm).
}
\item{threshold}{
The percentage (a value between 0 and 1) of scans that minimally need to contain a measurement.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Extracts background ions (x-axis: RT; y-axis: Intensity) from 2D mass spectrometry data (LC/MS or GC/MS). The raw data is loaded from the open standard file formats (mzML, mzXML or mzData) and all of the individidual mass traces (M/Z +/- ppm over the whole scan range) are retrieved. A mass trace is retained when it is present in more than the given percentage of scans (option 'threshold').

The resulting output file is in PeakML-format, containing a list of all the extracted background ions.

Remarks
1. At this time only centroid data is supported.
2. NetCDF is not supported as it misses necessary meta-information
3. Direct injection data will not yield correct results

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
