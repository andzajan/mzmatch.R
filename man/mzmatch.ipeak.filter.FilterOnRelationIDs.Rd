\name{mzmatch.ipeak.filter.FilterOnRelationIDs}
\alias{mzmatch.ipeak.filter.FilterOnRelationIDs}
\title{mzmatch.ipeak.filter.FilterOnRelationIDs}

\description{
Retrieves all entries annotated with one of the given list of relation-ids.
}

\usage{
mzmatch.ipeak.filter.FilterOnRelationIDs(JHeapSize=1425, i, relationids, o, rejected, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). The only allowed format is PeakML and when it is not set the input is read from standard in.
}
\item{relationids}{
Either a text-file with the relation-id's (1 id per line) or a PeakML file with entries containing the id's to filter on.
}
\item{o}{
Option for the ouput file(s). The file is written in the PeakML file format and peaks that passed the defined filter are saved here. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs. When multiple input files are defined this option needs to point to a directory.
}
\item{rejected}{
Option for the reject file. The file is written in the PeakML file format and peaks that have not passed the defined filter are saved here.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Retrieves all entries annotated with one of the given list of relation-ids. This list can either be specified with a text-file, where each line contains a single relation-id, or a PeakML file, where each entry containing a relation-id annotation is used to build up the relation-id list.

This tool is very practical when manualy selecting strong peaks from the extracted set of peaks. By running 'mzmatch.ipeak.sort.RelatedPeaks' two files can be generated. The result file containing all the entries, and the basepeaks file containing the most intense entry of each cluster of related derivatives. The basepeaks file is a drastical reduction of the complete file and can more easily be manually filtered with for example the PeakML Viewer. The resulting file can then be used as input to this tool, in order to retrieve the complete clusters of interest.
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
