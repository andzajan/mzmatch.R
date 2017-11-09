\name{mzmatch.ipeak.filter.SimpleFilter}
\alias{mzmatch.ipeak.filter.SimpleFilter}
\title{mzmatch.ipeak.filter.SimpleFilter}

\description{
Simple filter tool for PeakML files.
}

\usage{
mzmatch.ipeak.filter.SimpleFilter(JHeapSize=1425, i, o, rejected, databases, ppm, n, offset, mindetections, minscanid, maxscanid, minretentiontime, maxretentiontime, minmass, maxmass, minintensity, maxintensity, annotations, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). The only allowed format is PeakML and when it is not set the input is read from standard in. When multiple input files are set the output needs to be set to a directory.
}
\item{o}{
Option for the ouput file. The file is written in the PeakML file format and peaks that passed the defined filter are saved here. When this option is not set the output is written to the standard out. Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{rejected}{
Option for the reject file. The file is written in the PeakML file format and peaks that have not passed the defined filter are saved here.
}
\item{databases}{
Optional file containing molecules to filter on. The file should be in a tab-delimited format with name 	 formula. The mass is automatically calculated from the formula and used, together with the ppm option, to find all mass chromatograms within the mass window.
}
\item{ppm}{
Optional ppm-value, which is used in conjunction with the molecule-file. For each molecule the ppm-value is used to determine the mass window in which to look for mass chromatograms matching the molecule.
}
\item{n}{
Optional value indicating the maximum number peaks to filter from the file. When this value is not set all peaks are taken into acount.
}
\item{offset}{
Optional value for skipping the first number of peaks. When this value is not set the application starts at the beginning.
}
\item{mindetections}{
Works only for peak-sets; checks whether a set contains the minimum number of peaks.
}
\item{minscanid}{
Optional value for taking into account only those peaks above the given minimum scanid. When this value has not been set all peaks are taken into account.
}
\item{maxscanid}{
Optional value for taking into account only those peaks below the given maximum scanid. When this value has not been set all peaks are taken into account.
}
\item{minretentiontime}{
Optional value for taking into account only those peaks above the given minimum retentiontime. When this value has not been set all peaks are taken into account.
}
\item{maxretentiontime}{
Optional value for taking into account only those peaks below the given maximum retentiontime. When this value has not been set all peaks are taken into account.
}
\item{minmass}{
Optional value for taking into account only those peaks above the given minimum mass. When this value has not been set all peaks are taken into account.
}
\item{maxmass}{
Optional value for taking into account only those peaks below the given maximum mass. When this value has not been set all peaks are taken into account.
}
\item{minintensity}{
Optional value for taking into account only those peaks above the given minimum intensity. When this value has not been set all peaks are taken into account.
}
\item{maxintensity}{
Optional value for taking into account only those peaks below the given maximum intensity. When this value has not been set all peaks are taken into account.
}
\item{annotations}{
Optional value for taking into account only those peaks containing the given annotations.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Simple filter tool for PeakML files. The filter works on the top-most data of an entry (in the case of a peak-set the compounded information of that set is used) and can be used to filter out entries not conforming to the specified filter(s). As the filter only works on top level data, the easiest filters are based on the top level properties: mass, intensity, retention time and scanid. For each of these properties the appropriate filters have been defined. Additionally, as we have access to mass information, the tool also allows the user to load databases (in moleculedb format) and set a ppm-value to filter out those entries that match to a molecule in the databases.

For peaksets the filter for minimum number of detections removes all peaksets with a number of peaks below the set threshold.

Please note that the filter-commands are stacked when multiple are given. For example, when setting a minimum intensity and mindetections are both set, all entries not qualifying both criteria are removed.
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
