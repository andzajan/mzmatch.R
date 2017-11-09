\name{mzmatch.ipeak.util.Identify}
\alias{mzmatch.ipeak.util.Identify}
\title{mzmatch.ipeak.util.Identify}

\description{
Identifies the contents of the given PeakML file with the given databases.
}

\usage{
mzmatch.ipeak.util.Identify(JHeapSize=1425, i, o, ppm, databases, minrt, maxrt, rtwindow, massOverride, polarity, adducts, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file. The only allowed file format is PeakML and no limitations are set to its contents. When this is not set the input is read from the standard in.
}
\item{o}{
Option for the output file. This file is writen in the same PeakML file format as the input file, with the addition of identification annotations (tag: 'identification' containing the database id's).
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for matching the masses to those found in the supplied databases. This value is obligitory.
}
\item{databases}{
Option for the molecule databases to match the contents of the input file to. These files should adhere to the compound-xml format.
}
\item{minrt}{
Optional minimum retention time for excluding signals from the input-file. This is for example convenient for excluding lipids on an LC/MS setup with a HILIC column.
}
\item{maxrt}{
Optional maximum retention time for excluding signals from the input-file. This is for example convenient for including lipids on an LC/MS setup with a HILIC column.
}
\item{rtwindow}{
Optional retention time window for finding matches from the databases. If this is set and the database contains entries with previously registered retention times for molecules. When this value is not set, the stored retention times are ignored.
}
\item{rtwindowrelative}{
If TRUE, the rtwindow parameter does not specify an absolute time, but instead specifies a fraction of the molecule time. E.g. rtwindow=0.05 rtwindowfraction=true specifies a 5 percent window around the database retention time.
}
\item{massOverride}{
Option to say whether to override the calculated formula mass with the monoisotopic mass element from the database. If this option is specified, andthe mass element is not present, then the calculated formula mass is used.
}
\item{polarity}{
Option to discard molecules of the opposite polarity. If this option is specified, any molecules of the opposite polarity will not be checked. Molecules with no polarity specified will always be checked
}
\item{adducts}{
Option to specify which adducts to search for. Default values are "M+H,M-H" and depending on ionisation mode of the data set, what will be determined automatically. 
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Identifies the contents of the given PeakML file with the given databases. The databases are expected to have the format of the example below (a standard file format for this would be preferable) and contain all of the compounds to be tested for. Within the tool-chain, in the package 'mzmatch.ipeak.db', several tools are provided for converting downloadable files (usually on an ftp-server) from the major metabolite databases. Before identifying your files, please make sure you have the most recent version.

Matching is performed only on mass, which is taken from the topmost structure. In other words, if the PeakML file contains a list of mass chromatograms the mass of each individual mass chromatogram is matched to the database. However, if the PeakML file contains a list of matched mass chromatograms, the mean mass of the matched mass chromatograms is used for identification.

When a peak is positively identified the annotation 'identification' is extended with the unique database ID corresponding to the match. This will keep the clutter in the PeakML file to a minimum and the information associated to the tag up-to-date when the database is updated with a new version. Additionally, the tag provides a convenient wayof removing false identifications from the PeakML file from a UI environment.
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
