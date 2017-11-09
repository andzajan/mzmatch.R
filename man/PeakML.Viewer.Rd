\name{PeakML.Viewer}
\Rdversion{1.1}
\alias{PeakML.Viewer}
\title{PeakML.Viewer}
\description{
  This function is a wrapper function for "Java" tool PeakML Viewer. It will automatically download a jar file from the Sourceforge project page designed for OS on which R system is running. And will start the application with all compound data base files from mzmathc.R package associated to it.
}
\usage{
  PeakML.Viewer (arch="detect",install.path=NULL, JHeapSize=1024, uninstall=FALSE)
}
\arguments{
\item{arch}{
  If set to "detect" installer will try to automatically detect an Operating system and machine arhitecture to download a PeakML Viewer version which can run on given machine. Supported valuse are: "Windows_32", "Windows_64", Linux_32", "Linux_64", "OSX_32" and "OSX_64".
}
\item{install.path}{
  If set to "NULL" will isntall programm to "mzmatch.R" installation folder. User can specify desrired installation folder as a character vector.
}
\item{JHeapSize}{
  Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.
}
\item{uninstall}{
  If set to "TRUE" command will remove existing PeakML Viewer files from your system. This command can be used if you want to update a software to the newer version.
}
}
\details{
  A user interface application called PeakML Viewer enables quick visualisation, inspection and manipulation of the contents of a PeakML file (e.g. manual selection and/or export of peaks of interest).

  After loading a PeakML file, the ``entry'' view gives an overview of all the entries with the retention time, mass and intensity. The ordering from the original file is kept intact, making the results from sorting tools like mzmatch.ipeak.sort.RelatedPeaks accessible. An entry is highlighted in bold when it has been matched to a compound from a database with mzmatch.ipeak.util.Identify (determined by verifying whether the entry contains the system-defined annotation identification). By clicking on an entry, the associated traces will be displayed in the "graph" view and the identifications in the "identification" view (including the ppm-deviation and the putatively assigned molecular structure when available). There is an additional tab derivatives, which shows all the clustered related peaks with their identification as determined by mzmatch.ipeak.sort.RelatedPeaks (stored in the system defined annotations relation.id and relation.ship). The "filter" view allows the user to perform simple operations on the data (sorting and filtering), for zooming in on the entries of interest. The "trend" view gives an overview of the intensity levels for each entry (which in the case of a peakset can for example consist of multiple mass chromatograms). All the peaks belonging to the same set are grouped together in this plot and the mean, minimum and maximum values are displayed. The "sets" view shows all the measurements used and how they are organized in sets. With the check boxes all the peaks from a set or peaks individually can be switched on or off (i.e. not displayed). The "annotations" view gives an overview of all the annotations that are available for the current entry.
}
\value{
  This function returns no value.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.xcms.write}}
}

\keyword{mzmatch.R}
\keyword{ XCMS }
