\name{PeakML.bubbleplot}
\alias{PeakML.bubbleplot}
\title{PeakML.bubbleplot}
\description{
This function generates "bubble" plot of XC/MS experimental data. This graphical interpretation is suitable for quick assesing of data quality.
}
\usage{
PeakML.bubbleplot(xcmsset, peakidentifications, samples = NULL, samples2 = NULL, DB1string = NULL, bubbletype = "filled", plottype = "identified", outputfile = "bubbleplot.pdf", cor.threshold = 0.7, rtrange = NULL, mzrange = NULL, intensity.threshold.min = NULL, intensity.threshold.max = NULL, strict.relation.select = TRUE, scale = "max", logbase = 10, integral.val = "maxint", bubblesize = 3, smalllegend=FALSE)
}

\arguments{
  \item{xcmsset}{
The XCMS set structure containing the data to be written to the PeakML file.
}
  \item{peakidentifications}{
Matrix of peak attributes extracted by mzmatch software and PeakML.xcms.read function.
}
  \item{samples}{
Which samples should be plotted.
}
  \item{samples2}{
If 'plottype' is "correlation" sample name to correlate with should be given.
}
  \item{DB1string}{
Name of database file which should be used main data base for coloring
}
  \item{bubbletype}{
Values: "open" or "closed". 
}
  \item{plottype}{
"identified","compare" or "correlation"
}
  \item{outputfile}{
Name of output file in which graphical output shoud be generated (pdf format).
}
  \item{cor.threshold}{
If 'plottype' is 'correlation', threshold for correlated peaks can be set with this value.
}
  \item{rtrange}{
Numeric vector with two elements indicating RT range at which bublle plot should be plotted.
}
  \item{mzrange}{
Numeric vector with two elements indicating m/z range at which bublle plot should be plotted.
}
  \item{intensity.threshold.min}{
Minimum intensity threshold cut-off value.
}
  \item{intensity.threshold.max}{
Maximum intensity threshold cut-off value.
}
  \item{strict.relation.select}{
Restric relation classes for related peaks matching.
}
  \item{scale}{
Set scaling which to use before plotting.
}
  \item{logbase}{
If logaritmical scaling is used, log base can be defined.
}
  \item{integral.val}{
Which value to use for circle size, max intensity - 'maxo', or peak integral - 'into'
}
  \item{bubblesize}{
Scaling factor for bubble size, it depends on scaling which is used.
}
  \item{smalllegend}{
Defaults to FALSE. If is set to TRUE, legend is plotted only for unidentified,identified and their related features.
}
}

\value{
Graphical output in pdf file.
}
\references{
PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
Andris Jankevics (a.jankevics@rug.nl)
Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
\code{\link{PeakML.xcms.read}}
}

\keyword{ mzmatch.R }
\keyword{ XCMS }
\keyword{ peakML }
