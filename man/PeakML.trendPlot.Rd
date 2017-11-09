\name{PeakML.trendPlot}
\Rdversion{1.1}
\alias{PeakML.trendPlot}
\title{PeakML.trendPlot}
\description{
  A trend plot graphs the samples (probably in time series experiments) against a variable of interest (intensity, measured mass). It is used to examine whether or not there is a relationship between the signal examined and time.
}
\usage{
  PeakML.trendPlot(filename,sampleclasses,samplenames=NULL,outputfile,xaxisvalues=NULL,xaxislabels=NULL,log2scale=FALSE,collapse=FALSE,plotchromatograms=FALSE,xaxislabel=NULL,percentiles=c(15,20,35),classaverage=FALSE,legendposition="bottomright",yminvalue=NULL,xmaxvalue=NULL,colors= c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"),fontfamily="Helvetica",pdflayout=c(3,1))  
}

\arguments{
\item{filename}{
	The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function will fail.
}
\item{sampleclasses}{
	Character vector of length equal to count of measurements in PeakML data file. Every string in this vector should describe sample class.
}
\item{samplenames}{
	Defaults to NULL. If the order of samples should be changed before plotting measurement names in expected order should be provided. Order of class values also should be ordered according names "samplenames" vector.
}
\item{outputfile}{
	Name of pdf file which will be generated in current working directory of R.
}
\item{xaxisvalues}{
	Numeric vector of same length as number of measurements in PeakMl data file. If is set to 'NULL' sample index values are plotted on x axis.
}
\item{xaxislabels}{
	Character vector of same length as number of measurements in PeakMl data file. If is set to 'NULL' sample index values are plotted on x axis. This parameter can be supplied only if xaxisvalues are supplied. Xaxis will be labeled by character values instead of numeric ones in xaxisvalues vector.
}
\item{log2scale}{
	Applies log2 transformation on intensity values before plotting. Logical: FALSE or TRUE
}
\item{collapse}{
	Logical. If is set to TRUE, overlay plot by sample classes will be generated.
}
\item{plotchromatograms}{
	Logical. If is set to TRUE, pdf file with prefix "masschromatograms" containing mass traces of every peak will be generated.
}
\item{xaxislabel}{
	Character string which will be used to label x axis.
}
\item{percentiles}{
	Vector of the numeric values indicating percentage thresholds at which deviation values from mean intensity value should be plotted. If set to NULL, no percentiles are plotted. 
}
\item{classaverage}{
	Average intensity values for data points in each sample class will be calculated and plotted by sample class. Error bars indicates standart deviations from mean value. 
}
\item{legendposition}{
	At which side plot legend should be added. Possible values are character strings: "bottomleft", "bottomright", "topleft", "topright". 
}
\item{yminvalue}{
	At which value yaxis should start. If set to NULL, axis ranges are calculated automatically.
}
\item{xmaxvalue}{
	Maximum value for xaxis can be set, this is helpfull for placing legend on the plot. If set to NULL, axis ranges are calculated automatically.
}
\item{colors}{
	You can improve look of graphics by defining color palette which to use for plotting by using 657 pre-defined colors are available. Value should be character vector of the same length as the sampleclasses defined. Available values for color definitions can be obtained by R command: colors (). Graphical representation of R colors are also shown on this web-page: http://research.stowers-institute.org/efg/R/Color/Chart/
	Numeric values of colors can also be supplied instead of character values. 
	Default color palettes is: black, red, green3, blue, cyan, magenta, yellow, gray. That means that points for 8 different sample groups can be coloured in different color by default. 
}

\item{fontfamily}{
	The initial font family to be used, normally as a character string. In normal use it is one of "AvantGarde", "Bookman", "Courier", "Helvetica", "Helvetica-Narrow", "NewCenturySchoolbook", "Palatino" or "Times", and refers to the standard Adobe PostScript fonts families of those names which are included (or cloned) in all common PostScript devices. Refer to help (postscript) for more detailed description of available font families and specific characters. Default value: "Helvetica"

}
\item{pdflayout}{
	Numeric vector with two elements. The first number indicates in how many rows to place plots on every page of pdf document. The second number defines number of columns in pdf document.

}

}

\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.trendPlot function allows to monitor intensity trends between samples in PeakML data file. Sample order before plotting can be changed by using 'samplenames' argument. Average values plot (between sample classes) with error bars can be plotted if parameter "classaverage" is set to TRUE. 
Y axis range is calculated automatically, but minimum value can be override by 'yminvalue' parameter.

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
\keyword{ peakML }
