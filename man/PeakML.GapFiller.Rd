\name{PeakML.GapFiller}
\alias{PeakML.GapFiller}
\title{PeakML.GapFiller}
\description{
	For each sample, identify peak sets where that sample is not represented. For each of those peak sets, extract an EIC from RAW data file and add it to the peak set.
}
\usage{
  PeakML.GapFiller(filename, ionisation = "detect", Rawpath = NULL, outputfile,ppm=0,rtwin=0,nSlaves=1, fillAll=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    The name of the PeakML file to load. This file needs to exist, otherwise the function will fail.
  }
  \item{ionisation}{
    If value is set to "detect" ionization mode will be detected from peakml file. Values "positive", "negative", "neutral" are also supported.
  }
  \item{Rawpath}{
    The path to the raw data, which is needed to read EIC's for peaks which were missed by peak picker. If value is set to NULL, filepath will be detected from peakml file. 
  }
  \item{outputfile}{
    The name of the PeakML file to write after gap filling.
  }
  \item{ppm}{
    Mass window in ppm which to use when mass chromatogram is extracted from RAW data file, before exporting to peakml file. Default value is 0.
  }
  \item{rtwin}{
    The retention time window in seconds, defining the range where to look for matches. Default values is 0.
  }
  \item{nSlaves}{
    The number of parralles instancet to be run simultaneously. (EXPERIMENTAL, USE AT YOUR OWN RISK)
  }
  \item{fillAll}{
	If set to TRUE all chromatograms (alse these which are present in peakml file) will be reintegrated from the raw data file with RT and mass window of the givern peak group (min and max RT and mass values are taken from all existing peaks.) This gives an output where all chromatograms have a constant RT window.+
  }
}
\details{

	After combining a biological/technical replicates in the single file, there will always be peak sets that do not include peaks from every sample. This method produces EICs for those missing samples by extracting data from raw data in peak set region.
}
\value{
  This function returns no value, peakml file with appended EICs is generated.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl), Richard Scheltema (r.a.scheltema@rug.nl), Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
}

\seealso{
  \code{\link{PeakML.xcms.read}}
}

\keyword{ mzmatch.R }
\keyword{ XCMS }
