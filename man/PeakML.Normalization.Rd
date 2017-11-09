\name{PeakML.Normalization}
\alias{PeakML.Normalization}
\title{PeakML.Normalization}
\description{
	General tool to perform basic normalisation on a PeakML file resulting from the Combine process (ie each entry is the combination of peaks from multiple measurements). By default normalization to Total Ion Current (TIC) for every sample is carried out. Function supports also input of predefined values for every sample in peakml file.
}
\usage{
  PeakML.Normalization (filename,ionisation="detect",Rawpath=NULL,outputfile,values=NULL)
}

\arguments{
  \item{filename}{
    The name of the PeakML file to load. This file needs to exist, otherwise the function will fail.
  }
  \item{ionisation}{
    If value is set to "detect" ionization mode will be detected from peakml file. Values "positive", "negative", "neutral" are also supported.
  }
  \item{Rawpath}{
    The path to the raw data, which is needed to calculate a TIC's from the raw data files. If value is set to NULL, filepath will be detected from peakml file. 
  }
  \item{outputfile}{
    The name of the PeakML file to write after normalization.
  }
  \item{values}{
    If set to NULL normalisation will be done to TIC values (calculated automatically). User can supply vector of the same length as the numbes of samples in peakml file with normalization values instead of TIC values.
  }
}
\details{

	Correct pre-treatment of metabolomic data is an essential step before the application of multivariate statistical methods. This function implements normalisation to the TIC values which can deal with technical variability of the machine. User can also supply dedicated vector of normalisation values as the extra argument of the function call.
}
\value{
  This function returns no value, peakml file with normalized data will be created.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl), Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.xcms.read}}
}

\keyword{ mzmatch.R }
\keyword{ XCMS }
