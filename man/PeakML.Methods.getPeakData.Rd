\name{PeakML.Methods.getPeakData}
\Rdversion{1.1}
\alias{PeakML.Methods.getPeakData}
\title{PeakML.Methods.getPeakData}
\description{

	Returns peak data from mass chromatograms generated from PeakML files as a matrix of 11 columns. Each column represents average mz, min mz, max mz, averate RT (this is recalculated at maximum intensity of peak, warning -this approach differs from that one which is used in XCMS), min RT,  max RT, sum of intensity, max intensity, measurement id,  group id and peak count repectively.
	
}
\usage{
	PeakML.Methods.getPeakData(PeakMLtree, chromDataList)
}

\arguments{
	\item{PeakMLtree}{
		XML tree structure of the PeakML file.
	}
	\item{chromDataList}{
		List of mass chromatogramms as an output of PeakML.Methods.getChromData
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
  A numeric matrix containing 11 columns representing - "AVGMZ","MINMZ","MAXMZ","AVGRT","MINRT",  "MAXRT", "SUMINTENSITY", "MAXINTENSITY", "measurement id", "group id", "peakcount" ,
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getMassCorrection}, \link{PeakML.Methods.getProtonCoef}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
