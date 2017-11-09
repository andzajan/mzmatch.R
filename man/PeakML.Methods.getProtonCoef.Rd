\name{PeakML.Methods.getProtonCoef}
\Rdversion{1.1}
\alias{PeakML.Methods.getProtonCoef}
\title{PeakML.Methods.getProtonCoef}
\description{

	Returns the proton coefficient based on the ionisation mode. If ionisation argument is set to "detect", then the ionisation mode is detected from the first scan in the first sample. No ionisation switching is supported in mzmatch.R yet.
}
\usage{
	PeakML.Methods.getProtonCoef(PeakMLtree, ionisation="detect")
}

\arguments{
	\item{PeakMLtree}{
		XML tree of PeakML file content.
	}
	\item{ionisation}{
		Defaults to "detect", to detect the ionisation from the scan data. Other accepted values are "positive" (proton coef = 1), "negative" (proton coef. = -1) and empty string (proton coef=0).
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
  getProtonCoef returns list of two elements: proton coefficient, a numeric value and ionisation mode as a character string ("positive","negative").
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.GapFiller}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
