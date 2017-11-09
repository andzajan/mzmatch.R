\name{PeakML.Methods.DBidToCompoundName}
\Rdversion{1.1}
\alias{PeakML.Methods.DBidToCompoundName}
\title{PeakML.Methods.DBidToCompoundName}
\description{
 This function will convert identification ID's from data bases to actual compound names.
}
\usage{
PeakML.Methods.DBidToCompoundName (DBS,PeakMLdata,collapse=TRUE)
}
\arguments{
\item{DBS}{
	List of comma seperated data bases in XML format to be used.
}
\item{PeakMLdata}{
	R obejct in PeakMLdata format, normaly output from the PeakML.Read function.
}
\item{collapse}{
	If set to TRUE ID's matching several data bases will coconated as a comma separated text string.
}
}


\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.Methods.getGroupAnnotations function allows to access attributes defined in PeakML file.
}
\value{
	Matrix showing resolved ID's, compound formual and DB source.
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.Read}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
