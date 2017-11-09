\name{PeakML.DilutionTrendFilter}
\Rdversion{1.1}
\alias{PeakML.DilutionTrendFilter}
\title{PeakML.DilutionTrendFilter}
\description{
 This function will apply dilition trend filter to extract peak sets matching sample diltion trend as described in Jankevics et. al. 2012.
}
\usage{
PeakML.DilutionTrendFilter (filename,ionisation="detect",Rawpath=NULL,trendSets,p.value.thr=NULL,outputfile)
}
\arguments{
\item{filename}{
	Name ot the file which has to be processed in peakML format.
}
\item{ionisation}{
	Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
}
\item{trendSets}{
	Character vector containing a sample set names within input peakml file. Names of the sample sets should be given in the order matching dilution trend.
}
\item{p.value.thr}{
	Threhold for the p-value of correlation test to consider a both groups to be correlated. If set to NULL all peak sets are written to the output file with correltion and p.value annotations added to the peak groups.
}
\item{outputfile}{
	Name of the resulting peakML file.
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
	2 PeakML data files: file with peak sets matching given p-value threshold and peakML file with discarded peak sets. Annotationd dillution.corr (correlation value) and dillution.p.val (p-value of the correlation test) are added to the peak sets.
}
\references{
  Separating the wheat from the chaff: a prioritisation pipeline for the analysis of metabolomics datasets. 
Andris Jankevics, Maria Elena Merlo, Marcel de Vries, Roel J. Vonk, Eriko Takano, and Rainer Breitling 
Metabolomics 2012 8 (1), 29-36.
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis. Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
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
