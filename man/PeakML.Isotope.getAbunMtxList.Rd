\name{PeakML.Isotope.getAbunMtxList}
\Rdversion{1.1}
\alias{PeakML.Isotope.getAbunMtxList}
\title{PeakML.Isotope.getAbunMtxList}
\description{
	Generates a matrix containing the relative abundance of labelled vs unlabelled isotopes based on the the mean of the peaks in each groups.
}
\usage{
	PeakML.Isotope.getAbunMtxList(isotopeChroms, sampleGroups, useArea)
}

\arguments{
	\item{isotopeChroms}{
		The chrom data for isotopes see. PeakML.Isotope.getChromData.
	}
	\item{sampleGroups}{
		The sample groups.
	}
	\item{useArea}{
		TRUE will sum up all intensities to give area under the curve, FALSE gives maximum intensities.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getAbunMtxList returns a matrix containing the relative abundance of labelled vs unlabelled isotopes based on the the mean of the peaks in each groups, replicates along the columns and isotops along the rows.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Methods.getPhenoData} \link{PeakML.Methods.getChromData} \link{PeakML.Isotope.getIsotopes}}

}

\keyword{ mzmatch.R }
\keyword{ peakML }
