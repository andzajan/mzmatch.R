\name{PeakML.Isotope.getTrendList}
\Rdversion{1.1}
\alias{PeakML.Isotope.getTrendList}
\title{PeakML.Isotope.getTrendList}
\description{
	Returns the trend of labelled and unlabelled carbons in the compound.
}
\usage{
	PeakML.Isotope.getTrendList(intList, sampleGroups, useArea)
}

\arguments{
	\item{intList}{
		List of intensities. See PeakML.Isotope.getChromData
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
	PeakML.Isotope.getTrendList returns the trend of labelled and unlabelled carbons in the compound as a list of mean of the sum or max of internsities of all replicates in the form trendList[[peakGroup]][[sampleGroup]][[isotop]] <- meanIntRep
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
