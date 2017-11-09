\name{PeakML.Isotope.filteroutUnlabelled}
\Rdversion{1.1}
\alias{PeakML.Isotope.filteroutUnlabelled}
\title{PeakML.Isotope.filteroutUnlabelled}
\description{
	Used during untargetted isotope profiling. Filters out and removes all peak-sets that do not have sufficient (determined by the stringency) number of labelled peaks corresponding to the samples involved.
}
\usage{
	PeakML.Isotope.filteroutUnlabelled(isotopeList, numCarbons, fillGaps, sampleNames, stringency=30)
}

\arguments{
	\item{isotopeList}{
		The list of isotopes returned from PeakML.Isotope.getIsotopes.R function.
	}
	\item{numCarbons}{
		The number of carbons in the compound of interest. See PeakML.Methods.getCarban.R.
	}
	\item{fillGaps}{
		"ALLPEAKS" if all the peaks are to gap filled, requied to fill missing signals in some part of the curve. "MISSINGPEAKS", if only the missing peaks have to be filled Or "NONE" if no peaks have to be gap filled. Defaults to "ALLPEAKS".
	}
	\item{sampleNames}{
		The names of samples.
	}
	\item{stringency}{
		Percentage of labelled peaks that are expected to be present with respect to the total number of samples.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.filteroutUnlabelled filters out and removes all peak-sets that do not have sufficient (determined by the stringency) number of labelled peaks corresponding to the samples involved. Output is a list of isotopes minus those that do not pass the stringency criteria.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Methods.getPhenoData} \link{PeakML.Isotope.processUntargettedIsotopes}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
