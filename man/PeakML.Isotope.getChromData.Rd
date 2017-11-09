\name{PeakML.Isotope.getChromData}
\Rdversion{1.1}
\alias{PeakML.Isotope.getChromData}
\title{PeakML.Isotope.getChromData}
\description{
	Generates a list of lists containing mass, intensities and retentiontimes of all the isotopes in the isotopeList.
}
\usage{
	PeakML.Isotope.getChromData(isotopeList, chromDataList, phenoData, sampleGroups)
}

\arguments{
	\item{isotopeList}{
		List of isotopes see PeakML.Isotope.getIsotopes.
	}
	\item{chromDataList}{
		The chrom data list see PeakML.Methods.getChromData.
	}
	\item{phenoData}{
		See PeakML.Methods.getPhenoData.
	}
	\item{sampleGroups}{
		The sample groups.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.getChromData returns a list of lists containing mass, intensities and retentiontimes of all the isotopes in the isotopeList. Each sublist will be of the form rtList[[peakGroup]][[sampleGroup]][[isotop]][[replicate]] where the signals/values will be stored.
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
