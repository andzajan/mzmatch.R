\name{PeakML.Isotope.plotSamples}
\Rdversion{1.1}
\alias{PeakML.Isotope.plotSamples}
\title{PeakML.Isotope.plotSamples}
\description{
	Plots the isotopes of the specified compound as a pdf document. The pdf document in which the plots has to be placed has to be initiated prior to calling this function (see PeakML.Isotope.processTargettedIsotopes)
}
\usage{
	PeakML.Isotope.plotSamples(isotopeChroms, metName, metFormula, metMass, metComment, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon, label, exclude_from_plots)
}

\arguments{
	\item{isotopeChroms}{
		The chromatograms of the isotope. See PeakML.Isotope.getChromData.
	}
	\item{metName}{
		Name of the metabolite.
	}
	\item{metFormula}{
		Formula of the metabolite.
	}
	\item{metMass}{
		Mass of the compound.
	}
	\item{metComment}{
		Any comment that has to be inserted along with the metabolite name in the plot.
	}
	\item{stdRT}{
		The standard retention time of the compound if known.
	}
	\item{sampleType}{
		The type of sample, eg., positive ionisation, negetive ionisation.
	}
	\item{sampleGroups}{
		The sample groups
	}
	\item{plotOrder}{
		The order in which the isotopes has to be plotted.
	}
	\item{useArea}{
		TRUE will sum up all intensities to give area under the curve, FALSE gives maximum intensities.
	}
	\item{followCarbon}{
		The number of the carbon in a compound whose labelled or unlabelled trend has to be plotted.
	}
	\item{label}{
		The type of labelling used in the data being analysed.
	}
    \item{exclude_from_plots}{
		Exclude the groups in this list from the plots in the output.
	}

}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.plotSamples plots the isotopes of the specifed compound as a pdf document.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getTrendList} \link{PeakML.Isotope.processTargettedIsotopes} \link{PeakML.Isotope.getIsotopes} \link{PeakML.Isotope.getRatioMtxList}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
