\name{PeakML.Isotope.DB2Text}
\Rdversion{1.1}
\alias{PeakML.Isotope.DB2Text}
\title{PeakML.Isotope.DB2Text}
\description{
	Generates a R data frame containing masses, RTs, formulas, ppm, identification, name, classes and subclasses for all peaks in the PeakML file that were identified based on the database specified.
}
\usage{
	PeakML.Isotope.DB2Text(filename, databases)
}

\arguments{
	\item{filename}{
		The name of the PeakML file which contains the combined peaks of all samples that are filtered, gap filled and identified. See the advanced metabolomics data processing tutorial at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.html for more information.
	}
	\item{databases}{
		A vector containing databases (in numbers) based on which initial annotations were made. For more information, see 'Identify peaks from data bases' section at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.php
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.DB2Text returns a data frame containing masses, RTs, formulas, ppm, identification, name, classes and subclasses for all peaks in the PeakML file that were identified based on the database specified.
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
