\name{PeakML.Methods.getRawMat}
\Rdversion{1.1}
\alias{PeakML.Methods.getRawMat}
\title{PeakML.Methods.getRawMat}
\description{
	This function returns a matrix of peaks from raw LC/MS data file wihing given MZ an RT window.
}
\usage{
	PeakML.Methods.getRawMat(rawfile, scan_start, scan_finis, mz_start, mz_finis, correctedRT, uncorrectedRT)
}

\arguments{
	\item{rawfile}{
		R object of the raw data file created from openMSfile () command (mz.R package).
	}
	\item{scan_start}{
		Numeric value specifying starting scan from which to read data. 
	}
	\item{scan_finis}{
		Numeric value specifying end scan from which to read data. 
	}
	\item{mz_start}{
		Lowest m/z range to be returned.
	}
	\item{mz_finis}{
		Highest m/z range to be returned.
	}
	\item{correctedRT}{
		Vector of corrected RT values for each measurement scan. If no RT correction was applied raw RT's should be used.
	}
	\item{uncorrectedRT}{
		Vector of raw RT values for each measurement scan.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
	Matrix containing 5 clolumns: rawRT, correctedRT, scan_id, mz, intnesity.
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
}

\seealso{
  \code{\link{PeakML.Read}, \link{PeakML.GapFiller}, \link{PeakML.xcms.write.SingleMeasurement}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
