\name{PeakML.Write}
\Rdversion{1.1}
\alias{PeakML.Write}
\title{PeakML.Write}
\description{
  This function writes the given peak data to the file specified in the PeakML format. This file format
  supports the mass traces (features) and the sets-structure utilized by XCMS and as such provides a
  fully reproducible representation of the contents of an XCMS set structure.
}
\usage{
  PeakML.Write(peakMLdata=NULL, peakDataMtx=NULL, chromDataList=NULL, sampleNames=NULL, rawDataFullPaths=NULL, phenoData=NULL, corRT=NULL, rawRT=NULL, ionisation=NULL, GroupAnnotations=NULL, outFileName)
}
\arguments{
	\item{peakMLdata}{
	  The list of peak data as provided by the function PeakML.Read. If it is provided all other function arguments beside "outFileName" will be ignored.
	}
	\item{peakDataMtx}{
	  The peak data matrix as provided by the function PeakML.Read
	}
	\item{chromDataList}{
	  The list of chromatograms as provided by the function PeakML.Read
	}
	\item{sampleNames}{
	  The sample names as a character vector.
	}
	\item{rawDataFullPaths}{
	  Paths to the raw data files as a character vector.
	}
	\item{phenoData}{
	  The phenotypes as a character vector.
	}
	\item{corRT}{
	  The corrected retention times as list.
	}
	\item{rawRT}{
	  The raw retentiontimes as list.
	}
	\item{outFileName}{
	  The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
	  will fail.
	}
	\item{ionisation}{
	  Character value to be set on to "positive" or "negative" if the mass has to be corrected accordingly, or leave it as "neutral".
	}
	\item{GroupAnnotations}{
	  List of annotations which should be written to PeakML file. List names should contain annotation names. And each list should contain a vector of the length equal to number of peak groups. NA values should be used if annotation is not present in given group.
	}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.xcms.write function extends the XCMS package with write functionality for exporting data to other applications.
  For example, the results of an XCMS step can then be loaded in the user interface environment PeakML Viewer and manipulated
  there.
}
\value{
  This function returns no value, but writes a PeakML file from the the given data.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Achuthanunni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getPeakData}}
}

\keyword{mzmatch.R}
\keyword{ XCMS }
