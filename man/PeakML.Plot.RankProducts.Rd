\name{PeakML.Plot.RankProducts}
\Rdversion{1.1}
\alias{PeakML.Plot.RankProducts}
\title{PeakML.Plot.RankProducts}
\description{
	The rank product is a biologically motivated test for the detection of differentially expressed genes in replicated microarray experiments. It is a simple non-parametric statistical method based on ranks of fold changes. In addition to its use in expression profiling, it can be used to combine ranked lists in various application domains, including proteomics, metabolomics, statistical meta-analysis, and general feature selection. This function is using output from PeakML.RankProducts function to generate a pdf document, containing extracted chromatograms for each peak, show pfp and Erp values, box plots and all related peak set information. 
}
\usage{
PeakML.Plot.RankProducts (pdffile, PeakMLtree, RankProducts, identifications, sortby="Class1")
}
\arguments{
\item{pdffile}{
	Filename for a pdf file to be created. Text file containing tab separated values will be also created with the same file name.
}
\item{PeakMLtree}{
	Output from a "PeakML.Read" command, to be used to generate chromatograms plot, peak annotations and all relevant sample information.
}
\item{RankProducts}{
	Name of the ouptut object from a "PeakML.RankProducts" function.
}
\item{identifications}{
	R object containing compount idetifications, must be in the form as output from "PeakML.Methods.DBidToCompoundName" function.
}
\item{sortby}{
	Identifies if peaks in pdf output class should be sorted on Erp(expected rank product) value or "Class1" of "Class2".
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
	Output un pdf and text file format.
}
\references{
  Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P. (2004) Rank Products: A simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, FEBS Letters, 573:83--92
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Plot.RankProducts}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
