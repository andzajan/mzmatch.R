\name{PeakML.RankProducts}
\Rdversion{1.1}
\alias{PeakML.RankProducts}
\title{PeakML.RankProducts}
\description{
	The rank product is a biologically motivated test for the detection of differentially expressed genes in replicated microarray experiments. It is a simple non-parametric statistical method based on ranks of fold changes. In addition to its use in expression profiling, it can be used to combine ranked lists in various application domains, including proteomics, metabolomics, statistical meta-analysis, and general feature selection. This rank products test implements effective algorithm to carry out unpaired sample test. 
}
\usage{
PeakML.RankProducts (sampleGroups, inputTable, groupA, groupB, numberOfPermutations=1000, calculateProduct=TRUE, MinNumOfValidPairs=NA, replaceNA=FALSE, RandomPairs=NA)
}
\arguments{
\item{sampleGroups}{
	Character vector with a length equal to the number of sample in the data set. Normally output from the PeakML.Read() command. PeakMLtree$phenoData
}
\item{inputTable}{
	Table with a peak intensities. Each row should represent a single sample and each column a single peak. Normally outpu from PeakML.Methods.getCompleteTable command. 
}
\item{groupA}{
	 As typical peakml file can contain more than 2 sample classes, classes to compare should be specified with this function parameter. A character string giving a name fro Class 1 to be considered.
}
\item{groupB}{
	A character string with a name for the second sample class.
}
\item{numberOfPermutations}{
	Number of permutations used to construct null distribution, which is used to calculate p-values. Default value is 1000.
}
\item{calculateProduct}{
	If set to "TRUE", geometric mean of ranks is used during test (rank product), if set to "FALSE" sum of ranks is used (ranks sum).
}
\item{MinNumOfValidPairs}{
	Number of non-NA pairs to be detected, to include a peak in the test. Is set to "NA", n/2 is used, where n is a number of samples in the smallest group. If peak is discarded from the test, test values and p-values are set to "NA".
}
\item{replaceNA}{
	If set to "TRUE" all missing values in input data set are replaced with the lowest detected signal in whole data set. 
}
\item{RandomPairs}{
	Number of random assignment of the sample pairs. If set to NA n^2 is used, where n is a number of the samples in the smallest group.
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
	Rank products objetc, containing a ratios, pfp values and expected rank products. PeakML.Plot.RankProducts function can be used to generate pdf output file.
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
