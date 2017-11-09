\name{RCreateXMLDB}
\alias{RCreateXMLDB}
\title{RCreateXMLDB}
\description{
This exports data frame from R to xml database format supported by mzmatch program routines.
}
\usage{
RCreateXMLDB(data,outputfile)
}

\arguments{
  \item{data}{
Data frame which has a columns named: id, name, formula. These fields are required fields. If data frame has an extra columns named: inchi, description and synonyms then these fields also will be added to resulting xml file automatically.
}
  \item{outputfile}{
Name of the output XML file which should be generated.
}
}
\details{
Function can be used to generate XML format compound databases for mzmatch from spreedsheet like csv files (Excell spreedsheet for example). Input data should contain columns with names: id, name, formula. "id" must be unique identificator for the compound, "name" is a common name of the compound and "formula" is a molecular formula for the given compound. Additionaly columnns names "inchi" (should contain inchi strings for the given formulas), "description" (this can be used for labeling compound classes) and "synonyms". All other columns in input data set will be ignored and not included in ouptut XML file.
}
\value{
XML compounds database file.
}
\references{
PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
Andris Jankevics (a.jankevics@rug.nl)
Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
\code{\link{PeakML.xcms.read}}
}

\keyword{mzmatch.R}
\keyword{mzMatch}
