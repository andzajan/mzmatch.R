\name{mzmatch.XML.data.base.parser}
\alias{mzmatch.XML.data.base.parser}
\title{mzmatch.XML.data.base.parser}
\description{
This function will read content of mzmatch compound databese file in to R data frame structure.
}
\usage{
mzmatch.XML.data.base.parser (dbfile,elements=c("name","inchi","class","synonyms"))
}

\arguments{
  \item{dbfile}{
Database file in XML format which should be converted.
}
  \item{elements}{
Character vector of nodes which should be read from the XML compound database file.
}
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

\keyword{ peakML }
\keyword{ mzmatch.R }
\keyword{ mzMatch }
