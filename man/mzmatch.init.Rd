\name{mzmatch.init}
\Rdversion{1.1}
\alias{mzmatch.init}
\title{mzmatch.init}
\description{
  This function is used to define globals Java envirnoment parameters before startting mzmatch java library.
}
\usage{
  mzmatch.init(memorysize, version.1=TRUE)
}
\arguments{
\item{memorysize}{
  Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1024 megabytes (1GB).
}
\item{version.1}{
  For bacwkard compability it is possible to use mzmatch V1. If option version.1 is set to "FALSE", it will load the latest mzmatch version. PeakML files generated with mzmatch V.1 are not compatiable with V.2 and above and vice versa.
}
}
\details{
 To be able to call java library within R environment, the user should define a reserved amount of operating memory (RAM) for the Java virtual machine. Typically this value should be 1024 megabytes. On several systems initialization of the Java virtual machine can fail (if the machine has less than 1GB of RAM or more than 3GB, this sometimes happens on computers with the Windows XP operating system). The user can change the memory limits to for example 256,512,768,2048, etc.
}
\value{
  This function returns no value.
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

\keyword{ mzmathc.R }
\keyword{ mzMatch }
