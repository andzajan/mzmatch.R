\name{mzmatch.R.Setup}
\Rdversion{1.1}
\alias{mzmatch.R.Setup}
\title{mzmatch.R.Setup}
\description{
  This function is used to setup a data files folder, and define filenames, sample classes etc.
}
\usage{
  mzmatch.R.Setup (projectFolder=NULL, samplelist=NULL, outputfolder="peakml", sampleClassFromFolderNames = FALSE,
                 extractMeasurentDate = FALSE)
}
\arguments{
\item{projectFolder}{
  A name of the root folder where all data files and resulting output files are/will be stored. If set to 'NULL' user will be asked to navigate to proper folder with GUI.
}
\item{samplelist}{
  A tab limited text file containing at least 2 columns and sample related information. "filename" - full path and file name of the data file in mzXML. mzData or mzML format, "sampleClass" - character string defining a type of the sample and it's class.. If value is set to 'NULL' all "mzXML", "mzData" or "mzML" files present in the current data folder will be listed and user will be able to edit a sample class information through graphical interface.
}
\item{outputfolder}{
  Character string containin a name of the folder where output files in peakML format will be created.
}
\item{sampleClassFromFolderNames}{
  I data files are ordered in sub folders, a name of the subfolder automatically can be used as a class label.
}
\item{extractMeasurentDate}{
  Extract measuremnt date and time information from mzML data files. This information can be used as input for QC correction methods.
}


}
\details{
 To be able to call java library within R environment, the user should define a reserved amount of operating memory (RAM) for the Java virtual machine. Typically this value should be 1024 megabytes. On several systems initialization of the Java virtual machine can fail (if the machine has less than 1GB of RAM or more than 3GB, this sometimes happens on computers with the Windows XP operating system). The user can change the memory limits to for example 256,512,768,2048, etc.
}
\value{
  Object "sampleList" created in the R .GlobalEnv environment. A data frame containing at least 4 columns. "filename" - full path and file name of the data file in mzXML. mzData or mzML format, "sampleClass" - character string defining a type of the sample and it's class. "globalClass" - secondary sample class label (optional), can be used for example to define quality control samples and "real" data files. "outputfilenames" - a full path and file name of the outputfiles after peak picking step in peakml file format.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Andris Jankevics (andris.jankevics@manchester.ac.uk)
}

\seealso{
  \code{\link{mzmatch.init}}
}

\keyword{ mzmathc.R }
\keyword{ mzMatch }
