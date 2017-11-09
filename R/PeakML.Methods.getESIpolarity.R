PeakML.Methods.getESIpolarity <- function(filename)
{
	
	# Binary file, like cdf format, can't be parsed by XML.
	tryCatch(xmltree <- xmlParse(filename), 
		error= function(e) {cat ("ESI polarity mode can't be extracted from files which are not in mzXML or mzML formats. Polarity is set to neutral. \n")}
	)

	if (exists("xmltree",environment()))
	{
		xmlFileType <- xmlName(xmlRoot(xmltree))

		if (xmlFileType=="indexedmzML" | xmlFileType=="mzML")
		{
			#x reffers there to the default namespace of the document				
			negpolarity <- getNodeSet(xmlRoot(xmltree), "//x:cvParam[@name='negative scan']", "x")
			pospolarity <- getNodeSet(xmlRoot(xmltree), "//x:cvParam[@name='positive scan']", "x")

			if (length(pospolarity)==0) {polarity <- "negative"}
			if (length(negpolarity)==0) {polarity <- "positive"}
			if (length(pospolarity)>0 && length(negpolarity)>0)
			{ 
				cat ("mzML file contains scans for more than one ESI polarity mode. Such data files are not supported. \n")
				stop ()
			}		
		}

		if (xmlFileType=="indexedmzXML" | xmlFileType=="mzXML")
		{
			polarityf <- unique(unlist(getNodeSet(xmlRoot(xmltree), "//x:scan/@polarity", "x")))
			if (length(polarityf)>1)
			{
				cat ("mzXML file contains scans for more than one ESI polarity mode. Such data files are not supported. \n")
				stop ()	
			}

			# In case if no + or - is present, set it to neutral			
			polarity <- "neutral"			

			if (polarityf=="+") 
			{
				polarity <- "positive"
			}
			if (polarityf=="-")
			{
				polarity <- "negative"
			}
		}

		rm (xmltree)
		gc ()
		
	} else
	{
		polarity <- "neutral"
	}	
	polarity
}
