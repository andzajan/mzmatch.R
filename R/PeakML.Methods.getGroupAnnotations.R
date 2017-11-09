PeakML.Methods.getGroupAnnotations <- function(PeakMLtree=NULL,filename=NULL)
{
	if (!is.null(filename))
	{
		st <- system.time(PeakMLtree <- xmlInternalTreeParse(filename))
	}

	setidsrecord <- function (setid)
	{
		STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/peaks/peak[",setid,"]/annotations/annotation/label\"),xmlValue)",sep="")
		annotationnames <- eval (parse(text=STR))
		STR <- paste ("sapply(getNodeSet(PeakMLtree,\"/peakml/peaks/peak[",setid,"]/annotations/annotation/value\"),xmlValue)",sep="")
		annotationvalues <- eval (parse(text=STR))
		if (length(annotationvalues)!=0)
		{
			val <- cbind (setid,annotationnames,annotationvalues)
		} else
		{
			val <- NULL
		}
		val
	}

	nrpeakgroups <- as.numeric(sapply(getNodeSet(PeakMLtree,"/peakml/header/nrpeaks"),xmlValue))

	system.time(annots <- lapply (1:nrpeakgroups,setidsrecord))
	annots <- do.call (rbind,annots)

	annotnames <- unique (annots[,2])

	groupannotations <- matrix (ncol=length(annotnames),nrow=nrpeakgroups)

	rownums <- as.numeric(annots[,1])

	system.time({
	for (rn in 1:length(rownums))
	{
		groupannotations[rownums[rn],which(annotnames==annots[rn,2])] <- annots[rn,3]
	}
	})
	
	if (length(rownums)!=0)
	{
		attr <- vector ("list",length(annotnames))
		names(attr) <- annotnames
		for (annotid in 1:length(attr))
		{
			attr[[annotid]] <- groupannotations[,annotid]
		}
	} else
	{
		attr <- NULL
	}

	attr
}
