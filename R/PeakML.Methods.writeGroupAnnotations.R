PeakML.Methods.writeGroupAnnotations <- function(project, GroupAnnotations)
{
	annotnames <- names(GroupAnnotations)
	
	for (annot in 1:length(annotnames))
	{
		for (groupnumber in 1:length(GroupAnnotations[[annot]]))
		{
			if (!is.na(GroupAnnotations[[annot]][groupnumber]))
			{
				.jcall(project, returnSig="V",method="addGroupAnnotation",as.integer(groupnumber-1), as.character(annotnames[annot]),as.character(GroupAnnotations[[annot]][groupnumber]))
			}
		}
	}
}

