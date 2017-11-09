mzmatch.XML.data.base.parser <- function (dbfile,elements=c("name","inchi","class","synonyms"))
{
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	DBcontent <- xmlInternalTreeParse (dbfile)
	ids <- getNodeSet(DBcontent, "/compounds/compound/id")
	ids <- sapply (ids,xmlValue)
	formulas <- getNodeSet(DBcontent, "/compounds/compound/formula")
	formulas <- sapply (formulas,xmlValue)
	OUTPUT <- data.frame (id=ids,formula=formulas)
	masses <- rep(NA,length(formulas))
	for (formula in 1:length(formulas))
	{		
		masses[formula] <- .jcall (project,returnSig="D",method="formulaToMass",as.character(formulas[formula]))
	}
	OUTPUT$mass <- masses
	OUT <- matrix(ncol=length(elements),nrow=length(ids))
	for (element in 1:length(elements))
	{
		value <- getNodeSet(DBcontent, paste("/compounds/compound/",elements[element],sep=""))
		value <- sapply (value,xmlValue)
		if (length(value)!=0)
		{
			OUT[,element] <- value
		}
	}
	colnames(OUT) <- elements	
	OUTPUT <- cbind (OUTPUT,OUT)
}
