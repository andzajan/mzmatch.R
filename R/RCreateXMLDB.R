RCreateXMLDB <- function (data,outputfile)
{
	nodeGen <- function(nodenum)
	{
		compound <- xmlNode ("compound",xmlNode ("id", data[nodenum,which(childnames=="id")]))
		compound[[2]] <- xmlNode ("name", data[nodenum,which(childnames=="name")])
		compound[[3]] <- xmlNode ("formula", data[nodenum,which(childnames=="formula")])

		## Add all the rest nodes beside default ones
		notfilled <- which(!childnames%in%c("id","formula","name"))
		if (length(notfilled)!=0)
		{
			for (nd in 1:length(notfilled))
			{
				compound[[notfilled[nd]]] <- xmlNode (childnames[notfilled[nd]],data[nodenum,notfilled[nd]])
			}
		}
		compound
	}	
	childnames <- colnames (data)

	if ("id"%in%childnames & "name"%in%childnames & "formula"%in%childnames)
	{
		## For every row in input data table, make child nodes with subnodes
		##XMLfile <- newXMLDoc ("treedoc")
		compounds <- xmlNode("compounds")
			
	system.time(compoundlist <- lapply(1:nrow(data),nodeGen))
	compounds <- addChildren(compounds,kids=compoundlist)
	saveXML(doc=compounds,file=outputfile,prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n \n")
	} else
	{
		cat ("Values(column names) for fields: id,name and formula must be set.\n")
	}
}

