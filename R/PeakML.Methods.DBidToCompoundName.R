PeakML.Methods.DBidToCompoundName <- function (DBS,PeakMLdata,collapse=TRUE)
{	
	dbnameext <- function (i)
	{
		A <- unlist(strsplit(DBS[i],"/"))
		A <- A[length(A)]
		A <- sub(".xml","",A)
		A
	}

	annot.extract <- function(annot)
	{
		if (!is.na(PeakMLdata$GroupAnnotations$identification[annot]))
		{
			cat (annot,"\n")
			identifications <- unlist(strsplit(PeakMLdata$GroupAnnotations$identification[annot],", "))
			ppms <- unlist(strsplit(PeakMLdata$GroupAnnotations$ppm[annot],", "))
			adducts <- unlist(strsplit(PeakMLdata$GroupAnnotations$adduct[annot],", "))
			dat <- DBcont[which(as.character(DBcont[,1])%in%identifications),]
			dat[,2] <- sub ("\\[M1\\];\\[","", dat[,2])
			dat[,2] <- sub ("\\]n","", dat[,2])

			# Keep the same id row for adduct and ppm

			ppmout <- rep(NA, length(identifications))
			adductout <- rep(NA, length(identifications))
			
			for (idn in 1:length(identifications))
			{
				hit <- which(identifications==dat[idn,1])[1]
				adductout[idn] <- adducts[hit]
				ppmout[idn] <- ppms[hit]
			}
			dat$ppm <- round(as.numeric(ppmout),1)
			dat$adduct <- adductout

			if (collapse==TRUE)
			{
				unique.formulas <- unique(dat[,2])
				out <- NULL
				for (uniq in 1:length(unique.formulas))
				{
					hits <- which(dat[,2] == unique.formulas[uniq])
					out <- rbind(out,c(paste(dat[,1],collapse=", "),unique.formulas[uniq],dat[hits[1],3],paste(dat[,4],collapse=", "),paste(dat[,5],collapse=", ")))
				} 
			} else
			{
				out <- dat
			}
		} else
		{
			out <- NA
		}
		out
	}

	dbnames <- sapply(1:length(DBS),dbnameext)

	DBcont <- NULL
	for (i in 1:length(DBS))
	{
		DB <- mzmatch.XML.data.base.parser (dbfile=DBS[i],elements=c("name"))
		DB$db <- dbnames[i]
		DBcont <- rbind(DBcont,DB)
	}

	# For debug purposes only
	#for (i in 1:length(PeakMLdata$GroupAnnotations$identification))
	#{
	#	cat (i,"\n")
	#	annot.extract(i)
	#}

	id.resolved <- lapply (1:length(PeakMLdata$GroupAnnotations$identification),annot.extract)
	id.resolved
}



