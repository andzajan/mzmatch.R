PeakML.Isotope.DB2Text <- function (filename, databases){
	
	dbases <- dir(paste(find.package("mzmatch.R"),"/dbs",sep=""), full.names=TRUE)
	#lapply(dbases, function(x) tail(strsplit(x,"/")[[1]],1)) # incase want to work on both database name and number.
	dbases <- dbases[databases]
	
	## Read data base as table.
	DBcont <- NULL
	for (i in 1:length(dbases))
	{
		DB <- mzmatch.XML.data.base.parser (dbfile = dbases[i], elements = c("name", "class", "subclass"))
		dbname <- sub(paste(find.package("mzmatch.R"),"/dbs/",sep=""),"",dbases[i])
		dbname <- sub(".xml","",dbname)
		DB$db <- dbname
		DBcont <- rbind(DBcont,DB)
	}
	
	neededDB <- unique (DBcont[,7])
	
	Annotations <- PeakML.get.attributes(filename,attribute="getGroupAnnotation",annotations=c("ppm","identification"))
	PeakTable <- PeakML.get.attributes(filename, attribute="CompleteTable")
	idindexes <- vector("list",length(Annotations[[2]]))
	for (ind in 1:length(idindexes))
	{
		if(!is.na(Annotations[[2]][ind]))
		{
			ann <- unlist(strsplit(Annotations[[2]][ind],", "))
			indices <- which(DBcont[,1]%in%ann)
			idindexes[[ind]] <- indices
		}
	}
	## Check out whic idindexes are empty to exluce compounds which are not identified in DB's of the interest.

	identified <- sapply (1:length(idindexes),function (x) length(idindexes[[x]]))
	identified <- which (identified!=0)

	intensities <- PeakTable[[1]][,identified]
	masses <- PeakTable[[2]][,identified]
	masses <- apply(masses,2,mean,na.rm=TRUE)
	RTs <- round(PeakTable[[3]][,identified],0)
	RTs <- apply(RTs,2,mean,na.rm=TRUE)
	formulas <- rep(NA,length(masses))
	theormasses <- NULL
	compnames <- NULL
	dblist <- NULL
	classes <- NULL
	subclasses <- NULL

	for (index in 1:length(identified))
	{
		listnum <- identified[index]
		datatable <- DBcont[idindexes[[listnum]],]
		if (nrow(datatable)>1)
		{
			formula <- unique (datatable[,"formula"])
			formula <- sub("\\[M1];\\[","",formula)
			formula <- sub("\\]n","",formula)
			formula <- unique(formula)
			formula <- paste(formula,collapse=", ")
			mass <- paste(unique (datatable[,"mass"]),collapse=", ")
			compname <- paste(unique (datatable[,"name"]),collapse=", ")
			detectedin <- paste(unique (datatable[,"db"]),collapse=", ")
			class <- unique (datatable[,"class"])
			class <- sub("\\[","",class)
			class <- sub("\\]","",class)
			class <- paste(class,collapse=", ")
			subclass <- unique (datatable[,"subclass"])
			subclass <- sub("\\[","",subclass)
			subclass <- sub("\\]","",subclass)
			subclass <- paste(subclass,collapse=", ")
		} else
		{
			formula <- as.character(datatable$formula)
			formula <- sub("\\[M1];\\[","",formula)
			formula <- sub("\\]n","",formula)
			mass <- datatable$mass
			compname <- as.character(datatable$name)
			detectedin <- as.character(datatable$db)
			class <- as.character(datatable$class)
			subclass <- as.character(datatable$subclass)
			class <- sub("\\[","",class)
			class <- sub("\\]","",class)
			subclass <- sub("\\[","",subclass)
			subclass <- sub("\\]","",subclass)
		}
		if(length(formula)!=0){
			formulas[index] <- formula
		} else{
			formulas[index] <- NA
		}
		
		theormasses <- append(theormasses,mass)
		compnames <- append(compnames,compname)
		dblist <- append(dblist,detectedin)
		classes <- append(classes,class)
		subclasses <- append(subclasses,subclass)
	}

	ppmerr1 <- NULL
	for (massnr in 1:length(identified))
	{
		TheorM <- as.numeric(unlist(strsplit(theormasses[massnr],", ")))
		MeasM <- masses[identified][massnr]
		ppm <- paste(round((TheorM-MeasM)/TheorM*10^6,1),collapse=", ")
		ppmerr1 <- append(ppmerr1,ppm)
	}


	identtable <- data.frame(mass=masses, rt=RTs, formula=formulas, ppm=ppmerr1, id=Annotations[[2]][identified], name=compnames, class=classes, subclass=subclasses)
	#identtable <- data.frame(formula=formulas, rt=RTs, id=Annotations[[2]],identifications=compnames)
	inttable <- data.frame(t(PeakTable[[1]][,identified]),check.names=FALSE)
	identtable <- cbind(identtable, inttable)
	rownames (identtable) <- c(1:nrow(identtable))
	write.table(file=paste(filename,".csv",sep=""),identtable,sep="\t",row.names=FALSE)
	identtable
	#hwrite (identtable,paste(filename,".html",sep=""))
}

