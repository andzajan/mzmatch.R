PeakML.Isotope.ReadDB <- function (dbNames){
#	"""
#	Extends mzmatch.XML.data.base.parser to read a database based on its name
#	PRE: database name, extra elements to be read (for extensibility)
#	POST: data.frame with id name formula
#	"""
	
	DB <- NULL
	dbases <- dir(paste(find.package("mzmatch.R"),"/dbs",sep=""), full.names=TRUE) # mzmatch databases
	
	for (db in 1:length(dbases)){
		dbname <- strsplit(tail(strsplit(dbases[[db]], "/")[[1]], n=1), ".xml")[[1]][[1]]
		if (dbname %in% dbNames){
			dbase <- mzmatch.XML.data.base.parser (dbfile = dbases[[db]])
			dbase <- data.frame(id = dbase$id, name  = dbase$name, formula = dbase$formula)
			DB <- rbind(DB, dbase)
		}
	}
	DB
}
