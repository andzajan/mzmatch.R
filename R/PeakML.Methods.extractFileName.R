PeakML.Methods.extractFileName <- function(path)
{
	string <-unlist(strsplit (path,"/"))
	out <- string[length(string)]
	out
}
