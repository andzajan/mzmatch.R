PeakML.Methods.getRTWindow <- function(retTime, time){
	# PRE:
	#	original retention time in sec
	# 	required window size
	# POST:
	# 	list containing the window (min, max)

	rv <- vector("list",2)
	rv[[1]]<- retTime-time
	rv[[2]]<- retTime+time
	rv
}
