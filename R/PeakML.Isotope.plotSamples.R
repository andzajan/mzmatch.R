PeakML.Isotope.plotSamples <- function(isotopeChroms, metName, metFormula, metMass, metComment, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon, label, exclude_from_plots){

  cat (exclude_from_plots)
  cat("\n")
	element <- substr(label,1,1)
	if(is.null(stdRT)) stdRT <- NA
#	grads <- c("black",
#rgb(0.750,0.700,1.000),  rgb(0.400,0.320,0.800), rgb(0.150,0.060,0.600),
#rgb(0.700,0.900,1.000),  rgb(0.320,0.640,0.800), rgb(0.060,0.420,0.600),
#rgb(0.900,1.000,0.700),  rgb(0.640,0.800,0.320), rgb(0.420,0.600,0.060),
#rgb(1.000,0.850,0.700),  rgb(0.800,0.560,0.320), rgb(0.600,0.330,0.060),
#rgb(1.000,0.700,0.700),  rgb(0.800,0.320,0.320), rgb(0.600,0.060,0.060)
#	)

#	cbp <- c("black", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

	cbp <- c("#666666", "#1F78B4", "#6A3D9A", "#33A02C", "#B15928", "#FF7F00", "#E31A1C", "#A6CEE3", "#CAB2D6", "#B2DF8A", "#FFFF99", "#FDBF6F", "#FB9A99")

#	dist_grads <- c("black", rgb(1.000,0.750,0.500), rgb(1.000,0.500,0.000), rgb(1.000,1.000,0.600), rgb(0.700,1.000,0.550), rgb(0.200,1.000,0.000), rgb(0.650,0.930,1.000), rgb(0.100,0.700,1.000), rgb(0.800,0.750,1.000), rgb(0.400,0.300,1.000), rgb(1.000,0.600,0.750), rgb(1.000,1.000,0.200), rgb(0.900,0.100,0.200))

	clrs <- cbp
	numClrs <- length(clrs)

	axExFact <-1.1 # This is a buffer to get the y axis of the bar plots right
	par(oma=c(0,0,0,0))

	processRatioMtx <- function(ratioMtx){
		plotMtx <- matrix(nrow=nrow(ratioMtx), ncol=ncol(ratioMtx))
		rownames(plotMtx) <- row.names(ratioMtx)
		for (r in 1:nrow(ratioMtx)){
			sumr <- sum(ratioMtx[r,])
			for (c in 1:ncol(ratioMtx)){
				if (!is.na(ratioMtx[r,c])) {
					plotMtx[r,c] <- ratioMtx[r,c]/sumr
				} else {
					plotMtx[r,c] <- 0
				}
			}
		}
		plotMtx
	}


	mzList <- isotopeChroms[[1]]
	intList <- isotopeChroms[[2]]
	rtList <- isotopeChroms[[3]]

	numPeakGroups <- length(mzList)
	numSampleGroups <- length(mzList[[1]])
	numCarbons <- length(mzList[[1]][[1]])
	#numReplicates <- length(mzList[[1]][[1]][[1]])

	metName <- unlist(strsplit(as.character(metName), ", "))

	fillLabels <- c("UL" , paste(as.character("+"),c(1:numCarbons),sep=""))# create a list of isotop names

	fillColor <- clrs[1:numCarbons]

	if (numCarbons > numClrs) {
		n = numCarbons - numClrs
		rainClrs <- rainbow(n, s = 1, v = 1, start = 0, end = max(1,n - 1)/n, alpha = 1)
		fillColor <- c(clrs, rainClrs)
	}

	if (followCarbon > numCarbons){
		followCarbon <- numCarbons
	}

	for (peakGroup in 1:numPeakGroups){
		trendList <- PeakML.Isotope.getTrendList (intList, sampleGroups, useArea) [[peakGroup]]
		ratioMtx <- processRatioMtx (t(PeakML.Isotope.getRatioMtxList (intList, sampleGroups, useArea, metName[1]) [[peakGroup]]))

		par (mar=c(0,0,0,0))
		plot (c(1:10), c(1:10), xlab="", ylab="", pch="", axes=F)

		if(length(metName)==1){
			text(1.5, 6, metName[1], cex = 2, pos=4)
			text(1.5, 4, paste("Formula:", metFormula, metComment, "Mass:", round(metMass,3), "Std.RT:", stdRT, "Ion:", sampleType, sep = "  "), pos=4, cex = 1.5)
#			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, sep = "  "), pos=4, cex = 1.5)
		} else {
			text(1.5, 7, metName[1], cex = 2, pos=4)
			text(1.5, 5, paste(metName[2:length(metName)], collapse=","), cex = .5, pos=4)
			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, "Ion:", sampleType, sep = "  "), pos=4, cex = 1)
#			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, sep = "  "), pos=4, cex = 1)
		}
		legend(1.5, 3, fill=fillColor, fillLabels[1:numCarbons], bty="n", horiz=TRUE, x.intersp=.1)

		plot (1, 1, xlab="", ylab="", pch="", axes=F)
		text(1, 1, paste("G", as.character(peakGroup), sep=""), cex = 3)

		for (item in 1:length(plotOrder)){
			if(plotOrder[item] == "LEGEND"){

				#par(mar=c(0.5,0.5,0.5,0.5))
				par (mar=c(0,0,0,0))
				plot (1, 1, xlab="", ylab="", pch="", axes=F)
				legend("topleft",fill=fillColor, fillLabels[1:numCarbons], bty="n")

			} else if (plotOrder[item] == "TREND"){

				errBar <- function(x, y, upper, lower=upper, length=0.1,...){
					if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
					stop("vectors must be same length")
					arrows(x, y+upper, x, y-lower, angle=90, code=3, length=length, ...)
				}

				trendMtx <- PeakML.Isotope.getTrendMtx(trendList, sampleGroups)

                                if(!is.null(exclude_from_plots)){
                                  exIndex <- match(exclude_from_plots, colnames(trendMtx))
                                  trendMtx <-subset(trendMtx,select=-exIndex)
                                  lbls <- setdiff(sampleGroups, sampleGroups[exclude_from_plots])
                                }


                               # if(!is.null(exclude_from_plots)){
                               #   exclude <- which(sampleGroups %in% exclude_from_plots)
                               #   trendMtx <- trendMtx[,- exclude]
                               #   lbls <- setdiff(sampleGroups, sampleGroups[exclude_from_plots])
                    #
  #                              }

                                trendMtx.sum <- apply(trendMtx, 2, sum)
				trendMtx.sd <- apply(trendMtx, 2, sd)
				ylimit <- sum(max(trendMtx.sum), max(trendMtx.sd))
				# calculate std err
				#errMtx <- replace(trendMtx, trendMtx==0, NA)
				#stdErr <- apply(errMtx, 2, PeakML.Methods.getStdErr)

				#par (mar=c(2,4,2,0))
				par (mar=c(2,4,2,0), mgp=c(2,1,0))
#				par(mar=c(4,2.5,0.5,3))
				if (useArea == FALSE){
					#ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
					ylabel <- "mean peak height"
				} else {
					#ylabel <- paste(metName[1], "(mean peak area)", sep= " ")
					ylabel <- "mean peak area"
				}

				mp<-barplot(trendMtx, beside=FALSE, col=fillColor, ylab=ylabel, ylim=c(0,ylimit), border=NA, axisnames=FALSE, axes=FALSE)
				text(mp, par("usr")[3], labels = colnames(trendMtx), srt = 45, adj = 1, xpd =TRUE, cex=0.6)
				errBar(mp, trendMtx.sum, trendMtx.sd/numCarbons, lwd=.3)
				title("Trend plot", cex.main=0.8)
				axis(2)


			} else if (plotOrder[item] == "LABELLED"){

				if(is.na(followCarbon)) followCarbon <- 2

                        	if (followCarbon == 2){
					fcMtx <-PeakML.Isotope.getFCMtxAbun(trendList, sampleGroups, followCarbon, element)
                                  #fcMtx <-PeakML.Isotope.getFCMtxAbun(trendList, sg, followCarbon, element)
				} else {
					fcMtx <- PeakML.Isotope.getFCMtx(trendList, sampleGroups, followCarbon)
                                  #fcMtx <- PeakML.Isotope.getFCMtx(trendList, sg, followCarbon)
				}


                                if(!is.null(exclude_from_plots)){
                                                 exIndex <- match(exclude_from_plots, colnames(fcMtx))
                                                 fcMtx <-subset(fcMtx,select=-exIndex)
                                                 lbls <- setdiff(sampleGroups, sampleGroups[exclude_from_plots])
                                 }

				if (useArea == FALSE){
#					ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
					ylabel <- "mean peak height"
				} else {
#					ylabel <- paste(metName[1], "(mean peak area)", sep= " ")
					ylabel <- "mean peak area"
				}

                        	ylimit <- max(apply(fcMtx, 2, sum)) * axExFact

				if (!ylimit==0){

#					par (mar=c(2,4,2,0))
					par (mar=c(2,4,2,0), mgp=c(2,1,0))
#					par(mar=c(4,2.5,0.5,3))
					if (followCarbon==2){
						overlap <- "#993399" # expected NA / NA whichever is smaller
						expectedLarger <- "#FF9988" # expected natuabun is greater
						measuredLarger <- "#8877FF" # blue
						colvector <- c(overlap, expectedLarger, measuredLarger)

					} else {
						colvector <- fillColor[followCarbon]
					}

#					barplot(fcMtx, beside=FALSE, col=colvector, axisnames=FALSE, ylab=ylabel, ylim=c(0,ylimit), border=NA)
#					axis(1, las=3, at=c(1:length(sampleGroups)), labels=sampleGroups, lwd=0, cex.axis=.8)

					mp<-barplot(fcMtx, beside=FALSE, col=colvector, ylab=ylabel, ylim=c(0,ylimit), border=NA, axisnames=FALSE, axes=FALSE)
					text(mp, par("usr")[3], labels = colnames(fcMtx), srt = 45, adj = 1, xpd =TRUE, cex=0.6)
#					axis(1, at = mp, labels = FALSE)
					title(paste("Trend of ", followCarbon-1,"", element, " labelled isotopomer"), cex.main=0.8)
					#title(paste("Trend of natural abundance (", element, ")"), cex.main=0.8)
					axis(2)

					if (followCarbon==2){
						legend("topright",fill=colvector, c("Overlap", "<Expected", ">Expected"), bty="n")
					}
				} else {
#					par (mar=c(2,4,2,0))
					par(mar=c(2,4,2,0), mgp=c(2,1,0))
#					par(mar=c(0.5,0.5,0.5,0.5))
					plot (1, 1, xlab="", ylab="", pch="", axes=F)
				}

			} else if (plotOrder[item] == "RATIO") {

#				par (mar=c(2,4,2,0))
				par(mar=c(2,4,2,0), mgp=c(2,1,0))



				mp<-barplot(t(ratioMtx), beside=FALSE, col=fillColor, ylab="% area under peak", ylim = c(0,1), border=NA, axisnames=FALSE, axes=FALSE)
				text(mp, par("usr")[3], labels = row.names(ratioMtx), srt = 45, adj = 1, xpd =TRUE, cex=0.2)
#				axis(1, at = mp, labels = FALSE)
				title("Ratio", cex.main=0.8)
				axis(2)


                     } else if (plotOrder[item] == "TOTRATIO") {
				getRelAbunMtx <- function (trendList, sampleGroups, followCarbon){
					numCarbons <-  length(trendList[[1]])
					plotMtx <- matrix(nrow = 1, ncol = length(sampleGroups))
					dimnames(plotMtx) <- list(c("RelAbun"), sampleGroups)

					for (sam in 1:length(sampleGroups)){
						for (row in 1:length(rownames(plotMtx))){
							x <- trendList[[sam]][[followCarbon]]
							y <- sum(unlist(trendList[[sam]]))

							if (is.null(x)){
								x <- 0
							}
							if (is.null(y)){
								y <- 0
							}

							plotMtx[row, sam] <- x/y*100
						}
					}
					plotMtx
				}


				if(is.na(followCarbon)) followCarbon <- 2


				ftMtx <- getRelAbunMtx(trendList, sampleGroups, followCarbon)
				ftMtx[is.nan(ftMtx)] <- 0

                               if(!is.null(exclude_from_plots)){
                                  exIndex <- match(exclude_from_plots, colnames(ftMtx))
                                  ftMtx <-subset(ftMtx,select=-exIndex)
                                  lbls <- setdiff(sampleGroups, sampleGroups[exclude_from_plots])
                                }

				if (useArea == FALSE){
                    #					ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
					ylabel <- "% relative labelling"
				} else {
					ylabel <- "% relative labelling"
					#ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
				}


				ylimit <- max(apply(ftMtx, 2, sum)) * axExFact
				if (ylimit>100) ylimit <- 100
                #				if (is.nan(ylimit)){ylimit<-0}

				if (!ylimit==0){

                    #					par (mar=c(2,4,2,0))
					par (mar=c(2,4,2,0), mgp=c(2,1,0))
                    #					par(mar=c(4,2.5,0.5,3))

					colvector <- fillColor[followCarbon]

					mp<-barplot(ftMtx, beside=FALSE, col=colvector, ylab=ylabel, ylim=c(0,ylimit), border=NA, axisnames=FALSE, axes=FALSE)
					text(mp, par("usr")[3], labels=colnames(ftMtx), srt = 45, adj = 1, xpd =TRUE, cex=0.6)
                    #					axis(1, at = mp, labels = FALSE)
					title("Relative labelling pattern", cex.main=0.8)
					axis(2)
				} else {
                    #					par (mar=c(2,4,2,0))
					par(mar=c(2,4,2,0), mgp=c(2,1,0))
                    #					par(mar=c(0.5,0.5,0.5,0.5))
					plot (1, 1, xlab="", ylab="", pch="", axes=F)
				}




			} else if (plotOrder[item] == "EMPTY"){

				# This is a an empty plot for flexibility e.g. define empty if a plot is not needed
				plot (1, 1, xlab="", ylab="", pch="", axes=F)
			} else{
				if (! length(unlist(rtList[[peakGroup]][[item]])) == 0){
					maxRT <- max(unlist(rtList[[peakGroup]][[item]]),na.rm=TRUE)
					minRT <- min(unlist(rtList[[peakGroup]][[item]]),na.rm=TRUE)
					maxIN <- max(unlist(intList[[peakGroup]][[item]]),na.rm=TRUE)
					minIN <- 0
				}else{
					maxRT <- 0
					minRT <- 0
					maxIN <- 0
					minIN <- 0
				}

#				par(mar=c(4,4,0.5,0.5), mgp=c(1.5,0.5,0))
				par(mar=c(2,2,1,0), mgp=c(0.5,0.5,0))
				plot (1, 1, pch="", xlab="", ylab="", xlim=c(minRT,maxRT), ylim=c(minIN,maxIN), cex.axis=0.75)
				for (isotop in 1:numCarbons){
					numReplicates <- length(rtList[[peakGroup]][[item]][[isotop]])
					for (rep in 1:numReplicates){
						points(rtList[[peakGroup]][[item]][[isotop]][[rep]], intList[[peakGroup]][[item]][[isotop]][[rep]], type="l",col=fillColor[[isotop]])
					}
				}
				legend("topright", plotOrder[item], bty="n", cex=.5)
			}
		}
	}
}
