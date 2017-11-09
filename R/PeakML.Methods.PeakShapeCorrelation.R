PeakML.Methods.PeakShapeCorrelation <- function(chromatogram1,chromatogram2)
{
	# Calculate correlation based on peak shape
	RT1 <- chromatogram1[3,]
	INT1 <- chromatogram1[2,]
	RT2 <- chromatogram2[3,]
	INT2 <- chromatogram2[2,]

	# Check if RT window overlap
	RT1range <- c(min(RT1),max(RT1))
	RT2range <- c(min(RT2),max(RT2))
	if (RT1range[2]<RT2range[2])
	{
		if (RT1range[2]>= RT2range[1])
		{ overlap="Y"} else {overlap="N"}
	} else
	{
		if (RT2range[2]>= RT1range[1])
		{ overlap="Y"} else {overlap="N"}
	}

	if (length(RT1)>4 & length(RT2)>4 & overlap=="Y")
	{
		# Estimate scanrate from measured data
		scanrate <- round(min(c(RT1[2:length(RT1)]-RT1[1:(length(RT1)-1)],RT2[2:length(RT2)]-RT2[1:(length(RT2)-1)])),1)
		# create a grid of RT which should overlap
		RTs <- seq(min (c(RT1,RT2)),max(c(RT1,RT2)),scanrate)

		INTmat <- matrix(nrow=2,ncol=length(RTs))
		for (d in 1:ncol(INTmat))
		{
			HIT <- which(RT1>=(RTs[d]-scanrate/2) & RT1<=(RTs[d]+scanrate/2))
			if (length(HIT)!=0)
			{
				if (length(HIT)>1)
				{
					INTmat[1,d] <- max(INT1[HIT])
				} else
				{
					INTmat[1,d] <- INT1[HIT]
				}
			}
			HIT <- which(RT2>=(RTs[d]-scanrate/2) & RT2<=(RTs[d]+scanrate/2))
			if (length(HIT)!=0)
			{
				if (length(HIT)>1)
				{
					INTmat[2,d] <- max(INT2[HIT])
				} else
				{
					INTmat[2,d] <- INT2[HIT]
				}
			}
		}

		# Fill in NA values which are between measured values with average of two neighbors

		Nas <- which(is.na(INTmat[1,]))
		REM <- c(which(Nas==1),which(Nas==length(RTs)))
		if (length(REM)!=0)
		{
			Nas <- Nas[-c(REM)]
		}
		if (length(Nas)!=0)
		{
			for (d in 1:length(Nas))
			{
				if (!(is.na(INTmat[1,Nas[d]-1])) | !(is.na(INTmat[1,Nas[d]+1])))
				{
					INTmat[1,Nas[d]] <- mean(c(INTmat[1,Nas[d]-1],INTmat[1,Nas[d]+1]),na.rm=TRUE)
				} else
				{
					INTmat[1,Nas[d]] <- 0
				}
			}
		}

		Nas <- which(is.na(INTmat[2,]))
		REM <- c(which(Nas==1),which(Nas==length(RTs)))
		if (length(REM)!=0)
		{
			Nas <- Nas[-c(REM)]
		}

		if (length(Nas!=1))
		{
			for (d in 1:length(Nas))
			{
				if (!(is.na(INTmat[2,Nas[d]-1])) | !(is.na(INTmat[2,Nas[d]+1])))
				{
					INTmat[2,Nas[d]] <- mean(c(INTmat[2,Nas[d]-1],INTmat[2,Nas[d]+1]),na.rm=TRUE)
				} else
				{
					INTmat[2,Nas[d]] <- 0
				}
			}
		}

		INTmat[is.na(INTmat)]<-0
		RES <- cor.test(INTmat[1,],INTmat[2,])
		c(RES$estimate,RES$p.value)
	} else
	{
		c(0,1)
	}
}
