unsafeFindInterval2 <- function (x, vec, rightmost.closed = FALSE, all.inside = FALSE) 
{
    nx <- as.integer(length(x))
    nv <- as.integer(length(vec))
    index <- integer(nx)
    .C("find_interv_vec", xt = vec, n = nv, x = as.double(x), 
        nx = nx, as.logical(rightmost.closed), as.logical(all.inside), 
        index, DUP = FALSE, NAOK = TRUE, PACKAGE = "base")
    index
}

#unsafeFindInterval3 <- function (x, vec, rightmost.closed = FALSE, all.inside = FALSE) {
#    .Internal(findInterval(vec, x, rightmost.closed, all.inside))
#}

unsafeFindInterval3 <- function (x, vec, rightmost.closed = FALSE, all.inside = FALSE, left.open = FALSE) {
            .Internal(findInterval(vec=vec, x=x,
                                   rightmost.closed=rightmost.closed,
                                   all.inside=all.inside, left.open = left.open))

    }

PeakML.Methods.getRawMat <- function (allRawPeaks, scan_start,scan_finis, mz_start, mz_finis,correctedRT,uncorrectedRT)
{
    scans <- scan_start:scan_finis
    Rawpeaks <- allRawPeaks[scans]
	rawRT <- uncorrectedRT[scans]
	cRT <- correctedRT[scans]

	peakExtract <- function (ind)
	{
		MZpeakset <- Rawpeaks[[ind]]
        massData = MZpeakset[,1]
        if ( version$major == '2' ) {
            idxes = unsafeFindInterval2(c(mz_start, mz_finis), massData)
        } else {
            idxes = unsafeFindInterval3(c(mz_start, mz_finis), massData)
        }
        low = idxes[1] + (idxes[1] == 0 || mz_start != massData[idxes[1]])
        high = idxes[2]
        hit = low:high

        if (high - low < 0) {
            out <- NULL
        } else {
            dout <- MZpeakset[hit,]
			dout <- rbind(dout,NULL)
            maxIntens = which.max(dout[,2])
            out <- c(rawRT[ind],cRT[ind],scans[ind],dout[maxIntens,1],dout[maxIntens,2])
        }
		out
	}

	RawMat <- do.call(rbind,lapply (1:length(Rawpeaks),peakExtract))
	if (is.null(RawMat))
	{
		RawMat <- c(1,1,1,1,1)
	}
	RawMat <- rbind(RawMat,NULL)
	if (!is.null(RawMat))
	{
		colnames (RawMat) <- c("rawRT","correctedRT","scan_id","mz","intensity")
	}
	RawMat
}
