PeakML.Methods.getRTWindowFromString <- function (lb, ub){
  rtDiff <- ub - lb
  meanRT <- (ub+lb) / 2
  rtWin <- rtDiff/2
  c(meanRT, rtWin)
}
