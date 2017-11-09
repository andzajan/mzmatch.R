PeakML.Methods.baseCorrection <- function(signals, lambda=100){
	
					# Correct the baseline of noisy signals by estimating the trend based on asymmetric least squares
	# signals <- intensities
	# lambda <- smoothing parameter
 
	signals[which(signals==NA|signals==Inf|signals==NaN)]=0
        baseline <- evalWithTimeout({asysm(signals, lambda);}, timeout=10, onTimeout="silent");

        if(!is.null(baseline)){
          #baseline <- asysm(signals, lambda)
          corSignals <- signals - baseline
          corSignals <- corSignals - quantile(corSignals, probs=0.25)
          corSignals[which(corSignals < 0)]=0
        } else {
          corSignals <- signals
        }
        corSignals
}
