#' @title Simulate true DNA methylation data for a set of groups
#'
#' @description This function simulates true DNA methylation data. The true 
#' DNA methylation data is simulated using a mixture of three normal 
#' distributions to represent a probe that is not methylated, semi-methylated 
#' and methylated. To scale the values between 0 and 1, the values are 
#' transformed using an inverse-logit transformation. These scaled values 
#' represent the true proportion of methylation at each probe. 
#'
#' @param nProbes number of probes
#' @param nGroups number of groups. Default is 2.
#' @param nMixtures number of normal distributions to simulate. Default is 3. 
#' Can be changed, but user must then supply mixing proportions (propMixtures), 
#' mean (muMixtures) and variance (sigMixtures).
#' @param propMixtures weights or mixture proportions for the normal 
#' distributions. Should be a vector of length nMixtures. See 
#' \code{pickNormalMix} for default values.
#' @param muMixtures mean of the normal distributions. Should be a vector of 
#' length nMixtures. See \code{pickNormalMix} for default values. 
#' @param sigMixtures variance of the three normal distributions. Should be 
#' a vector of length nMixtures.  See \code{pickNormalMix} for default values. 
#' @param pDiff percent of probes different relative to Group 1. If 
#' \code{nGroups == 1}, pDiff should be 0. If \code{nGroups > 1}, the 
#' length of pDiff should be equal to \code{nGroups - 1}. Default 
#' is \code{pDiff = 0.05}.
#' @param pUp proportion of pDiff probes that are methylated relative to 
#' Group 1. If \code{nGroups == 1}, pUp is ignored. If \code{nGroups > 1}, the 
#' length of pUp should be equal to \code{nGroups - 1}. Default is 
#' \code{pUp = 0.80}.
#' @param verbose TRUE/FALSE argument specifying if verbose messages 
#' should be returned or not. Default is TRUE.
#' 
#' @return A list with elements
#' \item{objectType}{A string specifying the type of object.}
#' \item{nProbes}{Number of probes.}
#' \item{nGroups}{Number of groups.}
#' \item{trueParams}{Add me.}
#' \item{simObs}{Add me.}
#' \item{methDiffInd}{Add me.}
#' \item{methRange}{Add me.}
#' 
#' @author Stephanie Hicks
#' @export
#' @seealso \code{\link{simulateGExTruth}}
#' @examples
#' methTruth <- simulateMethTruth(nProbes = 2e4, nGroups = 2, 
#'                                  pDiff = 0.05, pUp = 0.80)
simulateMethTruth <- function(nProbes, nGroups = 2, nMixtures = 3, 
                              propMixtures = NULL, muMixtures = NULL, 
                              sigMixtures = NULL, pDiff = 0.05, pUp = 0.80, 
                              verbose = TRUE)
{
  if((nGroups > 1) && (length(pDiff) != (nGroups-1))){ 
    stop("[quantroSim]: pDiff must have length (nGroups - 1) if nGroups > 
         1")}
  
  if((nGroups == 1) && ((length(pDiff) != 1) || (pDiff != 0) )){
    warning("[quantroSim]: Cannot simulate differences between more than 
            one group if only one group exists. Using pDiff = 0.")
    pDiff = 0
  }
    
  if((nGroups > 1) && (length(pUp) != (nGroups-1))){ 
    stop("[quantroSim]: pUp must have length (nGroups - 1) if 
         nGroups > 1")}
  
  if((nGroups == 1) && (length(pUp) != 1)){
    stop("[quantroSim]: Cannot simulate differences between more than one 
          group if only one group exists. pUp should be length 1. ")}
  
  if( any(pDiff < 0) || any(pDiff > 1) ){
    stop("[quantroSim]: pDiff must be between 0 and 1.")}
  
  if((nMixtures != 3) && (is.null(propMixtures) || 
                          is.null(muMixtures) || is.null(sigMixtures)) ){
    stop("[quantroSim]: If number of mixtures (nMixtures) is not 3, user 
         must supply mixing proportions (propMixtures), mean (muMixtures) 
         and variance (sigMixtures).")
  }
  
  if( (!is.null(propMixtures) && (nMixtures != length(propMixtures)) ) ||
        (!is.null(muMixtures) && (nMixtures != length(muMixtures)) ) ||
        (!is.null(sigMixtures) && (nMixtures != length(sigMixtures)) ) ){
    stop("[quantroSim]: Length of propMixtures, muMixtures and sigMixtures 
         must be same as nMixtures.")
  }
  
  results <- pickNormalMix(N = nProbes, nMix = nMixtures, 
                           propMix = propMixtures, muMix = muMixtures, 
                           sigMix = sigMixtures)
  prop <- results$trueParams$propMix
  mu <- results$trueParams$muMix
  sig <- results$trueParams$sigMix
  
  if(verbose){
    mes <- "[quantroSim]: Simulating a mixture of %d Normal distributions 
            with mean (%s) and standard deviation (%s)"
    message(sprintf(mes, nMixtures, paste(mu, collapse=', '), 
                    paste(sig, collapse=', ')))
  }    
  
	if(nGroups == 1){ 
		methObsMat <- 1 / (1 + exp(-(results$simObs$normObs + 
                                   rnorm(nProbes, 0, 0.01))))
		output = list("objectType" = "simulateMethTruthObject", 
                  "nProbes" = nProbes, "nGroups" = nGroups, 
                  "trueParams" = results$trueParams, "simObs" = results$simObs,
		              "methDiffInd" = diffInd, "methRange" = methObsMat)
	}

	if(nGroups > 1){
    if(all(pDiff == 0)){
      normObs <- results$simObs$normObs
      normObsMat <- replicate(nGroups, normObs + rnorm(nProbes, 0, 0.01))
      methObsMat <- 1 / (1 + exp(-normObsMat))
      colnames(methObsMat) <- sapply(1:nGroups, function(x){
        paste0("Group_", formatC(x, flag=0, width=2))})
	    
      diffInd <- array(FALSE, dim = c(nProbes, nGroups - 1))
      output = list("objectType" = "simulateMethTruthObject", 
                    "nProbes" = nProbes, "nGroups" = nGroups, 
                    "trueParams" = results$trueParams,
                    "simObs" = results$simObs, "pDiff" = pDiff, 
                    "methDiffInd" = diffInd, "methRange" = methObsMat)
    }

    if(any(pDiff != 0)){
      trueZ <- results$simObs$trueMix
      
      normObs <- results$simObs$normObs
      normObsMat <- replicate(nGroups, normObs + rnorm(nProbes, 0, 0.01))
   
      downProbes <- (trueZ == 1)
      nDown <- length(which(downProbes))
      nUp <- length(which(!downProbes))
        
      sampleDownProbes = downProbesMat = replicate(nGroups, downProbes)
      diffInd = nowUpProbes = 
        nowDownProbes = array(NA, dim = c(nProbes, nGroups - 1))
      for(i in seq_len(nGroups-1)){
        nowUpProbes[,i] <- seq_len(nProbes) %in% 
            sample(which(downProbes), 
                   size = min(nProbes*pDiff[i]*pUp[i], nDown))

        nowDownProbes[,i] <- seq_len(nProbes) %in% 
            sample(which(!downProbes), 
                   size = min(nProbes*pDiff[i]*(1-pUp[i]) + 1, nUp))
        
        diffInd[,i] <- ifelse(nowUpProbes[,i] | nowDownProbes[,i], TRUE, FALSE)
        
        downProbesMat[,i+1] <- ifelse(nowUpProbes[,i], FALSE,
                                  ifelse(nowDownProbes[,i], TRUE, downProbes)) 
        
        normObsMat[nowUpProbes[,i],i+1] <- 
          pickNormalMix(N = length(which(nowUpProbes[,i])), 
                        nMix = nMixtures-1, propMix = c(0.2, 0.80), 
                        muMix = mu[2:3], sigMix = sig[2:3])$simObs$normObs
        
        normObsMat[nowDownProbes[,i],i+1] <-
          rnorm(length(which(nowDownProbes[,i])), 
                mean = mu[1], sd = sig[1])
      }
	
			methObsMat <- 1 / (1 + exp(-normObsMat))
			colnames(methObsMat) <- sapply(1:nGroups, function(x){
			    paste0("Group_", formatC(x, flag=0, width=2))})

			output = list("objectType" = "simulateMethTruthObject", 
                    "nProbes" = nProbes, "nGroups" = nGroups, 
                    "trueParams" = results$trueParams, 
                    "simObs" = results$simObs, "pDiff" = pDiff, 
                    "methDiffInd" = diffInd, "methRange" = methObsMat)
		} 
	}
	return(output)
}
