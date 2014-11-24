#' @title Simulate a mixture of normal distributions
#'
#' @description Helper function to simulate the normal mixtures with varying
#' hyperparameters (means and variances). 
#' 
#' @param N number of observations to simulate from the mixture of normals
#' @param nMix number of mixtures. Default is 3. 
#' @param propMix mixture proportions. If NULL (default), a mixture 
#' distribution with propMix = c(0.4, 0.2, 0.4) is simulated. Otherwise 
#' user must supply propMix. Length should equal nMix. 
#' @param muMix mean of nMix distributions. If NULL (default), a mixture 
#' distribution with muMix = c(-3, 1, 3) is simulated. Otherwise user must 
#' supply propMix. Length should equal nMix.
#' @param sigMix variance of nMix distributions. If NULL (default), a mixture 
#' distribution with sigMix = c(3, .4, 3) is simulated. Otherwise user must 
#' supply sigMix. Length should equal nMix.
#' 
#' @return A list with elements
#' \item{trueParams}{A data frame with the input parameters used to simulate
#' the mixture distribution (propMix, muMix and sigMix)}
#' \item{simObs}{A data frame containing trueMix (indicator specifying which 
#' distribution the observation came from) and normObs (the observed value 
#' from the mixture distribution).}
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' pickNormalMix(N = 10)
pickNormalMix <- function(N, nMix = NULL, propMix = NULL, 
                          muMix = NULL, sigMix = NULL)
{
		if(is.null(nMix)){ nMix <- 3}
		
		if(nMix == 3){
			if(is.null(propMix)){ propMix <- c(0.4, 0.2, 0.4) }
			if(is.null(muMix)){ muMix <- c(-3, 1, 3) }
			if(is.null(sigMix)){ sigMix <- c(3, .4, 3) }
		}
    
		if((nMix != 3) && (is.null(propMix) || is.null(muMix) || 
                         is.null(sigMix)) ){
		  stop("If number of mixtures (nMix) is not 3, user must supply mixing 
           proportions (propMix), mean (muMix) and variance (sigMix).")
		}
		
		if( (!is.null(propMix) && (nMix != length(propMix)) ) ||
		      (!is.null(muMix) && (nMix != length(muMix)) ) ||
		      (!is.null(sigMix) && (nMix != length(sigMix)) ) ){
		  stop("Length of propMix, muMix and sigMixtures must be same as nMix.")
		}
	
	trueParams <- data.frame("propMix" = propMix, "muMix" = muMix, 
                           "sigMix" = sigMix)
	z <- sample(1:nMix, N, replace = TRUE, prob = propMix)
	d = array(mvrnorm(N, mu = c(muMix), Sigma = diag(sigMix)), 
            dim = c(N, length(muMix)))

	normObs = array(0, dim = N)
	for(i in 1:nMix){ normObs[z == i] <- d[z == i, i] }
	simObs = data.frame("trueMix" = z, "normObs" = normObs)
	
	list("trueParams" = trueParams, "simObs" = simObs)
}
