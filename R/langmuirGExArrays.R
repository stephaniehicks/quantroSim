#' @title Langmuir parameters for gene expression using arrays
#'
#' @description These are the parameters used to simulate gene expression
#' arrays. The parmeters are used in the Langmuir adsorption model. 
#'
#' @param nProbes number of probes
#' @param nSamps number of samples in each group
#' @param nGroups number of groups
#' @param mua hyperparameter: mean of log normal distribution for a
#' @param siga hyperparameter: variance of log normal distribution for a 
#' @param mub hyperparameter: mean of log normal distribution for b
#' @param sigb hyperparameter: variance of log normal distribution for b
#' @param muOpt hyperparameter: mean of log normal distribution for optical 
#' noise
#' @param sigOpt hyperparameter: variance of log normal distribution for 
#' optical noise
#' @param muBG hyperparameter: mean of log normal distribution for d
#' @param sigBG hyperparameter: variance of log normal distribution for d
#' @param muERR hyperparameter: mean of log normal distribution for 
#' measurement error
#' @param sigERR hyperparameter: variance of log normal distribution for 
#' measurement error
#'
#' @return A list with elements
#' \item{gID}{Group number ID}
#' \item{a}{Simulated a values for the Langmuir Adsorption model}
#' \item{b}{Simulated b values for the Langmuir Adsorption model} 
#' \item{d}{Simulated d values for the Langmuir Adsorption model} 
#' \item{opt}{Simulated optical noise}
#' \item{errPM}{Simulated measurement error }
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' langmuirGExArrays(nProbes = 50000, nSamps = 4, nGroups = 2)
langmuirGExArrays <- function(nProbes, nSamps, nGroups, 
            mua = NULL, siga = NULL, mub = NULL, sigb = NULL, muOpt = NULL, 
 					  sigOpt = NULL, muBG = NULL, sigBG = NULL, muERR = NULL, sigERR = NULL)
{
	totSamps <- nSamps * nGroups

  # Default parameters
	if(is.null(siga)){ siga = 0.1 * diag(totSamps) }
	if(is.null(sigb)){ sigb = 0.1 * diag(totSamps) }
	if(is.null(sigOpt)){ sigOpt = 0.1 * diag(totSamps) }
	
	if(is.null(mua)){ mua = rep(20, totSamps) }
	if(is.null(mub)){ mub = rep(18, totSamps) }
	if(is.null(muOpt)){ muOpt = rep(5, totSamps) }
	
	if(is.null(muBG)){ muBG = rep(5, totSamps)}
	if(is.null(sigBG)){ sigBG = 0.1 * diag(totSamps) }	
	if(is.null(muERR)){ muERR = rep(0, totSamps) }
	if(is.null(sigERR)){ sigERR = 0.1 * diag(totSamps) }	
	
	# Sample parameters: a, b, d, err
	a <- sweep(2^(mvrnorm(nProbes, mu = rep(0, totSamps), 
                        Sigma = 0.01 * diag(totSamps))), 2,
	           2^(mvrnorm(n = 1, mu = mua, Sigma = siga)), FUN = "*")
	b <- sweep(2^(mvrnorm(nProbes, mu = rep(0, totSamps), 
                        Sigma = 0.01 * diag(totSamps))), 2,
	           2^(mvrnorm(n = 1, mu = mub, Sigma = sigb)), FUN = "*")	
	opt <- sweep(2^(mvrnorm(n = nProbes, mu = rep(0, totSamps), 
                          Sigma = 0.01 * diag(totSamps))), 2, 
				2^(mvrnorm(n = 1, mu = muOpt, Sigma = sigOpt)), FUN="*")
	d <- 2^(mvrnorm(n = nProbes, mu = muBG, Sigma = sigBG))
	errPM <- 2^(mvrnorm(n = nProbes, mu = muERR, Sigma = sigERR))
		
	list("gID" = rep(1:nGroups, each = nSamps), "a" = a, "b" = b, "d" = d,
	     "opt" = opt, "errPM" = errPM)	
}	

