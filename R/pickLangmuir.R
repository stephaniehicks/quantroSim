#' @title Simulate parameters from the Langmuir Adsorption model
#' 
#' @description This is a helper function that points to the
#' functions which contain the default parameters for each platform 
#' available in this package. This function can simulate 
#' parameters from the Langmuir model for each of the platforms available
#' \code{list.GEx.platforms} and \code{list.meth.platforms}. 
#' 
#' @param objectType Must be "simulateGExTruthObject" or 
#' "simulateMethTruthObject"
#' @param typePlatform Must specify a platform from \code{list.GEx.platforms}
#' or \code{list.meth.platforms}.
#' @param nProbes number of probes
#' @param nSamps number of samples in each group
#' @param nGroups number of groups
#' @param mua hyperparameter: mean of log normal distribution for a 
#' @param siga hyperparameter: variance of log normal distribution for a 
#' @param mub b hyperparameter: mean of log normal distribution for b
#' @param sigb b hyperparameter: variance of log normal distribution for b
#' @param muOpt hyperparameter: mean of log normal distribution for 
#' optical noise
#' @param sigOpt hyperparameter: variance of log normal distribution for 
#' optical noise
#' @param muBG hyperparameter: mean of log normal distribution for d
#' @param sigBG hyperparameter: variance of log normal distribution for d
#' @param muERR hyperparameter: mean of log normal distribution for 
#' measurement error
#' @param sigERR hyperparameter: variance of log normal distribution for 
#' measurement error
#'
#' @return The output of this function is the output directly reported from 
#' the functions \code{langmuirGExArrays} and \code{langmuirMethArrays}. 
#'
#' @author Stephanie Hicks
#' @export
#' @examples
#' ps <- pickLangmuir(objectType = "simulateMethTruthObject", 
#'                    typePlatform = "methArrays", nProbes = 1e5, 
#'                    nSamps = 5, nGroups = 2)
pickLangmuir <- function(objectType, typePlatform, nProbes, nSamps, nGroups, 
                  mua = NULL, siga = NULL, mub = NULL, sigb = NULL, 
                  muOpt = NULL, sigOpt = NULL, muBG = NULL, sigBG = NULL, 
                  muERR = NULL, sigERR = NULL)
{
	if(objectType == "simulateGExTruthObject"){
    if(typePlatform == "GExArrays"){ .langmuirFunction <- langmuirGExArrays }
    if(typePlatform == "RNASeq"){ .langmuirFunction <- langmuirRNASeq }
	}
	
	if(objectType == "simulateMethTruthObject"){
		if(typePlatform == "methArrays"){ .langmuirFunction <- langmuirMethArrays }
	}

  out <- .langmuirFunction(nProbes = nProbes, nSamps = nSamps, nGroups = nGroups, 
	                   mua = mua, siga = siga, mub = mub, sigb = sigb, 
	                   muOpt = muOpt, sigOpt = sigOpt, muBG = muBG, sigBG = sigBG, 
	                   muERR = muERR, sigERR = sigERR)
  return(out)
}
