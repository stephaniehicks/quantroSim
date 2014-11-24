#' @title Simulate observed DNA methylation data for a set of samples
#'
#' @description This function simulates observed DNA methylation data 
#' based on the set of platforms listed in \code{list.meth.platforms}. 
#' The DNA methylation is simulated using the Langmuir adsorption model. 
#' 
#' @param simulateMethTruthObject Must be an object created from 
#' \code{simulateMethTruth}. 
#' @param meth.platform Must specify a platform from 
#' \code{list.meth.platforms}. 
#' @param nSamps number of samples in each group
#' @param nMol number of molecules after amplification
#' @param verbose TRUE/FALSE option to print details about simulation. 
#' Default is TRUE.
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
#' @return A list of elements
#' \item{objectType}{a string specifying the type of object.}
#' \item{typePlatform}{platform used to simulate the gene expression.}
#' \item{nProbes}{number of probes}
#' \item{nSamps}{number of samples}
#' \item{nGroups}{number of groups}
#' \item{params}{list of elements reported from \code{pickLangmuirParameters}}
#' \item{pd}{phenoData containing information about the simulated samples}
#' \item{meth}{the simulated raw methylated data}
#' \item{unmeth}{the simulated raw unmethylated data}
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' methTruth <- simulateMethTruth(nProbes = 2e4, nGroups = 2, 
#'                                pDiff = 0.05, pUp = 0.80)
#' sim <- simulateMeth(methTruth, meth.platform = "methArrays", 
#'                     nSamps = 5, nMol = 10^6)
simulateMeth <- function(simulateMethTruthObject, meth.platform, nSamps, 
                          nMol = NULL, verbose = TRUE, mua = NULL, 
                          siga = NULL, mub = NULL, sigb = NULL,
                          muOpt = NULL, sigOpt = NULL, muBG = NULL, 
                          sigBG = NULL, muERR = NULL, sigERR = NULL)                         
{

	if(!exists("objectType", where = simulateMethTruthObject)){
		stop("Must supply a simulateMethTruthObject created from 
         simulateMethTruth().")
	}

  if( !(meth.platform %in% c("methArrays") ) ){
	  stop("Platform not available. Must specify a platform from 
         list.meth.platforms().")
	}
	
	if(simulateMethTruthObject$objectType == "simulateMethTruthObject"){
		rangeObject <- simulateMethTruthObject$methRange
		if(is.null(nMol)){ stop("Must provide nMol to simulate DNA methylation.") }
	} else { stop("The objectType is not supported. Must be an object from 
                simulateMethTruth().") }
	
	nProbes <- simulateMethTruthObject$nProbes
	nGroups <- simulateMethTruthObject$nGroups
	totSamps = nSamps * nGroups
  
  # DNA Methylation arrays
	if(meth.platform == "methArrays"){
    if(verbose){
      message("Simulating DNA methylation samples using the meth.platform: ", 
              meth.platform)  
    }
		ps <- pickLangmuir(objectType = simulateMethTruthObject$objectType, 
              typePlatform = meth.platform, nProbes = nProbes, 
              nSamps = nSamps, nGroups = nGroups, mua = mua, siga = siga, 
              mub = mub, sigb = sigb, muOpt = muOpt, sigOpt = sigOpt, 
              muBG = muBG, sigBG = sigBG, muERR = muERR, sigERR = sigERR)				

		## Generating new samples
		meth <- sapply(1:totSamps, function(x){ dat = ps$opt[,x] + ps$d[,x] + 
                langmuirTrans(rangeObject[,ps$gID[x]] * nMol, a=ps$a[,x],
                          b=ps$b[,x], d=0)*ps$errMeth[,x]									
		})
		rownames(meth) <- NULL
		colnames(meth) <- sapply(1:totSamps, function(x){
                          paste0("Sample_", formatC(x, flag=0, width=3))
                      })

		unmeth <- sapply(1:totSamps, function(x){ dat = ps$opt[,x] + ps$d[,x] + 
                langmuirTrans((1 - rangeObject[,ps$gID[x]]) * nMol, a=ps$a[,x], 
                           b=ps$b[,x], d=0) * ps$errUnmeth[,x]
		})
		rownames(unmeth) <- NULL
		colnames(unmeth) <- colnames(meth)

		pd <- data.frame("Sample_Name" = colnames(meth),
					           "Group" = sapply(ps$gID, function(x) 
                                paste0("Group_", formatC(x, flag=0, width=3))))

		output = list("objectType" = "simulateMethObject", 
                      "typePlatform" = meth.platform, "nProbes" = nProbes, 
                      "nSamps" = nSamps, "nGroups" = nGroups, 
                      "params" = ps, "pd" = pd, "meth" = round(meth, 0), 
                      "unmeth" = round(unmeth,0))
	}
	return(output)
}

