#' @title Simulate observed gene expression data for a set of samples
#'
#' @description This function simulates observed gene expression data 
#' based on the set of platforms listed in \code{list.GEx.platforms}. 
#' The gene expression arrays are simulated using the Langmuir adsorption 
#' model. 
#' 
#' @param simulateGExTruthObject Must be an object created from 
#' \code{simulateGExTruth}. 
#' @param GEx.platform Must specify a platform from \code{list.GEx.platforms}. 
#' @param nSamps number of samples in each group
#' @param usePCR TRUE/FALSE option. Default is FALSE. 
#' If TRUE, amplify RNA fragments with PCR.  
#' @param nPCRcycles number of PCR cycles. Default is 15. 
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
#' \item{usePCR}{A TRUE/FALSE object specifying if PCR was used to amplify 
#' the RNA counts.}
#' \item{nGenes}{number of genes}
#' \item{nProbes}{number of probes}
#' \item{nSamps}{number of samples}
#' \item{nGroups}{number of groups}
#' \item{params}{list of elements reported from \code{pickLangmuirParameters}}
#' \item{geneNames}{gene names (or probeset names)}
#' \item{nProbesPerGene}{number of probes per gene (or probeset)}
#' \item{pd}{phenoData containing information about the simulated samples}
#' \item{PM}{the simulated gene expression data}
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' geneTruth <- simulateGExTruth(nGenes = 2e4, nGroups = 2, 
#'                              pDiff = 0.05, foldDiff = 5)
#' sim <- simulateGEx(geneTruth, GEx.platform = "GExArrays", nSamps = 5)
simulateGEx <- function(simulateGExTruthObject, GEx.platform = "GExArrays", 
                        nSamps, usePCR = FALSE, nPCRcycles = 15, 
                        verbose = TRUE, mua = NULL, siga = NULL, mub = NULL, 
                        sigb = NULL, muOpt = NULL, sigOpt = NULL, muBG = NULL, 
                        sigBG = NULL, muERR = NULL, sigERR = NULL)
{
	if(!exists("objectType", where = simulateGExTruthObject)){
		stop("Must supply simulateGExTruthObject created from simulateGExTruth().")
	}

  if( !(GEx.platform %in% c("GExArrays") ) ){
	  stop("Platform not available. Must specify a platform from 
         list.GEx.platforms(). ")
	}
		
	if(simulateGExTruthObject$objectType == "simulateGExTruthObject"){
		rangeObject <- simulateGExTruthObject$geneRange
	} else { stop("The objectType is not supported. Must be an object from 
                simulateGExTruth().") } 
  
	nGenes <- simulateGExTruthObject$nGenes
	nGroups <- simulateGExTruthObject$nGroups
	totSamps <- nSamps * nGroups
	scaleCounts <- simulateGExTruthObject$scaleCounts

	# Gene Expression Arrays
	if(GEx.platform == "GExArrays"){
    if(verbose){
	    message("Simulating gene expression samples using the GEx.platform: ", 
              GEx.platform)
    }
	  # For each gene in simulateGExTruthObject, simulate probe level counts
    # number of probes representing each gene
    nProbesPerGene <- rbinom(nGenes, 20, prob = 0.6) 
	  nProbes <- sum(nProbesPerGene)
	  geneInfo <- data.frame("geneID" = 1:nGenes, rangeObject, 
                           "nProbes" = nProbesPerGene)
	  
    geneID = sapply(1:nGenes, function(x){
                    paste0("GENE", formatC(x, flag=0, width=5)) })
    totPCRfragments <- sapply(1:nGroups, function(x){ 
                              rep(rangeObject[,x], nProbesPerGene)})
	  probFragments <- runif(nProbes, 0, 1)
	  prePCRfragments <- matrix(rbinom(nProbes * nGroups, totPCRfragments,  
	                                   replicate(2, probFragments)), 
                              nProbes, nGroups)
	  if(usePCR){
	    message("Performing PCR amplification of RNA transcript counts with ", 
              nPCRcycles, " PCR cycles.")
	    probeRange <- (prePCRfragments * (1 + rbeta(nProbes, 1, 4))^(nPCRcycles))
	  } else {
	    message("No PCR amplification of RNA transcript counts.")
	    probeRange <- prePCRfragments
	  }
	  
		ps <- pickLangmuir(objectType = simulateGExTruthObject$objectType, 
              typePlatform = GEx.platform, nProbes = nProbes, nSamps = nSamps, 
              nGroups = nGroups, mua = mua, siga = siga, mub = mub, 
              sigb = sigb, muOpt = muOpt, sigOpt = sigOpt, muBG = muBG, 
              sigBG = sigBG, muERR = muERR, sigERR = sigERR)

		PM <- sapply(1:totSamps, function(x){ dat = ps$opt[,x] + ps$d[,x] + 
			        langmuirTrans(probeRange[,ps$gID[x]], a=ps$a[,x], 
                            b=ps$b[,x], d=0)*ps$errPM[,x] })

		rownames(PM) <- NULL
		colnames(PM) <- sapply(1:totSamps, function(x){ 
                        paste0("Sample_", formatC(x, flag=0, width=3)) })

		pd <- data.frame("Sample_Name" = colnames(PM), 
                     "Group" = sapply(ps$gID, function(x)
                                paste0("Group_", formatC(x, flag=0, width=3))))
		output <- list("objectType" = "simulateGExObject", 
                   "typePlatform" = GEx.platform, "usePCR" = usePCR, 
                   "nGenes" = nGenes, "nProbes" = nProbes, "nSamps" = nSamps, 
                   "nGroups" = nGroups, "params" = ps, "geneNames" = geneID, 
                   "nProbesPerGene" = nProbesPerGene, 
                   "scaleCounts" = scaleCounts, "pd" = pd, 
                   "probeRange" = probeRange, "PM" = round(PM, 0))
	}
	
	return(output)
}


