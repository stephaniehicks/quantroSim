#' @title Simulate true gene expression data for a set of groups
#'
#' @description This function simulates true gene expression data.  
#' The gene expression arrays are simulated using a zero inflated 
#' log-normal poisson distribution to represent the number or RNA 
#' molecules in a sample. 
#'
#' @param nGenes number of genes to simulate. 
#' @param nGroups number of groups. Default is 2. 
#' @param muHyp hyperparameters for mean of log normal. See 
#' \code{pickLogNormal} for default values. 
#' @param sigHyp hyperparameters for variance of log normal. See 
#' \code{pickLogNormal} for default values. 
#' @param zeroWeight proportion of zero-inflated counts. Default is 0.10.
#' @param pDiff percent of probes different relative to Group 1. If 
#' \code{nGroups == 1}, pDiff should be 0. If \code{nGroups > 1}, the 
#' length of pDiff should be equal to \code{nGroups - 1}. Default 
#' is \code{pDiff = 0.05}.
#' @param foldDiff fold difference relative to Group 1. If 
#' \code{nGroups == 1}, foldDiff is ignored. If \code{nGroups > 1}, the 
#' length of foldDiff should be equal to \code{nGroups - 1}. Default is 
#' \code{foldDiff = 5}.
#' @param scaleCounts To plot gene expression on log2 scale, pick a 
#' small amount to scale the counts by. Default is 1. 
#' @param verbose TRUE/FALSE argument specifying if verbose messages 
#' should be returned or not. Default is TRUE.
#'
#' @return A list with elements
#' \item{objectType}{A string specifying the type of object.}
#' \item{nGenes}{Number of genes simulated.}
#' \item{zeroWeight}{Proportion of zero-inflated counts used in simulation.}
#' \item{pDiff}{Percent of probes different relative to Group 1.}
#' \item{foldDiff}{fold difference relative to Group 1.}
#' \item{scaleCounts}{An positive real valued number specifying an additive 
#' scale for the RNA counts to be able to plot on the log2 scale.}
#' \item{genesDiffInd}{Boolean array referencing which probes are different
#' relative to Group 1.} 
#' \item{geneRange}{Data frame containing the simulated RNA molecules
#' for each gene and for each group.}
#' 
#' @author Stephanie Hicks
#' @export
#' @seealso \code{\link{simulateMethTruth}}
#' @examples
#' geneTruth <- simulateGExTruth(nGenes = 2e4, nGroups = 2, 
#'                                pDiff = 0.05, foldDiff = 5)
simulateGExTruth <- function(nGenes, nGroups = 2, muHyp = NULL, sigHyp = NULL, 
                             zeroWeight = 0.10, pDiff = 0.05, foldDiff = 5, 
                             scaleCounts = 1, verbose = TRUE)						

{
  if((nGroups > 1) && (length(pDiff) != (nGroups-1))){ 
    stop("[epigenomeSim]: pDiff must have length (nGroups - 1) if nGroups > 
         1")}

  if((nGroups == 1) && ((length(pDiff) != 1) || (pDiff != 0) )){
    warning("[epigenomeSim]: Cannot simulate differences between more than 
            one group if only one group exists. Using pDiff = 0.")
    pDiff = 0
  }
  
  if((nGroups > 1) && (length(foldDiff) != (nGroups-1))){ 
    stop("[epigenomeSim]: foldDiff must have length (nGroups - 1) if 
         nGroups > 1")}
  
  if((nGroups == 1) && (length(foldDiff) != 1)){
    stop("[epigenomeSim]: Cannot simulate differences between more than one 
          group if only one group exists. Fold diff should be length 1. ")}
    
  if( any(pDiff < 0) || any(pDiff > 1) ){
    stop("[epigenomeSim]: pDiff must be between 0 and 1.")}

  if( (zeroWeight < 0) || (zeroWeight > 1) ){
    stop("[epigenomeSim]: zeroWeight must be between 0 and 1.")}

  if(!is.null(muHyp) && (nGenes != length(muHyp))){ 
    stop("[epigenomeSim]: Hyperparameters for muHyp must be equal to number 
         of genes.")}
  
  if(!is.null(sigHyp) && (nGenes != length(sigHyp))){ 
    stop("[epigenomeSim]: Hyperparameters for sigHyp must be equal to number 
         of genes.")}

  poisMeanVals <- pickLogNormal(N = nGenes, muHyp = muHyp, sigHyp = sigHyp)

  rangeMu <- range(poisMeanVals)
  if(verbose){
      mes <- "[epigenomeSim]: Simulating RNA transcript counts using a Poisson 
          distribution with mean parameters from %.2f to %.2f"
      message(sprintf(mes, rangeMu[1], rangeMu[2]))
  }    


  # zero-inflated log-normal Poisson counts 
  zeros <- sample(c(0,1), nGenes, replace = TRUE, 
                  prob = c(zeroWeight, 1 - zeroWeight))
  poisCounts <- sapply(poisMeanVals, 
                       function(x){ rpois(1, lambda = x) })
  nRNAmol <- zeros
  nRNAmol[which(zeros == 1)] <- poisCounts[which(zeros == 1)]
  nR <- length(nRNAmol)
  
  if(nGroups == 1){
    output = list("objectType" = "simulateGExTruthObject", "nGenes" = nGenes, 
                  "nGroups" = nGroups, "zeroWeight" = zeroWeight,
                  "pDiff" = pDiff, "foldDiff" = foldDiff, 
                  "scaleCounts" = scaleCounts, "genesDiffInd" = NA, 
                  "geneRange" = as.matrix(nRNAmol))
  }

  if(nGroups > 1){
    geneRangeMat <- replicate(nGroups, nRNAmol + rpois(nR, 1))
    diffInd <- array(NA, dim = c(nR, nGroups - 1))
    for(i in seq_len(nGroups-1)){
      diffInd[,i] <- sample(c(TRUE, FALSE), nR, replace = TRUE, 
                        prob = c(pDiff[i], 1 - pDiff[i]))
      geneRangeMat[diffInd[,i],i+1] <- foldDiff[i] * nRNAmol[diffInd[,i]] 
    }
    colnames(geneRangeMat) <- sapply(1:nGroups, function(x) 
      paste0("Group_", formatC(x, flag=0, width=2)))

    output = list("objectType" = "simulateGExTruthObject", "nGenes" = nGenes,
                  "nGroups" = nGroups, "zeroWeight" = zeroWeight, 
                  "pDiff" = pDiff, "foldDiff" = foldDiff, 
                  "scaleCounts" = scaleCounts, "genesDiffInd" = diffInd, 
                  "geneRange" = geneRangeMat)
	}
	return(output)
}
