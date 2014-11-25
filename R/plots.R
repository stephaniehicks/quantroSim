#' @title Plot the \code{simulateGExTruthObject} from \code{simulateGExTruth}
#'
#' @description This function plots the simulated RNA molecules on the
#' log 2 scale from the \code{simulateGExTruthObject}. 
#' 
#' @param simulateGExTruthObject Must be an object created from 
#' \code{simulateGExTruth}. 
#'
#' @return The distribution of simulated RNA molecules for each gene in 
#' each group.
#' 
#' @import Biobase minfi quantro
#' 
#' @export
#' @examples
#' geneTruth <- simulateGExTruth(nGenes = 2e4, nGroups = 2, 
#'                               pDiff = 0.05, foldDiff = 5)
#' plotGExTruth(geneTruth)
plotGExTruth <- function(simulateGExTruthObject)
{ 
  if(!exists("objectType", where = simulateGExTruthObject)){
    stop("[quantroSim]: Must supply a gene expression object created 
         from simulateGExTruth().")
  }
  
  nGenes <- simulateGExTruthObject$nGenes
  nGroups <- simulateGExTruthObject$nGroups
  pDiff <- simulateGExTruthObject$pDiff
  geneRange <- simulateGExTruthObject$geneRange
  scaleCounts <- simulateGExTruthObject$scaleCounts 

  par(mfrow=c(1,1))
  percentZero <- sapply(1:nGroups, function(x){
    round((length(which(geneRange[,x] == 0)) / nGenes )* 100) })
  if(scaleCounts != 0){ geneRangeScaled <- geneRange + scaleCounts }  
  
  matdensity(log2(geneRangeScaled), col = 1 + 1:nGroups, 
             xlab = "Simulated RNA molecules on log2 scale", ylab = "density")
  legend('top', colnames(geneRange), lty = 1, col= 1 + 1:nGroups, bty ="n")  
  
  
}


#' @title Plot the simulated gene expression arrays or RNA-Sequencing
#'
#' @description Plot the simulated gene expression arrays or RNA-Sequencing. 
#' 
#' @description This function requires the simulateGExObject from 
#' \code{simulateGEx}. Boxplots and density plots are plotted showing the 
#' distribution of gene expression. 
#'
#' @param simulateGExObject Must be an object created from \code{simulateGEx}
#'
#' @author Stephanie Hicks
#' @export
#' @examples
#' geneTruth <- simulateGExTruth(nGenes = 2e4, nGroups = 2, 
#'                              pDiff = 0.05, foldDiff = 5)
#' sim <- simulateGEx(geneTruth, GEx.platform = "GExArrays", nSamps = 5)
#' plotGEx(sim)
plotGEx <- function(simulateGExObject)
{ 
  
  if(exists("objectType", where = simulateGExObject)){	
    if( !(simulateGExObject$typePlatform %in% c("GExArrays")) ){
      stop("Platform currently not supported. Must be a platform 
           list.GEx.platforms().")
    }
    nGenes <- simulateGExObject$nGenes
    scaleCounts <- simulateGExObject$scaleCounts
    nProbes <- simulateGExObject$nProbes
    nGroups <- simulateGExObject$nGroups
    PM <- simulateGExObject$PM
    pd <- simulateGExObject$pd
    gID <- simulateGExObject$params$gID
    }
  
  groupNames <- levels(simulateGExObject$pd$Group)
  
  par(mfrow=c(1,2))
  matboxplot(log2(PM), groupFactor = simulateGExObject$pd$Group, 
             main = "log2(PM Values)", col = 1 + 1:nGroups, range = 0, las=3)
  
  matdensity(log2(PM), groupFactor = simulateGExObject$pd$Group, 
             col = 1 + 1:nGroups, xlab = "log2(PM Values)", ylab = "density")
  legend('top', groupNames, lty = 1, col= 1 + 1:nGroups, bty ="n")  
  
  
  
}  




#' @title Plot the \code{simulateMethTruthObject} from \code{simulateMethTruth}
#'
#' @description This function plots the simulate DNA meylation data from the
#' \code{simulateMethTruthObject}.   
#'
#' @param simulateMethTruthObject Must be an object created from 
#' \code{simulateMethTruth}. 
#' 
#' @return Two plots: One is a histogram of mixture model of three normal 
#' distributions. The second is the density of the transformed methylation data
#' scaled between 0 and 1. 
#'
#' @export
#' @examples
#' methTruth <- simulateMethTruth(nProbes = 2e4, nGroups = 2, 
#'                                pDiff = 0.05, pUp = 0.80)
#' plotMethTruth(methTruth)
plotMethTruth <- function(simulateMethTruthObject)
{ 
  if(!exists("objectType", where = simulateMethTruthObject)){
    stop("[quantroSim]: Must supply a DNA Methylation object created from 
         simulateMethTruth().")
  }
  
  nProbes <- simulateMethTruthObject$nProbes
  nGroups <- simulateMethTruthObject$nGroups
  pDiff <- simulateMethTruthObject$pDiff
  methRange <- simulateMethTruthObject$methRange
  
  par(mfrow=c(1,2))
  plot(density(simulateMethTruthObject$simObs$normObs), xlab = "", 
       main = "Mixture of Normal distributions")
  
  matdensity(simulateMethTruthObject$methRange, col = 1 + 1:nGroups, 
             xlab = "Simulated Methylation", ylab = "density",
             xlim = c(-0.01,1.01))
  legend('top', colnames(simulateMethTruthObject$methRange), lty = 1, 
         col= 1 + 1:nGroups, bty ="n")  
}



#' @title Plot the simulated DNA methylation arrays
#'
#' @description This function requires the simulateMethObject from 
#' \code{simulateMeth}. Boxplots and density plots are plotted showing the 
#' distribution of DNA methylation.  
#' 
#' @param simulateMethObject Must be an object created from \code{simulateMeth} 
#'
#' @author Stephanie Hicks
#' @export
#' @examples
#' methTruth <- simulateMethTruth(nProbes = 2e4, nGroups = 2, 
#'                                pDiff = 0.05, pUp = 0.80)
#' sim <- simulateMeth(methTruth, meth.platform = "methArrays", 
#'                      nSamps = 5, nMol = 10^6)
#' plotMeth(sim)
plotMeth <- function(simulateMethObject)
{ 
  
  if(exists("objectType", where = simulateMethObject)){
    if( !(simulateMethObject$typePlatform %in% c("methArrays")) ){
      stop("Platform currently not supported. Must be a platform 
           list.meth.platforms().")
    }
    meth <- simulateMethObject$meth
    unmeth <- simulateMethObject$unmeth
    pd <- simulateMethObject$pd
    nProbes <- simulateMethObject$nProbes
    nGroups <- simulateMethObject$nGroups
    gID <- simulateMethObject$params$gID
    } else { 
      methDat <- simulateMethObject 
    }
  
  groupNames <- levels(simulateMethObject$pd$Group)
  
  par(mfrow=c(2,3))
  betaObject <- getBeta(getMethylSet(simulateMethObject), offset = 100)
  
  matboxplot(betaObject, groupFactor = simulateMethObject$pd$Group, 
             main = "Observed beta values", col = 1 + 1:nGroups)
  matboxplot(log2(meth), groupFactor = simulateMethObject$pd$Group, 
             main = "Observed methylation (log2)", col = 1 + 1:nGroups)
  matboxplot(log2(unmeth), groupFactor = simulateMethObject$pd$Group, 
             main = "Observed unmethylation (log2)", col = 1 + 1:nGroups)
  
  matdensity(betaObject, groupFactor = simulateMethObject$pd$Group, 
             col = 1 + 1:nGroups, xlab = "Observed beta values",
             ylab = "density", xlim = c(-0.01,1.01))
  legend('top', groupNames, lty = 1, col= 1 + 1:nGroups, bty ="n")  
  
  matdensity(log2(meth), groupFactor = simulateMethObject$pd$Group, 
             col = 1 + 1:nGroups, xlab = "Observed methylation (log2)", 
             ylab = "density")
  legend('top', groupNames, lty = 1, col= 1 + 1:nGroups, bty ="n")  
  
  matdensity(log2(unmeth), groupFactor = simulateMethObject$pd$Group, 
             col = 1 + 1:nGroups, xlab = "Observed unmethylation (log2)", 
             ylab = "density")
  legend('top', groupNames, lty = 1, col= 1 + 1:nGroups, bty ="n")  
}


