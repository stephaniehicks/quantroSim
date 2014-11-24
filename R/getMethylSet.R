#' @title Export the simulateMethObject to a MethylSet to be used in the 
#' R/Bioconductor package \code{minfi}.
#'
#' @description This function takes in the simulateMethObject and 
#' returns a MethylSet object which can be used in the minfi package. 
#' Can use functions from mini such as \code{getMeth}, \code{getUnmeth}, 
#' \code{getBeta}, \code{getM}, \code{getCN}. Because it is simulated data, 
#' there is no manifest and no preprocessing in the MethylSet returned. 
#' 
#' @param simulateMethObject Must be an object created from \code{simulateMeth} 
#'
#' @return A \code{MethylSet} object
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' methTruth <- simulateMethTruth(nProbes = 2e4, nGroups = 2)
#' sim <- simulateMeth(methTruth, meth.platform = "methArrays", 
#'                    nSamps = 5, nMol = 10^6)
#' mset <- getMethylSet(sim)
getMethylSet <- function(simulateMethObject)
{
  if(!exists("objectType", where = simulateMethObject)){
    stop("Must supply a simulateMethObject created from simulateMeth().")
  }  

  if(!(simulateMethObject$typePlatform %in% c("methArrays"))){
    stop("Platform not available. Must specify a platform from 
         list.meth.platforms().")
  }
  rownames(simulateMethObject$pd) <- colnames(simulateMethObject$meth)
  output <- MethylSet(Meth = simulateMethObject$meth, 
                      Unmeth = simulateMethObject$unmeth, 
                      phenoData = AnnotatedDataFrame(simulateMethObject$pd))
  return(output)
}

