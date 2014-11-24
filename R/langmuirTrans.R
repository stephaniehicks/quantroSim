#' @title Langmuir adsorption transformation
#'
#' @description To model DNA methylation and gene expression from arrays
#' this function represents the saturation reached in the arrays.  
#'
#' @param x expected number of methylated molecules after PCR and 
#' bisulfite-sequencing
#' @param a intensity from scanner
#' @param b scale parameter
#' @param d background noise
#'
#' @author Stephanie Hicks
#' @export
langmuirTrans <- function(x, a, b, d){
	if(length(a) == 1){ rep(a, length(x))}
	if(length(b) == 1){ rep(b, length(x))}
	if(length(d) == 1){ rep(d, length(x))}
	d + a * (x / (x + b))
}
