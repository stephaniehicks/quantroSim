#' @title Simulate the log normal values
#'
#' @description Helper function to simulate the log normal mean values with 
#' varying hyperparameters (means and variances). To be used as input for 
#' the means of the Poisson distribution. 
#' 
#' @param N number of log normal means to simulate
#' @param muHyp If NULL (default), hyperparameters for mean of log normal 
#' are simulated using Uniform(0, 6). Otherwise, user must supply muHyp.  
#' @param sigHyp If NULL (default), hyperparameters for variance of log 
#' normal are simulated using Uniform(0.1, 2). Otherwise, user must 
#' supply sigHyp.
#'
#' @return The output is a vector of log normal values. 
#' 
#' @author Stephanie Hicks
#' @export
#' @examples
#' pickLogNormal(N = 1000)
pickLogNormal <- function(N, muHyp = NULL, sigHyp = NULL)
{
	if(is.null(muHyp)){ muHyp <- runif(N, 0, 6) }
	if(is.null(sigHyp)){ sigHyp <- runif(N, 0.1, 2) }	
    
	apply(cbind(muHyp, sigHyp), 1, function(x){ 2^(rnorm(1, x[1], x[2]))})
}
