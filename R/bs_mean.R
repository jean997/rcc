#' Parametric bootstrapped means
#' @description This function implements the method described in Simon and Simon (2013)
#' @param z Parameter estimates
#' @param n Number of bootstrap replications
#' @param use.abs Rank z by absolute value.
#' @return A vector of mean estimates of the same length as \code{z}
#'@export
bs_mean <- function(z, n=200, use.abs=TRUE){
  p <- length(z)

  B <- replicate(n = n, expr = {
    w <- rnorm(p, mean=z)
    if(use.abs){
      k <- order(abs(w))
      sign(w[k])*(w[k]-z[k])
    }else{
      k <- order(w)
      w[k]-z[k]
    }
  })
  bias <- rowMeans(B)
  if(use.abs){
    j <- order(abs(z))
    k <- match(z, z[j])
    return(z - sign(z)*bias[k])
  }else{
    j <- order(z)
    k <- match(z, z[j])
    return(z - bias[k])
  }
}
