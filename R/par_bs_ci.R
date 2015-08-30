

#' Parametric bootstrapped confidence intervals to control RCC
#' @description This function implements algorithm 1 - simple
#' parametric bootstrap for normal estimates
#' @param z Parameter estimates
#' @param theta Used to produce oracle intervals for simulations whien E[z] is known.
#' @param alpha Confidence level
#' @param n Number of bootstrap replications
#' @param use.abs Rank z by absolute value.
#' @return A p by 2 matrix giving confidence intervals for each element of \code{z}
#'@export
par_bs_ci <- function(z, theta=z, alpha=0.1,
                      n=2000, use.abs=FALSE){
  p <- length(z)
  B <- replicate(n = n, expr = {
    w <- rnorm(p, mean=theta, sd=1)
    if(use.abs){
      k <- order(abs(w))
      sign(w[k])*(w[k]-theta[k])
    }else{
      k <- order(w)
      w[k]-theta[k]
    }})
  q1 <- alpha/2
  q2 <- 1-(alpha/2)
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})

  if(use.abs){
    j <- order(abs(z))
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
    which.neg <- which(z[j] < 0)
    my.ci[ which.neg , ] <- cbind(z[j][which.neg] + qs[1,which.neg], z[j][which.neg]+qs[2,which.neg])
  }else{
    j <- order(z)
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
  }
  jinv <- match(z, z[j])
  return(my.ci[jinv, ])
}
