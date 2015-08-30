
#' Parametric bootstrapped confidence intervals to control RCC
#' @description This function implements algorithm 4 - simple
#' parametric bootstrap for normal estimates with debiased mean estimation.
#' Means are debiased using bootstrapping as described in Simon and Simon (2013)
#' @param z Parameter estimates
#' @param alpha Confidence level
#' @param n Number of bootstrap replications
#' @param use.abs Rank z by absolute value.
#' @param lazy Assume the same bias in bootstrap samples as estimated in the original sample.
#' @return A p by 2 matrix giving confidence intervals for each element of \code{z}
#'@export
par_bs_db_ci <- function(z, alpha=0.1, n=2000, use.abs=FALSE, lazy=FALSE){
  p <- length(z)

  theta_tilde <- bs_mean(z, use.abs=use.abs)

  if(use.abs){
    j <- order(abs(z))
    bias <- (sign(z)*(z - theta_tilde))[j]
  }else{
    j <- order(z)
    bias <- (z-theta_tilde)[j]
  }
  B <- replicate(n = n, expr = {
    w <- rnorm(p, mean=theta_tilde, sd=1)
    if(use.abs){
      k <- order(abs(w))
      if(lazy) wadj <- w[k]-(sign(w[k])*bias)
      else wadj <- (bs_mean(z=w, n=200, use.abs=use.abs))[k]
      sign(w[k])*(wadj-theta_tilde[k])
    }else{
      k <- order(w)
      if(lazy) wadj <- w[k]-bias
      else wadj <- (bs_mean(z=w, n=200, use.abs=use.abs))[k]
      wadj-theta_tilde[k]
    }})
  q1 <- alpha/2
  q2 <- 1-(alpha/2)
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(0.05, 0.95))})
  if(use.abs){
    my.ci <- cbind(theta_tilde[j]-qs[2,], theta_tilde[j]-qs[1,])
    which.neg <- which(z[j] < 0)
    my.ci[ which.neg , ] <- cbind(theta_tilde[j][which.neg] + qs[1,which.neg], theta_tilde[j][which.neg]+qs[2,which.neg])
  }else{
    my.ci <- cbind(theta_tilde[j]-qs[2,], theta_tilde[j]-qs[1,])
  }
  jinv <- match(z, z[j])
  return(my.ci[jinv,])
}

