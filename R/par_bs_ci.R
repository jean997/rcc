

#' Parametric bootstrapped confidence intervals to control RCC
#' @description This function implements algorithm 1 - simple
#' parametric bootstrap for normal estimates
#' @param beta Parameter estimates
#' @param se Estimated standard error of z
#' @param theta Possibley shrunken estimates of 
#' @param level Confidence level
#' @param n Number of bootstrap replications
#' @param use.abs Rank based on abs(z) rather than z
#' @return A p by 2 matrix giving confidence intervals for each element of \code{z}
#'@export
par_bs_ci <- function(z, se = rep(1, length(z)), theta=z, level=0.1,
                      n=1000, use.abs=TRUE){
  p <- length(z)
  B <- replicate(n = n, expr = {
    w <- rnorm(p, mean=theta, sd=se)
    if(use.abs){
      k <- order(abs(w)/se)
      sign(w[k])*(w[k]-theta[k])
    }else{
      k <- order(w/se)
      w[k]-theta[k]
    }})

  q1 <- level/2
  q2 <- 1-(level/2)
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})

  if(use.abs){
    j <- order(abs(z)/se)
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
    which.neg <- which(z[j] < 0)
    my.ci[ which.neg , ] <- cbind(z[j][which.neg] + qs[1,which.neg], z[j][which.neg]+qs[2,which.neg])
  }else{
    j <- order(z/se)
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
  }
  jinv <- match(z, z[j])
  return(my.ci[jinv, ])
}
