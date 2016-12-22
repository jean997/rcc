

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
par_bs_ci_sm <- function(z, theta=z, alpha=0.1,
                      n=200, use.abs=FALSE){
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
  #qs[1,] <- isoreg(x=1:p, y=qs[1,])$yf
  #qs[2, ]<- isoreg(x=1:p, y=qs[2,])$yf
  #qs[1,] <- smooth.spline(x=1:p, y=qs[1,])$y
  #qs[2,] <- smooth.spline(x=1:p, y=qs[2,])$y
  if(use.abs){
    j <- order(abs(z))
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
    which.neg <- which(z[j] < 0)
    n.neg <- length(which.neg)
    which.pos <- which(z[j] >=0)
    n.pos <- length(which.pos)
    my.ci[ which.neg , ] <- cbind(z[j][which.neg] +
                                    qs[1,which.neg], z[j][which.neg]+qs[2,which.neg])
    for(k in 1:2){
      my.ci[which.pos, k] <- isoreg(x=1:n.pos, y=my.ci[which.pos,k])$yf
      my.ci[which.neg, k] <- -1*isoreg(x=1:n.neg, y=-my.ci[which.neg, k])$yf
    }
  }else{
    j <- order(z)
    my.ci <- cbind(z[j]-qs[2,], z[j]-qs[1,])
    my.ci[,1] <- isoreg(x=1:p, y=my.ci[,1])$yf
    my.ci[,2] <- isoreg(x=1:p, y=my.ci[,2])$yf
  }

  jinv <- match(z, z[j])
  return(my.ci[jinv, ])
}
