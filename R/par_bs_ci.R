

#' Parametric bootstrapped confidence intervals to control RCC
#' @description This function implements algorithm 1 - simple
#' parametric bootstrap for normal estimates
#' @param beta Parameter estimates
#' @param se Estimated standard error of z
#' @param rank.func A Function that takes as first argument beta/se and returns a list with items order and rank.
#' @param theta Possibley shrunken estimates of 
#' @param level Confidence level
#' @param n.rep Number of bootstrap replications
#' @param use.abs Rank based on abs(z) rather than z
#' @param ... Additional parameters to pass to rank.func
#' @return A p by 2 matrix giving confidence intervals for each element of \code{z}
#'@export
par_bs_ci <- function(z, se = rep(1, length(z)), 
                      rank.func=NULL, theta=z, level=0.9,
                      n.rep=1000, use.abs=TRUE, ...){
  dots <- list(...)
  ndots <- length(dots)
  
  if(is.null(rank.func)){
    my.rank.func <- function(stats){basic_rank(stats, use.abs=use.abs)}
  }else if(ndots > 0){
    my.rank.func <- function(stats){rank.func(stats, use.abs=use.abs, ...)}
  }else{
    my.rank.func <- function(stats){rank.func(stats, use.abs=use.abs)}
  }
  p <- length(z)
  B <- replicate(n = n.rep, expr = {
    w <- rnorm(p, mean=theta, sd=se)
    k <-my.rank.func(w/se)
    if(use.abs){
      sign(w[k$order])*(w[k$order]-theta[k$order])
    }else{
      w[k$order]-theta[k$order]
    }
  })

  q1 <- (1-level)/2
  q2 <- 1-(1-level)/2
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})

  j <- my.rank.func(z/se)
  my.ci <- cbind(z[j$order]-qs[2,], z[j$order]-qs[1,])
  if(use.abs){
    which.neg <- which(z[j$order] < 0)
    my.ci[ which.neg , ] <- cbind(z[j$order][which.neg] + qs[1,which.neg], z[j$order][which.neg]+qs[2,which.neg])
  }
  ci <- matrix(NA, nrow=p, ncol=2)
  jinv <- j$rank[!is.na(j$rank)]
  ci[!is.na(j$rank),] <- my.ci[jinv,]
  return(ci)
}

basic_rank <- function(stats, use.abs=TRUE){
  if(use.abs) stats <- abs(stats)
  p <- length(stats)
  j <- order(stats, decreasing=TRUE)
  rank <- match(1:p, j)
  return(list("order"=j, "rank"=rank))
}