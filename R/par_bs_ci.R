

#' Parametric bootstrapped confidence intervals to control RCC
#' @description This function implements algorithm 1 - simple
#' parametric bootstrap for normal estimates
#' @param beta Parameter estimates
#' @param se Estimated standard error of beta. Defaults to 1.
#' @param rank.func A Function that takes as first argument beta/se and returns a list with items order and rank.
#' @param theta Possibly shrunken estimates of E[beta]. Defaults to beta.
#' @param level Confidence level
#' @param n.rep Number of bootstrap replications
#' @param use.abs Base the rank on abs(beta) rather than beta
#' @param ... Additional parameters to pass to rank.func
#' @return A p by 2 matrix giving confidence intervals for each element of \code{beta}
#'@export
par_bs_ci <- function(beta, se = rep(1, length(beta)), 
                      rank.func=NULL, theta=beta, level=0.9,
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
  p <- length(beta)
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

  j <- my.rank.func(beta/se)
  my.ci <- cbind(beta[j$order]-qs[2,], beta[j$order]-qs[1,])
  if(use.abs){
    s <- sign(beta[j$order])
    s[beta[j$order]==0] <- 1
    which.neg <- which(s==-1)
    my.ci[ which.neg , ] <- cbind(beta[j$order][which.neg] + qs[1,which.neg], 
                                  beta[j$order][which.neg]+qs[2,which.neg])
  }else{
    s <- rep(1, p)
  }
  
  
  
  ci <- matrix(NA, nrow=p, ncol=2)
  jinv <- j$rank[!is.na(j$rank)]
  ci[!is.na(j$rank),] <- my.ci[jinv,]
  
  my.mean <- beta[j$order] - s*rowMeans(B)
  meanest <- rep(NA, p)
  meanest[!is.na(j$rank)] <- my.mean[jinv]
  return(list("ci"=ci, "mean"=meanest, "rank"=j$rank))
  
  
  return(ci)
}

basic_rank <- function(stats, use.abs=TRUE){
  if(use.abs) stats <- abs(stats)
  p <- length(stats)
  j <- order(stats, decreasing=TRUE)
  rank <- match(1:p, j)
  return(list("order"=j, "rank"=rank))
}