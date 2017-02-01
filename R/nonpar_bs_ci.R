#' Non-Parametric bootstrapped confidence intervals to control RCC
#' @param data An n by K matrix of data
#' @param analysis.func A function that performs the analysis. It must take exactly one argument (data) and return a
#' list or data frame including an item called 'estimate' and an item calle 'statistic'. 
#' These should both be vectors of length p, the number of parameters.
#' @param rank.func A Function that takes as first argument statistic as returned by analysis.func, as second argument use.abs
#'  and returns a list with items order and rank.
#' @param level Confidence level
#' @param n.rep Number of bootstrap replications
#' @param use.abs Rank based on abs(z) rather than z
#' @param ... Additional parameters to pass to rank.func
#' @return A p by 2 matrix giving confidence intervals 
#'@export
nonpar_bs_ci <- function(data, analysis.func, rank.func=NULL, level=0.9, res.orig=NULL,
                         n.rep=1000, use.abs=TRUE, parallel=FALSE, ...){
  n <- dim(data)[1]
  dots <- list(...)
  ndots <- length(dots)
  
  #Ranking Function
  if(is.null(rank.func)){
    my.rank.func <- function(stats){rcc:::basic_rank(stats, use.abs=use.abs)}
  }else if(ndots > 0){
    my.rank.func <- function(stats){rank.func(stats, use.abs=use.abs, ...)}
  }else{
    my.rank.func <- function(stats){rank.func(stats, use.abs=use.abs)}
  }
  
  #Analysis with original data
  if(is.null(res.orig)) res.orig <- analysis.func(data)
  rank.orig <- my.rank.func(res.orig$statistic)
  p <- length(res.orig$estimate)
  stopifnot(length(res.orig$statistic)==p)
  
  
  #Bootstrapping Function
  samples <- replicate(n=n.rep, expr={sample(1:n, size=n, replace=TRUE)})
  bs_func <- function(i){
    data.new <- data[samples[,i],]
    res <- analysis.func(data.new)
    rank.new <- my.rank.func(res$statistic)
    if(use.abs){
      s <- sign(res$estimate[rank.new$order])
      s[res$estimate[rank.new$order] ==0] <- 1
    }else{ 
      s <- rep(1, p)
    }
    return(s*(res$estimate[rank.new$order] - res.orig$estimate[rank.new$order]))
  }
  
  #Do the bootstrap
  save(samples, bs_func, data, analysis.func, my.rank.func, file="temp.RData")
  if(!parallel){
    B <- sapply(1:n.rep, FUN=function(i){
      bs_func(i)
    })
  }else{
    cores <- parallel::detectCores()-2
    cl <- parallel::makeCluster(cores, type="FORK")
    on.exit(parallel::stopCluster(cl))
    B <- parSapply(cl, 1:n.rep, FUN=function(i){
      bs_func(i)
    })
  }
  #Quantile
  q1 <- (1-level)/2
  q2 <- 1-(1-level)/2
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})
 
  #Pivot
  my.ci <- cbind(res.orig$estimate[rank.orig$order]-qs[2,], res.orig$estimate[rank.orig$order]-qs[1,])
  #Reflect intervals for negative estimates if using and abs. value ranking
  if(use.abs){
    s <- sign(res.orig$estimate[rank.orig$order])
    s[res.orig$estimate[rank.orig$order]==0] <- 1
    which.neg <- which(s==-1)
    my.ci[ which.neg , ] <- cbind(res.orig$estimate[rank.orig$order][which.neg] + qs[1,which.neg],
                                res.orig$estimate[rank.orig$order][which.neg]+qs[2,which.neg])
  }else{
    s <- rep(1, p)
  }
  
  #Re-order the confidence intervals to the orignial order. Note we only get intervals for ranked parameters
  ci <- matrix(NA, nrow=p, ncol=2)
  jinv <- rank.orig$rank[!is.na(rank.orig$rank)]
  ci[!is.na(rank.orig$rank),] <- my.ci[jinv,]
  
  #While were at it, calculate a debiased mean
  my.mean <- res.orig$estimate[rank.orig$order] - s*rowMeans(B)
  meanest <- rep(NA, p)
  meanest[!is.na(rank.orig$rank)] <- my.mean[jinv]
  return(list("ci"=ci, "mean"=meanest, "rank"=rank.orig$rank))
}