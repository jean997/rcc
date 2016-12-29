#' Calculate expected coverage of several alternative intervals for
#'  marginal coefficients in biomarker cutpoint example
#' @description This function can be used to generate simulations in section 3.1
#' @param n Number of samples
#' @param n.rep Number of replications
#' @param n.cutpoints Number of cutpoints between 0 and 0.0
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
biomarker_sim <- function(n=500, n.rep=1,  n.cutpoints=100, seed=NULL){
  stopifnot(n %%2 == 0)
  mean_outcome <- function(x, trt){
    if(x < 0.5) return(0)
    return((x-0.5)*trt)
  }
  
  if(!is.null(seed)) set.seed(seed)
  i <- 1
  simnames <- c("nonpar", "par", "wfb", "wfb2", "naive", "selInf1", "ash")
  COVERAGE <- array(dim=c(length(simnames), n.cutpoints, n.rep))
  WIDTH <- array(dim=c(length(simnames),  n.cutpoints, n.rep))
  
  while(i <= n.rep){
    cat(i, "..")
    #Generate data
    dat <- data.frame("trt"=rep(c(0, 1), each=n/2), "x"=runif(n=n))
    dat$mu <- apply(dat, MARGIN=1, FUN=function(z){mean_outcome(z[2],z[1])})
    dat$y <- rnorm(n=n, mean=dat$mu, sd=0.5)
    cutpoints <- seq(0, 0.9, length.out=n.cutpoints)
    stats <- lapply(cutpoints, FUN=function(xx){
      f <- lm(y~trt, data=dat[dat$x >= xx,])
      summary(f)$coefficients[2, 1:3]
    })
    stats <- data.frame(matrix(unlist(stats), byrow=TRUE, ncol=3))
    names(stats) <- c("beta", "se", "tstat")
    stats$truth <- ((1/8)-pmax(0, ((cutpoints-0.5)^2)/8))/(1-cutpoints)
    j <- order(abs(stats$tstat))
    
    #Non parametric bootstrap
    ci.nonpar <- biomarker_nonpar_ci(dat, stats, cutpoints, n.rep=500)
    COVERAGE[which(simnames=="nonpar"), ,i] <- (ci.nonpar[,1] < stats$truth & stats$truth < ci.nonpar[,2])[j]
    WIDTH[which(simnames=="nonpar"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[j]
    
    #Parametric bootstrap
    ci.par <- lr_bs_par_ci(stats$beta, stats$se, n.rep=500)
    COVERAGE[which(simnames=="par"), ,i] <- (ci.par[,1] < stats$truth & stats$truth < ci.par[,2])[j]
    WIDTH[which(simnames=="par"), , i] <- (ci.par[,2]-ci.par[,1])[j]
    
    #Parametric bootstrap with oracle mean (no covariance)
    ci.opar <- rcc:::lr_bs_par_ci(stats$beta, stats$se, n.rep=500, beta_hat_mean=stats$truth)
    COVERAGE[which(simnames=="opar"), ,i] <- (ci.opar[,1] < stats$truth & stats$truth < ci.opar[,2])[j]
    WIDTH[which(simnames=="opar"), , i] <- (ci.opar[,2]-ci.opar[,1])[j]
    
    #WFB CIs conditional on taking the top 10%
    ct <- quantile(abs(stats$tstat), probs=0.9)
    wfb <- lapply(stats$tstat, FUN=function(x){
      if(abs(x) < ct) return(c(NA, NA))
      ci <- try(rcc:::Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
      return(unlist(ci))
    })
    ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=n.cutpoints)
    ci.wfb[,1]<- ci.wfb[,1]*stats$se
    ci.wfb[,2]<- ci.wfb[,2]*stats$se
    COVERAGE[which(simnames=="wfb"), ,i] <- (ci.wfb[,1] < stats$truth & stats$truth < ci.wfb[,2])[j]
    WIDTH[which(simnames=="wfb"), , i] <- (ci.wfb[,2]-ci.wfb[,1])[j]
    
    #WFB CIs moving threshold
    wfb.2 <- lapply(stats$tstat, FUN=function(x, z){
      u <- abs(z[abs(z) < abs(x)])
      if(length(u) > 0) ct <-max(u)
      else ct <- min(abs(z))/2
      ci <- try(rcc:::Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
      return(unlist(ci))
    }, z=stats$tstat)
    ci.wfb.2 <- matrix(unlist(wfb.2), byrow=TRUE, nrow=n.cutpoints)
    ci.wfb.2[,1] <- ci.wfb.2[,1]*stats$se
    ci.wfb.2[,2] <- ci.wfb.2[,2]*stats$se
    COVERAGE[which(simnames=="wfb2"), ,i] <- (ci.wfb.2[,1] < stats$truth & stats$truth < ci.wfb.2[,2])[j]
    WIDTH[which(simnames=="wfb2"), , i] <- (ci.wfb.2[,2]-ci.wfb.2[,1])[j]
    
    #Naive ci
    ci.naive <- cbind(stats$beta-stats$se*qnorm(0.95), stats$beta + stats$se*qnorm(0.95))
    COVERAGE[which(simnames=="naive"), ,i] <- (ci.naive[,1] < stats$truth & stats$truth< ci.naive[,2])[j]
    WIDTH[which(simnames=="naive"), , i] <- (ci.naive[,2]-ci.naive[,1])[j]
    
    #Reid, Taylor, Tibshirani method (selectiveInference)
    M <- manyMeans(y=stats$tstat, k=0.1*n.cutpoints, alpha=0.1, sigma=1)
    ci.rtt1 <- matrix(nrow=n.cutpoints, ncol=2)
    ci.rtt1[M$selected.set, ] <- M$ci
    ci.rtt1[,1]<- ci.rtt1[,1]*stats$se
    ci.rtt1[,2]<- ci.rtt1[,2]*stats$se
    COVERAGE[simnames == "selInf1", ,i]<- (ci.rtt1[,1] < stats$trut & stats$trut < ci.rtt1[,2])[j]
    WIDTH[simnames=="selInf1", , i] <- (ci.rtt1[, 2] - ci.rtt1[,1])[j]
    
    ash.res <- ash(betahat = stats$beta, sebetahat = stats$se, mixcompdist = "normal")
    ci.ash <- ashci(ash.res, level=0.9, betaindex = 1:n.cutpoints )
    COVERAGE[simnames == "ash", ,i]<- (ci.ash[,1] < stats$trut & stats$trut < ci.ash[,2])[j]
    WIDTH[simnames=="ash", , i] <- (ci.ash[, 2] - ci.ash[,1])[j]
    
    
    i <- i+1
    
  }
  
  cat("\n")
  return(list("COVERAGE"=COVERAGE, "WIDTH"=WIDTH, "simnames"=simnames))
  
  
}

