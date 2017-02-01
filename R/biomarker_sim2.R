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
biomarker_sim2 <- function(effects, n=400, n.rep=1,  nmarkers=1000, 
                          seed=NULL, parallel=TRUE){
  stopifnot(n %%2 == 0)
  if(!is.null(seed)) set.seed(seed)
  
  simnames <- c("nonpar", "par", "wfb", "naive", "selInf1", "ash")
  
  #effects <- rnorm(n=nmarkers, sd=3)
  analysis.func <- function(data){
    y <- data[,2]
    trt <- data[,1]
    stats <- apply(data[, 3:(nmarkers+2)], MARGIN=2, FUN=function(ii){
      f <- lm(y~trt*ii)
      summary(f)$coefficients[4, 1:3]
    })
    stats <- data.frame(matrix(unlist(stats), byrow=TRUE, ncol=3))
    names(stats) <- c("estimate", "se", "statistic")
    return(stats)
  }
  
  COVERAGE <- array(dim=c(length(simnames), nmarkers, n.rep))
  WIDTH <- array(dim=c(length(simnames),  nmarkers, n.rep))
  i <- 1
  trt <- rep(c(0, 1), each=n/2)
  while(i <= n.rep){
    cat(i, "..")
    #Generate data
    w <- matrix(rnorm(n=n*nmarkers), nrow=n)
    Ey <- sapply(1:n, FUN=function(i){
      sum(effects*as.numeric(w[i,] > 0))*trt[i]
    })
    y <- rnorm(n=n, mean=Ey, sd=0.5)
    dat <- data.frame("trt"=trt, "y"=y)
    
    
    
    #Stats
    ix <- apply(w, MARGIN=2, FUN=function(wj){as.numeric(wj > 0)})
    mydata <- cbind(dat, ix)
    stats <- analysis.func(mydata)
    stats$truth <- effects
    j <- order(abs(stats$statistic), decreasing = TRUE)
    
    #Non parametric bootstrap
    ci.nonpar <- nonpar_bs_ci(data=mydata, analysis.func=analysis.func, n.rep=500, parallel=parallel, level=0.9)
    COVERAGE[which(simnames=="nonpar"), ,i] <- (ci.nonpar$ci[,1] < stats$truth & stats$truth < ci.nonpar$ci[,2])[j]
    WIDTH[which(simnames=="nonpar"), , i] <- (ci.nonpar$ci[,2] -ci.nonpar$ci[,1])[j]
    
    #Parametric bootstrap
    ci.par <- par_bs_ci(stats$estimate, stats$se, n=500, level=0.9, use.abs=TRUE)
    COVERAGE[which(simnames=="par"), ,i] <- (ci.par$ci[,1] < stats$truth & stats$truth < ci.par$ci[,2])[j]
    WIDTH[which(simnames=="par"), , i] <- (ci.par$ci[,2]-ci.par$ci[,1])[j]
    
    #Parametric bootstrap with oracle mean (no covariance)
    ci.opar <- par_bs_ci(stats$estimate, stats$se, n=500, theta=stats$truth, level=0.9, use.abs=TRUE)
    COVERAGE[which(simnames=="opar"), ,i] <- (ci.opar$ci[,1] < stats$truth & stats$truth < ci.opar$ci[,2])[j]
    WIDTH[which(simnames=="opar"), , i] <- (ci.opar$ci[,2]-ci.opar$ci[,1])[j]
    
    #WFB CIs conditional on taking the top 10%
    ct <- quantile(abs(stats$statistic), probs=0.9)
    wfb <- lapply(stats$statistic, FUN=function(x){
      if(abs(x) < ct) return(c(NA, NA))
      ci <- try(rcc:::Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
      return(unlist(ci))
    })
    ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=nmarkers)
    ci.wfb[,1]<- ci.wfb[,1]*stats$se
    ci.wfb[,2]<- ci.wfb[,2]*stats$se
    COVERAGE[which(simnames=="wfb"), ,i] <- (ci.wfb[,1] < stats$truth & stats$truth < ci.wfb[,2])[j]
    WIDTH[which(simnames=="wfb"), , i] <- (ci.wfb[,2]-ci.wfb[,1])[j]
    
    #Naive ci
    ci.naive <- cbind(stats$estimate-stats$se*qnorm(0.95), stats$estimate + stats$se*qnorm(0.95))
    COVERAGE[which(simnames=="naive"), ,i] <- (ci.naive[,1] < stats$truth & stats$truth< ci.naive[,2])[j]
    WIDTH[which(simnames=="naive"), , i] <- (ci.naive[,2]-ci.naive[,1])[j]
    
    #Reid, Taylor, Tibshirani method (selectiveInference)
    M <- manyMeans(y=stats$statistic, k=0.1*nmarkers, alpha=0.1, sigma=1)
    ci.rtt1 <- matrix(nrow=nmarkers, ncol=2)
    ci.rtt1[M$selected.set, ] <- M$ci
    ci.rtt1[,1]<- ci.rtt1[,1]*stats$se
    ci.rtt1[,2]<- ci.rtt1[,2]*stats$se
    COVERAGE[simnames == "selInf1", ,i]<- (ci.rtt1[,1] < stats$trut & stats$trut < ci.rtt1[,2])[j]
    WIDTH[simnames=="selInf1", , i] <- (ci.rtt1[, 2] - ci.rtt1[,1])[j]
    
    ash.res <- ash(betahat = stats$estimate, sebetahat = stats$se, mixcompdist = "normal")
    ci.ash <- ashci(ash.res, level=0.9, betaindex = 1:nmarkers, trace=FALSE )
    COVERAGE[simnames == "ash", ,i]<- (ci.ash[,1] < stats$trut & stats$trut < ci.ash[,2])[j]
    WIDTH[simnames=="ash", , i] <- (ci.ash[, 2] - ci.ash[,1])[j]
    i <- i+1
  }
  cat("\n")
  return(list("COVERAGE"=COVERAGE, "WIDTH"=WIDTH, "simnames"=simnames))
}

