
#' Calculate expected coverage of several alternative intervals for
#'  marginal coefficients in linear regression.
#' @description This function can be used to generate simulations in section 3.1
#' @param X.pop A matrix of coefficients for an entire population.
#' @param which.sample Indices of individuals who are sampled
#' @param n.rep Number of simulations to run
#' @param index Which SNPs contribute to the phenotype
#' @param beta Effect sizes for contributing SNPs.
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
linreg_sim <- function(X.pop, which.sample=1:100, n.rep=1,
                      index= c(), beta=c(0), seed=NULL, parallel=FALSE){
  if(!is.null(seed)) set.seed(seed)
	n <-length(which.sample)
  p <- dim(X.pop)[2]
  X.sample <- X.pop[which.sample, ]
  
  lr_func <- function(data){
    y <- data[,1]
    X <- data[, -1]
    ests <- rcc:::many_lr(y, X, parallel=FALSE)
    df <- data.frame("estimate"=ests$beta_hat, "se"=ests$se_hat, "statistic"=ests$beta_hat/ests$se_hat)
    return(df)
  }
 
  simnames <- c("nonpar", "par", "wfb", "wfb2", "naive", "selInf1", "ash")
	COVERAGE <- array(dim=c(length(simnames), p, n.rep))
	WIDTH <- array(dim=c(length(simnames),  p, n.rep))
	i <- 1
	while(i <= n.rep){
		cat(i, " ")
    #Calculate population effects
    z <- lr_pheno_effects_population(X.pop=X.pop, index=index,
                                     beta=beta, sd.err=1,
                                     which.sample=which.sample, parallel=parallel)
    y <- z$y
    effects <- z$effects

    #Run linear regressions for each collumn of X.sample
    data <- cbind(y, X.sample)
		res.orig <- lr_func(data)
		j <- order(abs(res.orig$statistic), decreasing=TRUE)
    rank <- match(1:p, j)
		cat("Got effect sizes.\n")

		#Non parametric bootstrap
		ci.nonpar <- nonpar_bs_ci(data, analysis.func = lr_func,n.rep=1000, res.orig=res.orig, level = 0.9, parallel=parallel)
		ci.nonpar <- ci.nonpar$ci
		COVERAGE[which(simnames=="nonpar"), ,i] <- (ci.nonpar[,1] < effects & effects < ci.nonpar[,2])[j]
		WIDTH[which(simnames=="nonpar"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[j]

		#Parametric bootstrap
		ci.par <- par_bs_ci(res.orig$estimate, res.orig$se,n.rep=1000, level=0.9)
		COVERAGE[which(simnames=="par"), ,i] <- (ci.par[,1] < effects & effects < ci.par[,2])[j]
		WIDTH[which(simnames=="par"), , i] <- (ci.par[,2] -ci.par[,1])[j]

		#WFB CIs conditional on taking the top 10%
		t <- res.orig$estimate/res.orig$se
		ct <- abs(t)[j][floor(0.1*p) + 1]
		wfb <- lapply(t, FUN=function(x){
			if(abs(x) < ct) return(c(NA, NA))
	 		ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
			return(ci)
			})
		ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=p)
		ci.wfb[,1] <- ci.wfb[,1]*res.orig$se
		ci.wfb[,2] <- ci.wfb[,2]*res.orig$se
		COVERAGE[which(simnames=="wfb"), ,i] <- (ci.wfb[,1] < effects & effects < ci.wfb[,2])[j]
		WIDTH[which(simnames=="wfb"), , i] <- (ci.wfb[,2] -ci.wfb[,1])[j]

		#WFB CIs sliding threshold
		wfb.2 <- lapply(t, FUN=function(x){
			u <- abs(t[abs(t) < abs(x)])
			if(length(u) > 0) ct <-max(u)
				else ct <- min(abs(t))/2
			ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
			if(class(ci) == "try-error") return(c(NA, NA))
			return(unlist(ci))
			})
		ci.wfb.2 <- matrix(unlist(wfb.2), byrow=TRUE, nrow=p)
		ci.wfb.2[,1] <- ci.wfb.2[,1]*res.orig$se
		ci.wfb.2[,2] <- ci.wfb.2[,2]*res.orig$se
		COVERAGE[which(simnames=="wfb2"),,i] <- (ci.wfb.2[,1] < effects & effects < ci.wfb.2[,2])[j]
		WIDTH[which(simnames=="wfb2"), , i] <- (ci.wfb.2[,2] -ci.wfb.2[,1])[j]

		#Naive ci
		ci.naive <- cbind(res.orig$estimate - res.orig$se*qt(0.95, df=p-1),
		                  res.orig$estimate + res.orig$se*qt(0.95, df=p-1))
		COVERAGE[which(simnames=="naive"), ,i] <- (ci.naive[,1] < effects & effects < ci.naive[,2])[j]
		WIDTH[which(simnames=="naive"), , i] <- (ci.naive[,2] -ci.naive[,1])[j]

		#Reid, Taylor, Tibshirani method (selectiveInference)
		M <- manyMeans(y=t, k=0.1*p, alpha=0.1, sigma=1)
		ci.rtt1 <- matrix(nrow=p, ncol=2)
		ci.rtt1[M$selected.set, ] <- M$ci
		ci.rtt1[,1] <- ci.rtt1[,1]*res.orig$se
		ci.rtt1[,2] <- ci.rtt1[,2]*res.orig$se
		COVERAGE[simnames == "selInf1", ,i]<- (ci.rtt1[,1] < effects & effects < ci.rtt1[,2])[j]
		WIDTH[simnames=="selInf1", , i] <- (ci.rtt1[, 2] - ci.rtt1[,1])[j]

		
		#ashr credible intervals
		ash.res <- ash(betahat = res.orig$estimate, sebetahat = res.orig$se, mixcompdist = "normal")
		ci.ash <- ashci(ash.res, level=0.9, betaindex = 1:p, trace=FALSE)
		COVERAGE[simnames == "ash", ,i]<- (ci.ash[,1] < effects & effects < ci.ash[,2])[j]
		WIDTH[simnames=="ash", , i] <- (ci.ash[, 2] - ci.ash[,1])[j]
		
		
		i <- i+1

	}

	cat("\n")
	return(list("COVERAGE"=COVERAGE,  "WIDTH"=WIDTH, "simnames"=simnames))
}

