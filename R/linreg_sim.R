
#' Calculate expected coverage of several alternative intervals for
#'  marginal coefficients in linear regression.
#' @description This function can be used to generate simulations in section 3.1
#' @param X.pop a matrix of coefficients for an entire population.
#' @param which.sample Indices of individuals who are sampled
#' @param n.rep Number of simulations to run
#' @param which.X0 Which coefficients contribute to the phenotype
#' @param beta Effect size (including an intercept). Phenotypes whill be generated as
#' \code{y = beta \%*\% X0} where
#' \code{X0 = cbind(rep(1, length(which.sample)), X.pop[which.sample, which.X0])}
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
linreg_sim <- function(X.pop, which.sample, n.rep=10,
                      which.X0 = c(), beta=c(0)){

	n <-length(which.sample)
  p <- dim(X.pop)[2]
  X <- X.pop[which.sample, ]
  #X0 is matrix of variants with effects for sampled individuals
	if(length(which.X0)==0){
		X0 <- matrix(1, nrow=n, ncol=1)
	}else{
		X0 <- cbind(rep(1, n), X[,which.X0])
	}

  simnames <- c("nonpar", "par", "wfb", "wfb2", "naive")
	COVERAGE <- array(dim=c(length(simnames), p, n.rep))
	WIDTH <- array(dim=c(length(simnames),  p, n.rep))
	i <- 1
	while(i <= n.rep){
		cat(i, " ")
    #Calculate population effects
    z <- lr_pheno_effects_population(X.pop=X.pop, which.X0=which.X0,
                                     beta=beta, sd.err=1,
                                     which.sample=which.sample)
    y <- z$y
    effects <- z$effects

    #Run linear regressions for each collumn of X.sample
		my.bhat <-  many_lr(y, X)
		j <- order(abs(my.bhat$beta_hat/my.bhat$se_hat))

		#Non parametric bootstrap
		ci.nonpar <- lr_bs_nonpar_ci(y,X,my.bhat$beta_hat,my.bhat$se_hat,n.rep=1000)
		ci.nonpar <- ci.nonpar$ci
		COVERAGE[which(simnames=="nonpar"), ,i] <- (ci.nonpar[,1] < effects & effects < ci.nonpar[,2])[j]
		WIDTH[which(simnames=="nonpar"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[j]

		#Parametric bootstrap
		ci.par <- lr_bs_par_ci(my.bhat$beta_hat,my.bhat$se_hat,n.rep=1000)
		COVERAGE[which(simnames=="par"), ,i] <- (ci.par[,1] < effects & effects < ci.par[,2])[j]
		WIDTH[which(simnames=="par"), , i] <- (ci.par[,2] -ci.par[,1])[j]

		#WFB CIs conditional on taking the top 10%
		t <- my.bhat$beta_hat/my.bhat$se_hat
		ct <- quantile(abs(t), probs=0.9)
		wfb <- lapply(t, FUN=function(x){
			if(abs(x) < ct) return(c(NA, NA))
	 		ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
      if(class(ci) == "try-error") return(c(NA, NA))
			return(ci)
			})
		ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=p)
		ci.wfb[,1] <- ci.wfb[,1]*my.bhat$se_hat
		ci.wfb[,2] <- ci.wfb[,2]*my.bhat$se_hat
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
		ci.wfb.2[,1] <- ci.wfb.2[,1]*my.bhat$se_hat
		ci.wfb.2[,2] <- ci.wfb.2[,2]*my.bhat$se_hat
		COVERAGE[which(simnames=="wfb2"),,i] <- (ci.wfb.2[,1] < effects & effects < ci.wfb.2[,2])[j]
		WIDTH[which(simnames=="wfb2"), , i] <- (ci.wfb.2[,2] -ci.wfb.2[,1])[j]

		#Naive ci
		ci.naive <- cbind(my.bhat$beta_hat - my.bhat$se_hat*qt(0.95, df=p-1),
		                  my.bhat$beta_hat + my.bhat$se_hat*qt(0.95, df=p-2))
		COVERAGE[which(simnames=="naive"), ,i] <- (ci.naive[,1] < effects & effects < ci.naive[,2])[j]
		WIDTH[which(simnames=="naive"), , i] <- (ci.naive[,2] -ci.naive[,1])[j]

		i <- i+1

	}

	cat("\n")
	return(list("COVERAGE"=COVERAGE,  "WIDTH"=WIDTH, "simnames"=simnames))
}

