#' Calculate expected coverage of several alternative intervals for
#'  marginal coefficients in linear regression.
#' @description This function can be used to generate simulations in section 3.1
#' @param Sigma Covariance matrix from which data will be generated.
#' @param n Number of samples
#' @param n.rep Number of replications
#' @param eS Optional eigen decomposition of Sigma
#' @return A list with elements:
#' \describe{
#' \item{\code{COVERAGE}}{A 5 by p by n array giving the coverage at each simulation.
#' \code{COVERAGE[i, j, k]} is either 0 or 1 indicating weather the interval formed using
#' method i for rank j in simulation k covered its respective parameter}
#' \item{\code{WIDTH}}{ A 5 by p by n array giving the width of each interval in each simulation.}
#'  \item{\code{simnames}}{Name of each type of the 5 types of confidence intervals}
#' }
#'@export
cor_sim <- function(Sigma, n, n.rep=10,  eS=NULL){

  p <- dim(Sigma)[1]
	if(is.null(eS)) eS <- eigen(Sigma)

	i <- 1
	D <- diag(1/sqrt(diag(Sigma)))
	S <- D%*% Sigma %*% D
	r <- S[upper.tri(S)]
	nr <- length(r)

  simnames <- c("nonpar", "par", "wfb", "wfb2", "naive", "selInf1")
	COVERAGE <- array(dim=c(length(simnames), nr, n.rep))
	WIDTH <- array(dim=c(length(simnames),  nr, n.rep))

	while(i <= n.rep){
		X <- easy_mvrnorm(n=n, mu=rep(0, p), eS=eS)
		cat(i, " ")
		W_hat <- cor(X)
		w_hat <- W_hat[upper.tri(W_hat)]
		j <- order(abs(w_hat))

		#Non parametric bootstrap
		ci.nonpar <- cor_bs_nonpar_ci(X, w_hat, n.rep=1000)
		ci.nonpar <- ci.nonpar$ci
		COVERAGE[which(simnames=="nonpar"), ,i] <- (ci.nonpar[,1] < r & r < ci.nonpar[,2])[j]
		WIDTH[which(simnames=="nonpar"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[j]

		#Parametric bootstrap
		ci.par <- cor_bs_par_ci(w_hat, n, n.rep=1000)
		COVERAGE[which(simnames=="par"), ,i] <- (ci.par[,1] < r & r < ci.par[,2])[j]
		WIDTH[which(simnames=="par"), , i] <- (ci.par[,2]-ci.par[,1])[j]

		#Parametric bootstrap with oracle mean (no covariance)
		ci.opar <- cor_bs_par_ci(w_hat, n, n.rep=1000, w_hat_mean=r)
		COVERAGE[which(simnames=="opar"), ,i] <- (ci.opar[,1] < r & r < ci.opar[,2])[j]
		WIDTH[which(simnames=="opar"), , i] <- (ci.opar[,2]-ci.opar[,1])[j]

		#WFB CIs conditional on taking the top 10%
		z_hat <- r2z(w_hat)
		z_hat_scaled <- sqrt(n-3)*z_hat
		ct <- quantile(abs(z_hat_scaled), probs=0.9)
		wfb <- lapply(z_hat_scaled, FUN=function(x){
			if(abs(x) < ct) return(c(NA, NA))
			ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
			if(class(ci) == "try-error") return(c(NA, NA))
			return(unlist(ci))
			})
		ci.wfb <- matrix(unlist(wfb), byrow=TRUE, nrow=nr)
		ci.wfb<- ci.wfb/sqrt(n-3)
		ci.wfb[,1] <- z2r(ci.wfb[,1])
		ci.wfb[,2] <- z2r(ci.wfb[,2])
		COVERAGE[which(simnames=="wfb"), ,i] <- (ci.wfb[,1] < r & r < ci.wfb[,2])[j]
		WIDTH[which(simnames=="wfb"), , i] <- (ci.wfb[,2]-ci.wfb[,1])[j]
		if(any(!is.na(WIDTH[which(simnames=="wfb"), 1:4455, i]))){
			save(X, j, ci.wfb, file="wfberror.RData")
			return(0)
		}

		#WFB CIs moving threshold
		wfb.2 <- lapply(z_hat_scaled, FUN=function(x, z){
			u <- abs(z[abs(z) < abs(x)])
			if(length(u) > 0) ct <-max(u)
				else ct <- min(abs(z))/2
			ci <- try(Shortest.CI(x, ct=ct, alpha=0.1), silent=TRUE)
			if(class(ci) == "try-error") return(c(NA, NA))
			return(unlist(ci))
			}, z=z_hat_scaled)
		ci.wfb.2 <- matrix(unlist(wfb.2), byrow=TRUE, nrow=nr)
		ci.wfb.2<- ci.wfb.2/sqrt(n-3)
		ci.wfb.2[,1] <- z2r(ci.wfb.2[,1])
		ci.wfb.2[,2] <- z2r(ci.wfb.2[,2])
		COVERAGE[which(simnames=="wfb2"), ,i] <- (ci.wfb.2[,1] < r & r < ci.wfb.2[,2])[j]
		WIDTH[which(simnames=="wfb2"), , i] <- (ci.wfb.2[,2]-ci.wfb.2[,1])[j]

		#Naive ci
		ci.naive <- cor_ci(W_hat, n)
		COVERAGE[which(simnames=="naive"), ,i] <- (ci.naive[,1] < r & r < ci.naive[,2])[j]
		WIDTH[which(simnames=="naive"), , i] <- (ci.naive[,2]-ci.naive[,1])[j]

		#Reid, Taylor, Tibshirani method (selectiveInference)
		M <- manyMeans(y=z_hat_scaled, k=0.1*nr, alpha=0.1, sigma=1)
		ci.rtt1 <- matrix(nrow=nr, ncol=2)
		ci.rtt1[M$selected.set, ] <- M$ci
		COVERAGE[simnames == "selInf1", ,i]<- (ci.rtt1[,1] < r & rs < ci.rtt1[,2])[j]
		WIDTH[simnames=="selInf1", , i] <- (ci.rtt1[, 2] - ci.rtt1[,1])[j]


		i <- i+1

	}

	cat("\n")
	return(list("COVERAGE"=COVERAGE, "WIDTH"=WIDTH, "simnames"=simnames))


}

