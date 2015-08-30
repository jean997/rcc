
linreg_sims <- function(X, X.pop, n.rep=10,
                        which.X0 = c(), beta=c(0)){

	i <- 1
	n <- dim(X)[1]
	p <- dim(X)[2]

	simnames <- c("nonpar", "par", "wfb", "wfb2", "naive")
	WIDTH <- array(dim=c(9, p, n.rep))
	while(i <= n.rep){
		cat(i, " ")

		if(!is.null(X.pop)){
			z <- lr_pheno_effects_population(X.pop=X.pop, which.X0=which.X0, beta=beta, sd.err=1, which.sample=1:n)
			y <- z$y
			effects <- z$effects
		}else{
			y <- lr_sim_pheno(X0, beta=beta, sd.err=1)
		}

		my.bhat <-  many_lr(y, X)
		j <- order(abs(my.bhat$beta_hat/my.bhat$se_hat))


		#Non parametric
		ci.nonpar <- lr_bs_nonpar_ci(y, X, my.bhat$beta_hat, my.bhat$se_hat, n.rep=1000, skip.na=TRUE)
		mean.nonpar <- ci.nonpar$mean
		ci.nonpar <- ci.nonpar$ci
		C_NPAR[,i] <- (ci.nonpar[,1] < effects & effects < ci.nonpar[,2])[j]
		WIDTH[which(wnames=="nonpar"), , i] <- (ci.nonpar[,2] -ci.nonpar[,1])[j]


		ci.par <- lr_bs_par_ci(my.bhat$beta_hat, my.bhat$se_hat, n.rep=1000)
		C_PAR[,i] <- (ci.par[,1] < effects & effects < ci.par[,2])[j]
		WIDTH[which(wnames=="par"), , i] <- (ci.par[,2] -ci.par[,1])[j]

		ci.par2 <- lr_bs_par_ci(my.bhat$beta_hat, my.bhat$se_hat, n.rep=1000, bhat_mean=mean.nonpar)
		C_PAR2[,i] <- (ci.par2[,1] < effects & effects < ci.par2[,2])[j]
		WIDTH[which(wnames=="par2"), , i] <- (ci.par2[,2] -ci.par2[,1])[j]

		ci.opar <- lr_bs_par_ci(my.bhat$beta_hat, my.bhat$se_hat, n.rep=1000, bhat_mean=effects)
		C_OPAR[,i] <- (ci.opar[,1] < effects & effects < ci.opar[,2])[j]
		WIDTH[which(wnames=="opar"), , i] <- (ci.opar[,2] -ci.opar[,1])[j]

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
		C_WFB[,i] <- (ci.wfb[,1] < effects & effects < ci.wfb[,2])[j]
		WIDTH[which(wnames=="wfb"), , i] <- (ci.wfb[,2] -ci.wfb[,1])[j]

		t <- (my.bhat$beta_hat/my.bhat$se_hat)
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
		C_WFB2[,i] <- (ci.wfb.2[,1] < effects & effects < ci.wfb.2[,2])[j]
		WIDTH[which(wnames=="wfb2"), , i] <- (ci.wfb.2[,2] -ci.wfb.2[,1])[j]

		bhat.stefan <- my.bhat$beta_hat + stefan.bias(my.bhat$beta_hat)
		ci.stefan <- lr_bs_par_ci(my.bhat$beta_hat, my.bhat$se_hat, n.rep=1000, bhat_mean=bhat.stefan)
		C_ST[,i] <- (ci.stefan[,1] < effects & effects < ci.stefan[,2])[j]
		WIDTH[which(wnames=="stefan"), , i] <- (ci.stefan[,2] -ci.stefan[,1])[j]

		ci.naive <- cbind(my.bhat$beta_hat - my.bhat$se_hat*qt(0.95, df=p-1), my.bhat$beta_hat + my.bhat$se_hat*qt(0.95, df=p-2))
		C_N[,i] <- (ci.naive[,1] < effects & effects < ci.naive[,2])[j]
		WIDTH[which(wnames=="naive"), , i] <- (ci.naive[,2] -ci.naive[,1])[j]

		ci.bf <- cbind(my.bhat$beta_hat - my.bhat$se_hat*qt(1-(0.05/p), df=p-1), my.bhat$beta_hat + my.bhat$se_hat*qt(1-(0.05/p), df=p-2))
		C_BF[,i] <- (ci.bf[,1] < effects & effects < ci.bf[,2])[j]
		WIDTH[which(wnames=="bf"), , i] <- (ci.bf[,2] -ci.bf[,1])[j]
		i <- i+1

	}

	cat("\n")
	return(list("par2"=C_PAR2, "stefan"=C_ST, "naive"=C_N, "nonpar"=C_NPAR, "par"=C_PAR, "bf"=C_BF, "wfb_10"=C_WFB, "n"=n, "p"=p, "Sigma"=Sigma, "which.X0"=which.X0, "beta"=beta, "wfb_2" =C_WFB2, "oracle_par"=C_OPAR, "width"=WIDTH, "wnames"=wnames))


}


plot_sims <- function(sims, top.x=0.2, main="", which.set=c(1, 2)){
  p <- dim(sims$par)[1]
	k <- max(1, floor(p-(p*top.x)))
	naive <- rowMeans(sims$naive)[k:p]
	oracle_par <- rowMeans(sims$oracle_par)[k:p]
	par2 <- rowMeans(sims$par2)[k:p]
	stefan <- rowMeans(sims$stefan)[k:p]
	nonpar <- rowMeans(sims$nonpar)[k:p]
	par <- rowMeans(sims$par)[k:p]
	wfb <- rowMeans(sims$wfb_10)[k:p]
	wfb2 <- rowMeans(sims$wfb_2, na.rm=TRUE)[k:p]
	#ylim <- range(c(naive, nonpar, par, wfb, wfb2), na.rm=TRUE)
	ylim <- range(c(naive, par2, stefan, nonpar, par, wfb, wfb2), na.rm=TRUE)
	plot(k:p, k:p, type="n", xlab="Order Statistic",ylim=ylim, ylab="Coverage", main=main, col="green")
	if(2 %in% which.set){
		points(k:p, wfb2, col="darkgreen", pch=3)
		points(k:p, wfb, col="green", pch=3)
	}
	if(1 %in% which.set){
		points(k:p, naive, col="black")
		points(k:p, par, col="red", pch=4)
		#points(k:p, par2, col="orange")
		#points(1:p, rowMeans(sims$stefan), col="purple")
		points(k:p, nonpar, col="blue", pch=4)
	}
	#points(k:p, oracle_par, col="black")
	abline(h=0.9)
	legend("bottom", legend=c("Std", "Non-Par", "Par", "WFB", "WFB2"), col=c("black", "blue", "red", "green", "darkgreen"), pch=c(1, 4, 4, 3, 3))

}

plot_width <- function(sims, top.x=0.2, main="", which.set=c(1, 2)){
  p <- dim(sims$par)[1]
  k <- max(1, floor(p-(p*top.x)))
	N <- sims$wnames
  naive <- rowMeans(sims$width[which(N=="naive"), , ])[k:p]
  oracle_par <- rowMeans(sims$width[which(N=="opar"), , ])[k:p]
  par2 <- rowMeans(sims$width[which(N=="par2"), , ])[k:p]
  par <- rowMeans(sims$width[which(N=="par"), , ])[k:p]
  nonpar <- rowMeans(sims$width[which(N=="nonpar"), , ])[k:p]
  stefan <- rowMeans(sims$width[which(N=="stefan"), , ])[k:p]
  wfb <- rowMeans(sims$width[which(N=="wfb"), , ], na.rm=TRUE)[k:p]
  wfb2 <- rowMeans(sims$width[which(N=="wfb2"), , ], na.rm=TRUE)[k:p]
  #ylim <- range(c(naive, nonpar, par, wfb, wfb2), na.rm=TRUE)
  ylim <- range(c(naive, par2, nonpar, par, wfb, wfb2), na.rm=TRUE)
	plot(k:p, k:p, type="n", xlab="Order Statistic",ylim=ylim, ylab="Avg. Width", main=main, col="green")
	if(2 %in% which.set){
		points(k:p, wfb2, col="darkgreen", pch=3)
		points(k:p, wfb, col="green", pch=2)
	}
	if(1 %in% which.set){
		points(k:p, naive, col="black")
		points(k:p, par, col="red", pch=4)
		#points(k:p, par2, col="orange")
		#points(1:p, rowMeans(sims$stefan), col="purple")
		points(k:p, nonpar, col="blue", pch=4)
	}
	legend("bottom", legend=c("Std", "Non-Par", "Par", "WFB", "WFB2"), col=c("black", "blue", "red", "green", "darkgreen"), pch=c(1, 4, 4, 3, 3))
}
