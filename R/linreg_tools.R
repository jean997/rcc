#'@import Matrix
#'@import MASS

#Generate data according to a linear model
#X0 is probably a small subset of the total number of predictors
lr_sim_pheno <- function(X0, beta, sd.err){
	n <- dim(X0)[1]
	err <- rnorm(n, 0, sd.err)
	beta <- matrix(beta, nrow=ncol(X0), ncol=1)
	y <- X0 %*% beta + err
}


#Run many linear regressions
many_lr <- function(y, X){
	p <- dim(X)[2]
	n <- dim(X)[1]
	Z <- lapply(1:p, FUN=function(i){
			if(all(X[,i] == X[1, i])) return(c(0, 0))
			my.X <- cbind(rep(1, n), X[,i])
			Xsum <- sum(X[,i])
			Xsqsum <- sum(X[,i]^2)
			XtXinv <- solve(matrix( c(n, Xsum, Xsum, Xsqsum), nrow=2, byrow=T) )
			Xty <- matrix(c(sum(y), sum(y*X[,i])), nrow=2)
			beta_hat <- XtXinv%*%Xty
			res <- y - my.X %*% beta_hat
			s_hat <-  1/(n-2)*sum(res^2)
			V <- XtXinv%*%t(my.X) %*% diag(rep(s_hat, n)) %*% my.X %*% XtXinv
			c(beta_hat[2], sqrt(V[2,2]))
		})
	Z <- matrix(unlist(Z), nrow=2, ncol=p, byrow=FALSE)
	return(list("beta_hat"=Z[1,], "se_hat"=Z[2,]))
}

#Generate bootstrap confidence intervals from observed y and X
#Confidence intervals come out in the original order
lr_bs_nonpar_ci <- function(y, X, beta_hat, se, n.rep=2000, alpha=0.1){
	p <- dim(X)[2]
	n <- dim(X)[1]
	B <- replicate(n=n.rep, expr={
			S <- sample(1:n, size=n, replace=TRUE)
			y.new <- y[S]
			X.new <- X[S, ]
			b <- many_lr(y.new, X.new)
			#order by t-statistic
			t.stat <- b$beta_hat/b$se_hat
			t.stat[ b$beta_hat ==0 & b$se_hat==0] <- 0
			k <- order(abs(t.stat))
			s <- sign(t.stat[k]); s[t.stat[k] ==0] <- 1
			s*(b$beta_hat[k] - beta_hat[k])
		})
	q1 <- alpha/2
	q2 <- 1-alpha/2
	qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})
  j <- order(abs(beta_hat/se))

  my.ci <- cbind(beta_hat[j]-qs[2,], beta_hat[j]-qs[1,])
	#Compensate for efect of abs value
	which.neg <- which(beta_hat[j] < 0)
	my.ci[ which.neg , ] <- cbind(beta_hat[j][which.neg] + qs[1,which.neg],
	                              beta_hat[j][which.neg]+qs[2,which.neg])

	jinv <- match(beta_hat, beta_hat[j])
	my.ci <- my.ci[jinv,]

	my.mean <- beta_hat[j] - sign(beta_hat[j])*rowMeans(B)
	my.mean <- my.mean[jinv]
  return(list("ci"=my.ci, "mean"=my.mean))
}

#Parametric bootstrap confidence intervals
lr_bs_par_ci <-  function(beta_hat, se, n.rep=2000, beta_hat_mean=NULL, alpha=0.1){
	p <- length(beta_hat)
	if(is.null(beta_hat_mean)) beta_hat_mean <- beta_hat

	B <- replicate(n = n.rep, expr = {
    w <- rnorm(p, mean=beta_hat_mean, sd=se)
		t.stat <- w/se
		k <- order(abs(t.stat))
    sign(t.stat[k])*(w[k]-beta_hat_mean[k])
		})

	q1 <- alpha/2
	q2 <- 1-alpha/2
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(q1, q2))})

  j <- order(abs(beta_hat/se))
	my.ci <- cbind(beta_hat[j]-qs[2,], beta_hat[j]-qs[1,])
	which.neg <- which(beta_hat[j] < 0)
	my.ci[ which.neg , ] <- cbind(beta_hat[j][which.neg] + qs[1,which.neg], beta_hat[j][which.neg]+qs[2,which.neg])
	jinv <- match(beta_hat, beta_hat[j])
	my.ci <- my.ci[jinv,]
  return(my.ci)
}

#Calculate population effects
lr_pheno_effects_population <- function(X.pop, which.X0, beta, sd.err, which.sample, use.exp=FALSE){
	stopifnot(length(beta) == length(which.X0) + 1)
	n <- dim(X.pop)[1]
	X0 <- cbind(rep(1, n), X.pop[, which.X0])
	if(use.exp){
		err <- rexp(n, rate=1) -1
	}else{
		err <- rnorm(n, 0, sd.err)
	}
	beta <- matrix(beta, nrow=ncol(X0), ncol=1)
	y <- X0 %*% beta + err
	effects.pop <- many_lr(y, X.pop)
	return(list("effects"=effects.pop$beta_hat, "y"=y[which.sample]))
}

lr_pheno_effects_sample <- function(X.sample, which.X0, beta, sd.err){
	stopifnot(length(beta) == length(which.X0) + 1)
	n <- dim(X.sample)[1]
	X0 <- cbind(rep(1, n), X.sample[, which.X0])
	err <- rnorm(n, 0, sd.err)
	beta <- matrix(beta, nrow=ncol(X0), ncol=1)
	y <- X0 %*% beta
	effects.sample <- many_lr(y, X.sample)
	y <- y+ err
	return(list("effects"=effects.sample$beta_hat, "y"=y))
}


#Parametric bootstrap confidence intervals
OLDlr_bs_par_ci <-  function(beta_hat, se, n.rep=2000, bhat_mean=NULL){
	p <- length(beta_hat)
	if(is.null(bhat_mean)) bhat_mean <- beta_hat

	B <- replicate(n = n.rep, expr = {
    w <- rnorm(p, mean=bhat_mean, sd=se)
		t.stat <- w/se
		k <- order(abs(t.stat))
    w[k]-bhat_mean[k]})

  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(0.05, 0.95))})

  j <- order(abs(beta_hat/se))
	my.ci <- cbind(beta_hat[j]-qs[2,], beta_hat[j]-qs[1,])
	k <- match(beta_hat, beta_hat[j])
	my.ci <- my.ci[k,]
  return(my.ci)
}

