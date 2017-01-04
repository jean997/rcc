#'@import tidyr
#'@import ggplot2
#'@import psychometric
#'@import selectiveInference
#'@import gridExtra
#'@import ashr

#Run many linear regressions
many_lr <- function(y, X, parallel=FALSE){
	p <- dim(X)[2]
	n <- dim(X)[1]

  my_lr_func <- function(i){
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
    return(c(beta_hat[2], sqrt(V[2,2])))
  }	
	if(!parallel){
	  Z <- lapply(1:p, FUN=function(i){
			my_lr_func(i)
		})
	}else{
	  cores <- detectCores()-1
	  cl <- makeCluster(cores, type="FORK")
	  on.exit(stopCluster(cl))
	  Z <- parLapply(cl, 1:p, fun=function(i){
	    my_lr_func(i)
	  })
	}
	Z <- matrix(unlist(Z), nrow=2, ncol=p, byrow=FALSE)
	df <- data.frame("beta_hat"=Z[1,], "se_hat"=Z[2,])
	return(df)
}

#Generate bootstrap confidence intervals from observed y and X
#Confidence intervals come out in the original order
lr_bs_nonpar_ci <- function(y, X, beta_hat, se, n.rep=2000, alpha=0.1, parallel=FALSE){
	p <- dim(X)[2]
	n <- dim(X)[1]
	if(!parallel){
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
	}else{
	  cores <- detectCores()-1
	  cl <- makeCluster(cores, type="FORK")
	  on.exit(stopCluster(cl))
	  B <- parSapply(cl, 1:n.rep, FUN=function(ii){
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
	}
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

#Calculate population effects
lr_pheno_effects_population <- function(X.pop, index, beta, sd.err, which.sample, parallel=TRUE){
	stopifnot(length(beta) == length(index))
	n <- dim(X.pop)[1]
	err <- rnorm(n, 0, sd.err)
	if(length(index)==0){
	  y <- err	
	}else{
	  beta <- matrix( beta, nrow=length(index), ncol=1)
	  y <- X.pop[, index] %*% beta + err
	}
	effects.pop <- many_lr(y, X.pop, parallel=parallel)
	return(list("effects"=effects.pop$beta_hat, "y"=y[which.sample]))
}
