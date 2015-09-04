easy_mvrnorm <- function (n = 1, mu, eS, tol=1e-6){
  p <- length(mu)
  if (!all(dim(eS$vectors) == c(p, p)))
    stop("incompatible arguments")
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  if (n == 1)
    drop(X)
  else t(X)
}

#Pairwise correlations: CI using fisher r to z transformation
#Implemented in the CIr function from psychometric
cor_ci <- function(W, n){
  p <- dim(W)[1]
  k <- p*(p-1)/2
  r <- W[upper.tri(W)]
  CI <- matrix(unlist(lapply(r, FUN=CIr, n=n)), nrow=k, byrow=TRUE)
  return(CI)
}

#Generate bootstrap confidence intervals from observed X
#Confidence intervals come out in the original order
#Get confidence intervals for z and convert to r
cor_bs_nonpar_ci <- function(X, w_hat, n.rep=2000){
  p <- dim(X)[2]
  n <- dim(X)[1]
  z_hat <- r2z(w_hat)
  B <- replicate(n=n.rep, expr={
    S <- sample(1:n, size=n, replace=TRUE)
    X.new <- X[S, ]
    W <- cor(X.new)
    z <- r2z(W[upper.tri(W)])
    k <- order(abs(z))
    sign(z[k])*(z[k]-z_hat[k])
  })
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(0.05, 0.95))})
  j <- order(abs(z_hat))
  my.ci.z <- cbind(z_hat[j]-qs[2,], z_hat[j]-qs[1,])
  which.neg <- which(z_hat[j] < 0)
  my.ci.z[ which.neg , ] <- cbind(z_hat[j][which.neg] + qs[1,which.neg], z_hat[j][which.neg]+qs[2,which.neg])
  #Re order
  k <- match(z_hat, z_hat[j])
  my.ci.z <- my.ci.z[k,]
  #Convert to r
  my.ci <- cbind(z2r(my.ci.z[,1]), z2r(my.ci.z[,2]))
  return(list("ci"=my.ci, "mean"=rowMeans(B)))
}


#Parametric bootstrap confidence intervals
#Use the approximate normal distribution or z = 1/2*log(1+x / 1-x) with sd 1/sqrt(n-3)
cor_bs_par_ci <-  function(w_hat, n, w_hat_mean=NULL,  n.rep=2000){
  if(is.null(w_hat_mean)) w_hat_mean <- w_hat
  z_hat <- r2z(w_hat)
  z_hat_mean <- r2z(w_hat_mean)
  nz <- length(z_hat)

  B <- replicate(n = n.rep, expr = {
    z <- rnorm(nz, mean=z_hat_mean, sd=1/sqrt(n-3))
    k <- order(abs(z))
    sign(z[k])*(z[k]-z_hat_mean[k])
  })

  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(0.05, 0.95))})

  j <- order(abs(z_hat))
  my.ci.z <- cbind(z_hat[j]-qs[2,], z_hat[j]-qs[1,])
  which.neg <- which(z_hat[j] < 0)
  my.ci.z[ which.neg , ] <- cbind(z_hat[j][which.neg] + qs[1,which.neg], z_hat[j][which.neg]+qs[2,which.neg])
  #Re order
  k <- match(z_hat, z_hat[j])
  my.ci.z <- my.ci.z[k,]
  #Convert to r
  my.ci <- cbind(z2r(my.ci.z[,1]), z2r(my.ci.z[,2]))
  return(my.ci)
}
