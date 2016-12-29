biomarker_nonpar_ci <- function(dat, stats, cutpoints, n.rep=500, level=0.1){
  n <- dim(dat)[1]
  B <- replicate(n=n.rep, expr={
    S <- sample(1:n, size=n, replace=TRUE)
    dat.new <- dat[S, ]
    stats.new <- lapply(cutpoints, FUN=function(xx){
      if(sum(dat$x >=xx) <= 2) return(rep(NA, 3))
      f <- lm(y~trt, data=dat.new[dat.new$x >= xx,])
      coef <- summary(f)$coefficients
      if(nrow(coef)==1) return(rep(NA, 3))
      return(coef[2, 1:3])
    })
    stats.new <- data.frame(matrix(unlist(stats.new), byrow=TRUE, ncol=3))
    names(stats.new) <- c("beta", "se", "tstat")
    k <- order(abs(stats.new$tstat))
    sign(stats.new$beta[k])*(stats.new$beta[k]-stats$truth[k])
  })
  qs <- apply(B, MARGIN=1, FUN=function(x){quantile(x, probs=c(level/2, 1-(level/2)), na.rm=TRUE)})
  j <- order(abs(stats$tstat))
  my.ci <- cbind(stats$beta[j]-qs[2,], stats$beta[j]-qs[1,])
  which.neg <- which(stats$beta[j] < 0)
  my.ci[ which.neg , ] <- cbind(stats$beta[j][which.neg] + qs[1,which.neg],
                              stats$beta[j][which.neg]+qs[2,which.neg])
  #Re order
  k <- match(stats$beta, stats$beta[j])
  my.ci <- my.ci[k,]
  return(my.ci)
}