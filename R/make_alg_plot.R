
make_alg_plots <- function(seed, idx = 1, n=3, lwd =2, text.size=7,
                           axis.text.size=15, axis.labs.size=20,
                           min.dist.mc = 0.3, min.dist.delta=0.3){
  set.seed(seed)
  X <- replicate(n=500, expr = rnorm(n=1000, mean=example_params[,2], sd=1))
  ix <- apply(X, MARGIN=2, FUN=function(x){
    order(x, decreasing = TRUE)[idx]
  })

  delta <- X[cbind(ix, 1:500)]-example_params[ix,2]

  qs <- quantile(delta, probs=c(0.05, 0.95))

  delta <- data.frame(delta)

  ret <- list()

  #Montecarlo plots
  for(i in 1:n){
    dat <- data.frame(cbind(X[,i], example_params[,2]))
    names(dat)=c("sample", "theta")
    t1hat <- sort(X[,i], decreasing=TRUE)[idx]
    t1 <- example_params[order(X[,i], decreasing=TRUE)[idx],2]
    ss <- paste0(i, ",", idx)
    ss2 <- paste0("(", idx, ")")
    my.at <- c(t1, t1hat)
    if(abs(t1hat-t1) < min.dist.mc){
      b <- min.dist.mc- abs(t1hat-t1)
      if(t1hat < t1) my.at <- c(t1 + b, t1hat -b)
        else my.at <- c(t1-b, t1hat + b)
    }
    h <- ggplot(dat) + theme_bw() +
      geom_density(aes(theta), lty=0, fill="orange", alpha=0.3) +
      geom_density(aes(sample), lty=0, fill="purple", alpha=0.3) +
      geom_segment(x=t1hat, xend=t1hat, y=0, yend=1) + geom_segment(x=t1, xend=t1, y=0, yend=1, lty=2) +
      geom_segment(x=t1hat, xend=t1, y=0.05, yend=0.05, lty=3) +
      scale_x_continuous(breaks=my.at,
                         limits=range(X[, 1:n]),
                         labels=c( as.expression(bquote(theta["s"[.(i)]~.(ss2)])),
                                  as.expression(bquote(vartheta[.(ss)])))) +
      theme(axis.title=element_blank(), axis.text.y = element_blank(),
            panel.grid.minor=element_blank(), panel.grid.major = element_blank(),
            axis.ticks = element_blank(), rect=element_blank(),
            axis.text.x=element_text(size=axis.text.size)) #+
      #annotate(geom="text", x=t1 + (t1hat-t1)/2, y=0.12,
      #         label=paste0("tilde(delta)['", i, ",[", idx, "]']") , parse=TRUE, size=text.size)
    if(i==1) h <- h + ggtitle("Step 1: Sample from G")
    ret[[i]] <- h

  }



  ss <- paste0("[", idx, "]")
  my.at <- c(0, qs[1], qs[2])
  if(abs(qs[1]) < min.dist.delta){
    b <- min.dist.delta - abs(qs[1])
    my.at[2] <- qs[1]+ sign(qs[1])*b
  }
  if(abs(qs[2]) < min.dist.delta){
    b <- min.dist.delta - abs(qs[2])
    my.at[3] <- qs[2]+ sign(qs[2])*b
  }
  ret[[i+1]] <- ggplot(delta) + geom_density(aes(x=delta, y=..scaled..), fill="blue", alpha=0.3, lty=0) +
    theme_bw() +
    geom_vline(xintercept=0) + ggtitle(bquote("Step 2a: Quantiles of "~tilde(delta)[.(ss)])) +
    geom_vline(xintercept=qs[1], lty=2) + geom_vline(xintercept=qs[2], lty=2) +
    geom_segment(x=0, xend=qs[1], y=0.1, yend=0.1, col="violetRed", lty=1, lwd=lwd) +
    geom_segment(x=0, xend=qs[2], y=0.13, yend=0.13, col="violetRed", lty=1, lwd=lwd) +
    geom_segment(x=0, xend=qs[2], y=0.13, yend=0.13, col="violetred1", lty=2, lwd=lwd) +
    scale_x_continuous(breaks=my.at,
                     limits=c(min(0, min(delta)), max(delta)),
                     labels=c("0", as.expression(bquote(tilde(H)[.(ss)]^"-1"~"(0.05)")),
                              as.expression(bquote(tilde(H)[.(ss)]^"-1"~"(0.95)")))) +
  theme(axis.title=element_blank(), axis.text.y = element_blank(),
        panel.grid.minor=element_blank(), panel.grid.major = element_blank(),
        axis.ticks = element_blank(), rect=element_blank(),
        axis.title.x = element_text(size=axis.labs.size),
        axis.text.x=element_text(size=axis.text.size))
        #axis.text.x=element_text(hjust=c(0.9, 0, 0), size=axis.text.size))


  x = rnorm(n=1000, mean=example_params[,2], sd=1)
  t1hat <- sort(x, decreasing=TRUE)[idx]
  t1 <- example_params[order(x, decreasing=TRUE)[idx],2]
  ci <- c(t1hat - qs[2], t1hat-qs[1])
  cinorm <- t1hat + c(-1, 1)*qnorm(0.95)
  x <- data.frame(x)
  R = range(c(0, ci, t1, t1hat, cinorm))

  si <- paste0("s(", idx, ")")
  ssi <- paste0("(", idx, ")")
  ret[[i+2]] <- ggplot(x) + theme_bw() +
    geom_hline(yintercept = 0) + ylim(c(0, 1)) + ggtitle("Step 2b: Pivot")+
    geom_segment(x=t1hat, xend=t1hat, y=0, yend = 1) +
    geom_segment(x=t1, xend=t1, y=0, yend = 1, lty=2) +
    #Oracle Montecarlo interval
    geom_segment(x=t1hat, xend=t1hat-qs[1], y=0.1, yend=0.1, col="violetRed", lty=1, lwd=lwd) +
    geom_segment(x=t1hat, xend=t1hat-qs[2], y=0.13, yend=0.13, col="violetRed", lty=1, lwd=lwd) +
    geom_segment(x=t1hat, xend=t1hat-qs[2], y=0.13, yend=0.13, col="violetred1", lty=2, lwd=lwd) +
    annotate(geom="text", x=t1hat-qs[1], y=0.13, col="violetred", label=")", size=axis.labs.size) +
    annotate(geom="text", x=t1hat-qs[2], y=0.13, col="violetred", label="(", size=axis.labs.size) +
    annotate(geom="text", x=(ci[1] + 0.25*abs(ci[2]-ci[1])), y=0.2,
             label=paste0("CI['s(", idx, ")']^'oracle'"), parse=TRUE, size=text.size) +
    geom_segment(x=cinorm[1], xend=cinorm[2], y=0.7, yend=0.7, color="blue", lty=1, lwd=lwd) +
    annotate(geom="text", x=cinorm[1], y=0.7, col="blue", label="(", size=axis.labs.size) +
    annotate(geom="text", x=cinorm[2], y=0.7, col="blue", label=")", size=axis.labs.size) +
    annotate(geom="text", x=(cinorm[1] + 0.75*abs(cinorm[2]-cinorm[1])), y=0.77,
             label=paste0("CI['s(", idx, ")']^'naive'"), parse=TRUE, size=text.size) +
    scale_x_continuous(breaks=c(0,  t1,  t1hat), limits=R,
                     labels=c("0", as.expression(bquote(theta[.(si)])),
                      as.expression(bquote(hat(theta)[.(ssi)])))) +
    theme(axis.title=element_blank(), axis.text.y = element_blank(), rect = element_blank(),
        panel.grid.minor=element_blank(), panel.grid.major = element_blank(),
        axis.ticks.x= element_blank(), axis.text.x=element_text(size=axis.text.size))

  return(ret)
}

