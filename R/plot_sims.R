#'Plot coverage of example simulations
#'@param R A list produced by \code{example_sim}
#'@param simnames Names of CIs to plot. The order they are given will be the plotting order.
#'@param legend.names Names to print in legend
#'@param cols,shapes colors and shapes to use.
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@return A ggplot object
#'@export
plot_coverage <- function(R, simnames, legend.names,
                                  cols, shapes, main="", proportion=0.2,
                                    y.axis.off=FALSE, set.range=FALSE,
                                    legend.position=c(0.28, 0.4)){
  which.keep <- which(R$simnames %in% simnames)
  p <- dim(R$COVERAGE)[2]
  k <- floor((1-proportion)*p)
  k <- max(k, 1)
  avg.coverage <- matrix(nrow=(p-k+1), ncol=length(which.keep))
  nms <- c()
  for(i in 1:length(which.keep)){
    w <- which.keep[i]
    nms <- c(nms, R$simnames[w])
    avg.coverage[,i] <- rowMeans(R$COVERAGE[w, , ], na.rm=TRUE)[k:p]
    #na.rm=TRUE for WFB method which sometimes has errors and doesn't generate intervals
  }
  avg.coverage <- data.frame(avg.coverage)
  names(avg.coverage) <- nms
  avg.coverage$Rank <- k:p
  avg.coverage.long <- gather(avg.coverage, "Method", "RCC", -Rank)
  if(set.range) ylim=range(avg.coverage.long$RCC)
  #Re-order factor levels
  avg.coverage.long$Method <- factor( as.character(avg.coverage.long$Method),
                                      levels=simnames)
  m <- avg.coverage.long$Method
  o <- c()
  for(i in 1:length(simnames)){
    o <- c(o, which(m==simnames[i]))
  }
  avg.coverage.long <- avg.coverage.long[o, ]
  h <- ggplot(avg.coverage.long, aes(x=Rank)) + geom_hline(aes(yintercept=0.9)) +
    geom_point(aes(y=RCC,  color=Method, shape=Method)) +
    theme_bw(base_size = 14) +
    scale_shape_manual(values=shapes, labels=legend.names) +
    scale_color_manual(values=cols, labels=legend.names) +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

  h <- h + theme(legend.position=legend.position,
                      legend.key = element_rect(color="white", fill="white", size = 0.6, linetype=0),
                      legend.background = element_rect(fill="transparent"),
                      legend.title= element_text(size=10),
                      legend.text=element_text(size=9))
  if(y.axis.off) h <- h + theme(axis.title.y = element_blank())
  if(!main == "") h <- h + ggtitle(main) + theme(plot.title=element_text(size=20))
  if(set.range) return(h + ylim(ylim))
  return(h+ ylim(c(0, 1)) )

}

#'Plot average width of example simulations
#'@param R A list produced by \code{example_sim}
#'@param simnames Names of CIs to plot. The order they are given will be the plotting order.
#'@param legend.names Names to print in legend
#'@param cols,shapes colors and shapes to use.
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@return A ggplot object
#'@export
plot_width <- function(R, simnames, cols, shapes, legend.names,
                               main="", proportion=0.2,
                               y.axis.off=FALSE,
                               legend.position="none"){
  which.keep <- which(R$simnames %in% simnames)
  p <- dim(R$WIDTH)[2]
  k <- floor((1-proportion)*p)
  k <- max(k, 1)

  ylim=range(apply(R$WIDTH[which.keep, , ], MARGIN=1,  FUN=function(x){ z <- rowMeans(x); return(z[k:p])}), na.rm=TRUE)

  avg.width <- matrix(nrow=(p-k+1), ncol=length(which.keep))
  nms <- c()
  for(i in 1:length(which.keep)){
    w <- which.keep[i]
    nms <- c(nms, R$simnames[w])
    avg.width[,i] <- rowMeans(R$WIDTH[w, , ])[k:p]
  }
  avg.width <- data.frame(avg.width)
  names(avg.width) <- nms
  avg.width$Rank <- k:p
  avg.width.long <- gather(avg.width, "Method", "Average Width", -Rank)
  avg.width.long$Method <- factor( as.character(avg.width.long$Method),
                                   levels=simnames)
  m <- avg.width.long$Method
  o <- c()
  for(i in 1:length(simnames)){
    o <- c(o, which(m==simnames[i]))
  }
  avg.width.long <- avg.width.long[o, ]

  h <- ggplot(avg.width.long, aes(x=Rank)) +
    geom_point(aes(y=`Average Width`, group=Method,  color=Method, shape=Method)) +
    theme_bw(base_size = 14) + ylim(ylim) +
    scale_shape_manual(values=shapes, labels=legend.names) +
    scale_color_manual(values=cols, labels=legend.names) +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

  h <- h + theme(legend.position=legend.position,
                 legend.key = element_rect(color="white", fill="white", size = 0.6, linetype=0),
                 legend.background = element_rect(fill="transparent"),
                 legend.title= element_text(size=10),
                 legend.text=element_text(size=9))

  if(y.axis.off) h <- h + theme(axis.title.y = element_blank())
  if(!main == "") h <- h + ggtitle(main) + theme(plot.title=element_text(size=20))
  return(h)

}

