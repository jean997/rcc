#'Plot coverage of example simulations
#'@param R A list produced by \code{example_sim}
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@return A ggplot object
#'@export
plot_example_coverage <- function(R, main="", proportion=0.2,
                                    y.axis.off=FALSE,
                                    legend.position=c(0.28, 0.4)){
  which.keep <- which(R$simnames %in% c("naive", "wfb", "wfb2", "par", "oracle"))
  p <- dim(R$COVERAGE)[2]
  k <- floor((1-proportion)*p)
  k <- max(k, 1)
  avg.coverage <- matrix(nrow=(p-k+1), ncol=length(which.keep))
  nms <- c()
  for(i in 1:length(which.keep)){
    w <- which.keep[i]
    nms <- c(nms, R$simnames[w])
    avg.coverage[,i] <- rowMeans(R$COVERAGE[w, , ])[k:p]
  }
  avg.coverage <- data.frame(avg.coverage)
  names(avg.coverage) <- nms
  avg.coverage$Rank <- k:p
  avg.coverage.long <- gather(avg.coverage, "Method", "RCC", -Rank)

  avg.coverage.long$Method <- factor( as.character(avg.coverage.long$Method),
                                      levels=c("naive", "wfb", "wfb2", "oracle", "par"))
  m <- avg.coverage.long$Method
  o <- c( which(m=="naive"), which(m=="wfb2"), which(m=="par"), which(m=="wfb"), which(m=="oracle") )
  avg.coverage.long <- avg.coverage.long[o, ]
  h <- ggplot(avg.coverage.long, aes(x=Rank)) + geom_hline(aes(yintercept=0.9)) +
    geom_point(aes(y=RCC,  color=Method, shape=Method)) +
    theme_bw(base_size = 14) + ylim(c(0, 1)) +
    scale_shape_manual(values=c(1:4, 0), labels=c("Marginal", "WFB", "WFB-Sliding",
                                                  "Oracle", "Parametric Bootstrap")) +
    scale_color_manual(values=c("black", "gold4", "forestgreen", "blue", "deeppink3"), labels=c("Marginal", "WFB", "WFB-Sliding",
                                                                                                "Oracle", "Parametric Bootstrap")) +
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

#'Plot average width of example simulations
#'@param R A list produced by \code{example_sim}
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@return A ggplot object
#'@export
plot_example_width <- function(R, main="", proportion=0.2,
                                y.axis.off = TRUE,
                               legend.position="none"){
  which.keep <- which(R$simnames %in% c("naive", "wfb", "wfb2", "par", "oracle"))
  p <- dim(R$WIDTH)[2]
  k <- floor((1-proportion)*p)
  k <- max(k, 1)

  ylim=range(apply(R$WIDTH[which.keep, , ],
                   MARGIN=1,  FUN=function(x){ z <- rowMeans(x); return(z[k:p])}), na.rm=TRUE)

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
                                   levels=c("naive", "wfb", "wfb2", "oracle", "par"))
  m <- avg.width.long$Method
  o <- c( which(m=="naive"), which(m=="wfb2"), which(m=="par"), which(m=="wfb"), which(m=="oracle") )
  avg.width.long <- avg.width.long[o, ]

  h <- ggplot(avg.width.long, aes(x=Rank)) +
    geom_point(aes(y=`Average Width`, group=Method,  color=Method, shape=Method)) +
    theme_bw(base_size = 14) + ylim(ylim) +
    scale_shape_manual(values=c(1:4, 0), labels=c("Marginal", "WFB", "WFB-Sliding",
                                                  "Oracle", "Parametric Bootstrap")) +
    scale_color_manual(values=c("black", "gold4", "forestgreen", "blue", "deeppink3"),
                       labels=c("Marginal", "WFB", "WFB-Sliding",
                                "Oracle", "Parametric Bootstrap")) +
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

