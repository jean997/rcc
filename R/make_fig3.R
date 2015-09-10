

make_fig3 <- function(){
  titles <- c("Low Correlation", "High Correlation")
  plot.list.c <- list()
  Z <- list()
  Z[[1]] <- rcc:::getobj("chr1_sim_results.RData")
  Z[[2]] <- rcc:::getobj("dnarep_sim_results.RData")

  p <- 4950
  k <- floor(0.8*p)
  mincov <- min(rowMeans(Z[[1]]$COVERAGE[1, , ])[k:p], na.rm=TRUE)
  for(i in 2:5){
    mincov <- min(mincov, min(rowMeans(Z[[1]]$COVERAGE[i, , ])[k:p], na.rm=TRUE))
  }
  for(i in 1:2){

    main =  titles[i]

    if(i ==1) y.axis.off=FALSE
      else y.axis.off <- TRUE
    plot.list.c[[i]] <- plot_coverage(Z[[i]], main=main, proportion=0.2,
                                      cols=c("black", "forestgreen", "deeppink3", "gold4", "red"),
                                      shapes= c(1, 3, 0, 2, 8),
                                      simnames=c("naive", "wfb2", "par", "wfb", "nonpar"),
                                      legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap"),
                                      y.axis.off=y.axis.off, legend.position = "none", set.range=c(mincov, 1))
  }
    plot.list.w <- list()
    for(i in 1:2){
      main =  ""
      if(i == 1) y.axis.off=FALSE
        else y.axis.off = TRUE
      plot.list.w[[i]] <- plot_width(Z[[1]], main=main, proportion=0.2,
                                     cols=c("black", "forestgreen", "deeppink3", "gold4", "red"),
                                     shapes= c(1, 3, 0, 2, 8),
                                     simnames=c("naive", "wfb2", "par", "wfb", "nonpar"),
                                     legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap"),
                                     y.axis.off=y.axis.off, legend.position = "none")
    }
    h <- arrangeGrob(plot.list.c[[1]], plot.list.c[[2]],
                     plot.list.w[[1]], plot.list.w[[2]],
                     ncol = 2)
    ggsave(filename = "corr_fig.ps", plot=h, dpi=600, height=8, width=8)

}
