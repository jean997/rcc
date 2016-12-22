

make_fig3 <- function(){
  titles <- c("Low Correlation", "High Correlation")
  plot.list.c <- list()
  Z <- list()
  Z[[1]] <- rcc:::getobj("sim_results/chr1_sim_results.RData")
  Z[[2]] <- rcc:::getobj("sim_results/dnarep_sim_results.RData")

  p <- 4950
  k <- floor(0.8*p)
  mincov <- min(rowMeans(Z[[1]]$COVERAGE[1, , ])[k:p], na.rm=TRUE)
  for(i in 2:6){
    mincov <- min(mincov, min(rowMeans(Z[[1]]$COVERAGE[i, , ])[k:p], na.rm=TRUE))
  }
  for(i in 1:2){

    main =  titles[i]

    if(i ==1) y.axis.off=FALSE
      else y.axis.off <- TRUE
    plot.list.c[[i]] <- plot_coverage(Z[[i]], main=main, proportion=0.2,
                                      cols=c("black", "forestgreen", "deeppink3", "red", "darkorchid", "gold4"),
                                      shapes= c(1, 3, 0, 2, 8, 4), span=0.05,
                                      ltys= c( 2, 4, 1, 3, 6, 5),
                                      simnames=c("naive", "wfb2", "par",  "nonpar", "selInf1", "wfb"),
                                      legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap", "Selective Inference"),
                                      y.axis.off=y.axis.off, legend.position = "none", y.range=c(mincov, 1))
  }
    plot.list.w <- list()
    for(i in 1:2){
      main =  titles[i]
      if(i == 1) y.axis.off=FALSE
        else y.axis.off = TRUE
      plot.list.w[[i]] <- plot_width(Z[[1]], main=main, proportion=0.2,
                                     cols=c("black", "forestgreen", "deeppink3", "red", "darkorchid", "gold4"),
                                     shapes= c(1, 3, 0, 2, 8, 4), span=0.05,
                                     ltys= c( 2, 4, 1, 3, 6, 5),
                                     simnames=c("naive", "wfb2", "par",  "nonpar", "selInf1", "wfb"),
                                     legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap", "Selective Inference"),
                                     y.axis.off=y.axis.off, legend.position = "none")
    }
    return(list("covereage"=plot.list.c, "width"=plot.list.w))
}
