
make_fig2 <- function(){
  titles <- c("No Effects", "Ten Large Effects", "Ten Small Dispersed Effects", "Combination of Effects")
  beta.type <- c("none", "ten_equal", "ten_norm", "mixture")
  plot.list.c <- list()
  for(i in 1:4){
    main =  titles[i]
    fn <- paste0("sim_results/linreg_", beta.type[i], ".RData")
    Z <- rcc:::getobj(fn)
    if(i ==1) y.axis.off=FALSE
      else y.axis.off <- TRUE
    plot.list.c[[i]] <- plot_coverage(Z, main=main, proportion=0.2,
                                              cols=c("black", "forestgreen", "deeppink3", "red", "darkorchid", "gold4"),
                                              shapes= c(1, 3, 0, 2, 8, 4), span=0.2, y.range=c(-0.05, 1.05),                                             
                                              ltys= c( 2, 4, 1, 3, 6, 5),
                                              simnames=c("naive", "wfb2", "par",  "nonpar", "selInf1", "wfb"),
                                              legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap", "Selective Inference"),
                                        y.axis.off=y.axis.off, legend.position = "none")
  }
  plot.list.w <- list()
  for(i in 1:4){
    if(i ==1 | i ==3) y.max=2
      else if(i==2) y.max=3.2
        else y.max=NULL
    fn <- paste0("sim_results/linreg_", beta.type[i], ".RData")
    Z <- rcc:::getobj(fn)
    main =  titles[i]
    if(i == 1) y.axis.off=FALSE
      else y.axis.off = TRUE
    plot.list.w[[i]] <- plot_width(Z, main=main, proportion=0.2,
                                      cols=c("black", "forestgreen", "deeppink3", "red", "darkorchid", "gold4"),
                                      shapes= c(1, 3, 0, 2, 8, 4), span=0.2,
                                      ltys= c( 2, 4, 1,3,6, 5), y.max=y.max,
                                      simnames=c("naive", "wfb2", "par", "nonpar", "selInf1", "wfb"),
                                      legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Non-Parametric Bootstrap", "Selective Inference"),
                                      y.axis.off=y.axis.off, legend.position = "none")
  }
  return(list("covereage"=plot.list.c, "width"=plot.list.w))
}
