
make_fig1 <- function(){
  titles <- c("(1)", "(2)", "(3)", "(4)")
  plot.list.c <- list()
  for(i in 1:4){
    main =  titles[i]
    Z <- sim.list[[i]]
    if(i ==1) y.axis.off=FALSE
      else y.axis.off <- TRUE
    plot.list.c[[i]] <- plot_coverage(Z, main=main, proportion=0.2,
                                              cols=c("black", "forestgreen", "deeppink3", "gold4", "blue"),
                                              shapes= c(1, 3, 0, 2, 4),
                                              simnames=c("naive", "wfb2", "par", "wfb", "oracle"),
                                              legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Oracle"),
                                        y.axis.off=y.axis.off, legend.position = "none")
  }
  plot.list.w <- list()
  for(i in 1:4){
    Z <- sim.list[[i]]
    main =  ""
    if(i == 1) y.axis.off=FALSE
      else y.axis.off = TRUE
    plot.list.w[[i]] <- plot_width(Z, main=main, proportion=0.2,
                                      cols=c("black", "forestgreen", "deeppink3", "gold4", "blue"),
                                      shapes= c(1, 3, 0, 2, 4),
                                      simnames=c("naive", "wfb2", "par", "wfb", "oracle"),
                                      legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Oracle"),
                                      y.axis.off=y.axis.off, legend.position = "none")
  }
  h <- arrangeGrob(plot.list.c[[1]], plot.list.c[[2]], plot.list.c[[3]], plot.list.c[[4]],
                 plot.list.w[[1]], plot.list.w[[2]], plot.list.w[[3]], plot.list.w[[4]],
                 ncol = 4)
  ggsave(filename = "example_fig.ps", plot=h, dpi=600, height=8, width=16)
}


make_fig1_si <- function(){
  sim.list <- rcc:::getobj("sim_results/sims.RData")
  titles <- c("All zero", "All small, non-zero", "A few large, non-zero", "A few small, non-zero")
  plot.list.c <- list()
  for(i in 1:4){
    main =  titles[i]
    Z <- sim.list[[i]]
    if(i ==1){
      y.axis.off=FALSE
      #lp <- c(0.28, 0.4)
      lp <- "none"
    }else{
      y.axis.off <- TRUE
      lp <- "none"
    }
    plot.list.c[[i]] <- plot_coverage(Z, main=main, proportion=0.2,
                                      cols=c("black", "forestgreen", "deeppink3","blue", "darkorchid", "gold4"),
                                      shapes= c(1, 3, 0, 3, 2, 4), span=0.2, y.range=c(-0.05, 1.05),
                                      ltys=c( 2, 4, 1, 3, 6, 5),
                                      simnames=c(  "naive", "wfb2", "par", "oracle", "selInf1", "wfb"),
                                      legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "Oracle", "SelInf", "WFB"),
                                      y.axis.off=y.axis.off, legend.position = lp)
  }
  plot.list.w <- list()
  for(i in 1:4){
    Z <- sim.list[[i]]
    main =  ""
    if(i == 1) y.axis.off=FALSE
    else y.axis.off = TRUE
    plot.list.w[[i]] <- plot_width(Z, main=main, proportion=0.2,
                                   cols=c("black", "forestgreen", "deeppink3","blue", "darkorchid", "gold4"),
                                   shapes= c(1, 3, 0, 3, 2, 4), span=0.2,
                                   ltys=c( 2, 4, 1, 3, 6, 5),
                                   simnames=c(  "naive", "wfb2", "par", "oracle", "selInf1", "wfb"),
                                   legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "Oracle", "SelInf", "WFB"),
                                   y.axis.off=y.axis.off, legend.position = "none")
  }
  h <- arrangeGrob(plot.list.c[[1]], plot.list.c[[2]], plot.list.c[[3]], plot.list.c[[4]],
                   plot.list.w[[1]], plot.list.w[[2]], plot.list.w[[3]], plot.list.w[[4]],
                   ncol = 4)
  ggsave(filename = "example_fig2.pdf", plot=h, dpi=600, height=8, width=16)
  return(h)
}



