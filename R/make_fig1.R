
make_fig1 <- function(){
  titles <- c("(1)", "(2)", "(3)", "(4)")
  plot.list.c <- list()
  for(i in 1:4){
    main =  titles[i]
    Z <- sim.list[[i]]
    if(i ==1) y.axis.off=FALSE
      else y.axis.off <- TRUE
    plot.list.c[[i]] <- plot_example_coverage(Z, main=main, proportion=0.2,
                                        y.axis.off=y.axis.off, legend.position = "none")
  }
  plot.list.w <- list()
  for(i in 1:4){
    Z <- sim.list[[i]]
    main =  ""
    if(i == 1) y.axis.off=FALSE
      else y.axis.off = TRUE
    plot.list.w[[i]] <- plot_example_width(Z, main=main, proportion=0.2,
                                    legend.position = "none", y.axis.off = y.axis.off)
  }
  h <- arrangeGrob(plot.list.c[[1]], plot.list.c[[2]], plot.list.c[[3]], plot.list.c[[4]],
                 plot.list.w[[1]], plot.list.w[[2]], plot.list.w[[3]], plot.list.w[[4]],
                 ncol = 4)
  ggsave(filename = "example_fig.ps", plot=h, dpi=600, height=8, width=16)
}
