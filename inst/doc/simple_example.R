## ---- echo = FALSE, warning=FALSE----------------------------------------
library(rcc)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- show='hold'--------------------------------------------------------
titles <- c("All Zero", "All N(0, 1)", "100 Effects=3", "100 Effects N(0, 1)")
for(i in 1:4){
  hist(example_params[,i], breaks=15, main=titles[i], xlab="Theta")
}

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(6587900)
#  sim.list <- list()
#  for(i in 1:4){
#    sim.list[[i]] <- example_sim(example_params[,i], n=200, use.abs=FALSE)
#  }

## ------------------------------------------------------------------------
names(sim.list[[1]])
sim.list[[1]]$simnames
dim(sim.list[[1]]$COVERAGE)
dim(sim.list[[1]]$WIDTH)

## ---- fig.show='hold', fig.width=6.5-------------------------------------
plot_coverage(sim.list[[2]], proportion=0.2,
              cols=c("black", "forestgreen", "deeppink3", "gold4", "blue"),
              shapes= c(1, 3, 0, 2, 4),
              simnames=c("naive", "wfb2", "par", "wfb", "oracle"),
              legend.names = c("Marginal", "WFB-Sliding", "Parametric Bootstrap", "WFB", "Oracle"), legend.position = "right")

