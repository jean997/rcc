#' Run one linear regression simulation using builtin parameters.
#' @description This function is a wrapper for linreg_sim and calls it using built-in data
#' `macs_params`. `macs_params` is a list containing the effect sizes and variables
#' used to create the simulations in section 3.1. This function writes data to an output
#' file and is intended to be run in parallell.
#' @param beta.type One of "none", "ten_equal", "ten_norm", or "mixture"
#' @param which.run An integer that will get included in the name of the output file.
#' @param file.prefix Simulations results will be saved to file.prefix_`beta.type`_n`which.run`.RData
#' @param seed Optional seed to set
#'
#' @return This function saves an object produced by `linreg_sims` to a file.
#'@export
run1_linreg_sim <- function(beta.type, which.run, file.prefix, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  eff <- 2
  eff.sd <- 1

  n=100

  X.pop <- population_genotypes

  p <- 500

  which.type <- which(c("none", "ten_equal", "ten_norm", "mixture") == beta.type)
  stopifnot(length(which.type)==1)
  data("macs_params", package="rcc")
  which.X0 <- macs_params$include[[which.type]]
  beta <- macs_params$betas[[which.type]]

  z <- linreg_sim(X.pop=X.pop, which.sample = 1:n, n.rep=1,
                  which.X0=which.X0, beta=beta)
  save(z, file=paste0(file.prefix,  "_", beta.type,  "_n", which.run, ".RData"))
}
