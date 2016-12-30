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
  data("population_genotypes", package="rcc")
  X.pop <- population_genotypes
  which.type <- which(c("none", "ten_equal", "ten_norm", "mixture") == beta.type)
  stopifnot(length(which.type)==1)
  data("linreg_params", package="rcc")
  
  z <- linreg_sim(X.pop=X.pop, which.sample = 1:100, n.rep=1,
                  index=linreg_params[[which.type]]$index, 
                  beta=linreg_params[[which.type]]$effect, seed=seed)
  save(z, file=paste0(file.prefix,  "_", beta.type,  "_n", which.run, ".RData"))
}
