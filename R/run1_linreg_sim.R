
run1_linreg_sim <- function(beta.type, which.run, file.prefix){
  eff <- 2
  eff.sd <- 1

  n=100

  X.pop <- population_genotypes

  p <- 500

  which.type <- which(c("none", "ten_equal", "ten_norm", "mixture") == beta.type)
  stopifnot(length(which.type)==1)
  which.X0 <- macs_params$include[[which.type]]
  beta <- macs_params$betas[[which.type]]

  set.seed(5989615)
  all.seeds <- floor(runif(n=400, min=1000, max=1e7))
  set.seed(all.seeds[which.run])
  z <- linreg_sim(X.pop=X.pop, which.sample = 1:n, n.rep=1,
                  which.X0=which.X0, beta=beta)
  save(z, file=paste0(file.prefix,  beta.type,  "_n", which.run, ".RData"))
}
