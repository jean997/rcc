

create_macs_params <- function(){
  betas <- list()
  include <- list()
  #No effects
  betas[[1]] <- c(0)
  include[[1]] <- c()
  #10 large effects
  betas[[2]] <- c(0, rep(2, 10))
  include[[2]] <- seq(from=1, to=500, by=50)
  #10 normal effects
  set.seed(100)
  betas[[3]]<- c(0, rnorm(10, mean=0, sd=1))
  include[[3]] <-  seq(from=1, to=500, by=50)
  #Mixture of effects
  betas[[4]] <- rnorm(500, mean=0, sd=0.25)
  big.effects <- seq(from=1, to=500, by=50)
  betas[[4]][big.effects] <- rnorm(10, mean=2, sd=1)
  include[[4]] <- 1:500
  return(list("betas"=betas, "include"=include))
}
