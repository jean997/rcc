

create_linreg_params <- function(){
  set.seed(100)
  params <- list()
  #No effects
  params[[1]] <- data.frame(matrix(nrow=0, ncol=2))
  names(params[[1]]) <- c("index", "effect")
  #10 large effects
  params[[2]] <- data.frame(cbind(seq(from=1, to=500, by=50), rep(2, 10)))
  names(params[[2]]) <- names(params[[1]])
  #10 normal effects
  params[[3]]<- data.frame(cbind(seq(from=1, to=500, by=50), rnorm(10, mean=0, sd=1)))
  names(params[[3]]) <- names(params[[1]])
  #Mixture of effects
  params[[4]] <- data.frame(cbind(1:500, rnorm(500, mean=0, sd=0.25)))
  big.effects <- seq(from=1, to=500, by=50)
  params[[4]][big.effects, 2] <-  rnorm(10, mean=2, sd=1)
  names(params[[4]]) <- names(params[[1]])
  return(params)
}
