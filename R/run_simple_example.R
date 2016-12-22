
run_simple_example <- function(){
  set.seed(6587900)
  sim.list <- list()
  for(i in 1:4){
    sim.list[[i]] <- example_sim(example_params[,i], n=200, use.abs=FALSE)
  }
  return(sim.list)
  #save(sim.list, file="data/example_sims.RData")
}
