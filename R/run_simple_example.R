
run_simple_example <- function(){
  set.seed(6587900)
  sim.list <- list()
  for(i in 1:4){
    sim.list[[i]] <- example_sim(DELTA[,i], n=200, use.abs=FALSE)
  }
  save(sim.list, file="data/example_sims.RData")
}
