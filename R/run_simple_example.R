
run_simple_example <- function(outfile.prefix){
  set.seed(6587900)
  for(i in 1:4){
    fn <- paste(outfile.prefix, ".i.RData")
    Z <- example_sim(DELTA[,i], n=200, use.abs=FALSE)
    save(Z, file=fn)
  }
}
