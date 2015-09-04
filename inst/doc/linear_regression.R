## ---- echo = FALSE, warning=FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(rcc)

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(5989615)
#  all.seeds <- floor(runif(n=400, min=1000, max=1e7))
#  
#  for(beta.type in c("none", "ten_equal", "ten_norm", "mixture")){
#      for(i in 1:400){
#        run1_linreg_sim(beta.type="none", which.run=i, file.prefix="linreg", seed=all.seeds[i])
#      }
#  }

## ---- eval=FALSE---------------------------------------------------------
#  for(beta.type in c("none", "ten_equal", "ten_norm", "mixture")){
#      COVERAGE <- array(dim=c(5, 500, 400))
#      WIDTH <- array(dim=c(5, 500, 400))
#      for(i in 1:400){
#        z <- rcc:::getobj(paste0("linreg_", beta.type, "_n", i, ".RData"))
#        COVERAGE[, , i] <- z$COVERAGE
#        WIDTH[, , i] <- z$WIDTH
#      }
#      simnames <- z$simnames
#      Z <- list("COVERAGE"=COVERAGE, "WIDTH"=WIDTH, "simnames"=simnames)
#      save(Z, file=paste0("linreg_", beta.type, ".RData"))
#  }

