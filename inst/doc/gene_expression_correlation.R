## ---- echo = FALSE, warning=FALSE----------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(rcc)

## ------------------------------------------------------------------------
Sigma <- cov(t(chr1_genes))
dim(Sigma)

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(19519)
#  chr1_sim_results <- cor_sim(Sigma=Sigma, n=30, n.rep=200)

## ------------------------------------------------------------------------
Sigma <- cov(t(dna_rep_genes))
dim(Sigma)

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(19500)
#  dna_rep_sim_results <- cor_sim(Sigma=Sigma, n=30, n.rep=200)

