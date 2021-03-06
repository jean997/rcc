% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/par_bs_ci.R
\name{par_bs_ci}
\alias{par_bs_ci}
\title{Parametric bootstrapped confidence intervals to control RCC}
\usage{
par_bs_ci(beta, se = rep(1, length(beta)), rank.func = NULL, theta = beta,
  level = 0.9, n.rep = 1000, use.abs = TRUE, ...)
}
\arguments{
\item{beta}{Parameter estimates}

\item{se}{Estimated standard error of beta. Defaults to 1.}

\item{rank.func}{A function that takes as first argument
the t-statistics beta/se and returns a list with items order and rank.
See rcc:::basic_rank for an example. If NULL, the basic_rank function 
will be used which ranks based on the size of the
test statistics.}

\item{theta}{Possibly shrunken estimates of E[beta]. Defaults to beta.}

\item{level}{Confidence level}

\item{n.rep}{Number of bootstrap replications}

\item{use.abs}{Base the rank on abs(beta) rather than beta}

\item{...}{Additional parameters to pass to rank.func}
}
\value{
A data frame giving original estimates and standard errors, 
confidence intervals, debiased point estimates, and rank for each parameter.
}
\description{
This function implements the parametric bootstrap 
(see Section 2.3 of the referenced paper). The user supplies point estimates, standard
errors and optionally, a ranking function.
}
\examples{
#generate 100 fake parameter estimates
theta <- c(rep(0, 90), rnorm(n=10)) #vector of means
beta <- rnorm(n=100, mean=theta, sd=1)
cis <- par_bs_ci(beta=beta, n.rep=500) #calculate parametric bootstrap confidence intervals
head(cis)
}
