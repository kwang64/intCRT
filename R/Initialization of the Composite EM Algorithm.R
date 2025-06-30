### Initializes the coefficient and baseline hazard terms
em_init <- function(p, tau){

  # Setting the coefficient vector to zero and the baseline hazard function to 1/rho.s for each stratum
  beta.init <- rep(0, p)
  lambda.init <- lapply(tau, FUN=function(x)rep(1/length(x), length(x)))

  return(list(beta=beta.init, lambda=lambda.init))
}
