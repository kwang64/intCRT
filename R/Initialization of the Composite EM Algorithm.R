#' Initialize the coefficient and the baseline hazard terms
#'@param p dimension of the coefficient vector
#'@param tau list of length S containing tau.s for each stratum
#'@return A list with two components
#'\describe{
#'  \item{beta}{A list of vector 0s with length of \code{p}, representing initialized coefficients }
#'  \item{lambda}{a list of vectors where each vector has an initial baseline hazard value (1/length of tau.s) } }
#'@export

em_init <- function(p, tau){

  # Setting the coefficient vector to zero and the baseline hazard function to 1/rho.s for each stratum
  beta.init <- rep(0, p)
  lambda.init <- lapply(tau, FUN=function(x)rep(1/length(x), length(x)))

  return(list(beta=beta.init, lambda=lambda.init))
}
