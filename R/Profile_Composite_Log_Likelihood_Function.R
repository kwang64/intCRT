### Returns the estimated baseline hazard function(s) that maximize(s) the profile composite log-likelihood
#' Estimate Baseline Hazard Functions via Profile Composite Likelihood
#'
#' Returns the estimated baseline hazard function(s) that maximize the profile composite log-likelihood for a given fixed effect vector \eqn{\beta}.
#'
#' @param beta Numeric vector. Current fixed effect estimates of length \eqn{p}.
#' @param l List of length \eqn{S}. Each element \code{l[[s]]} is a numeric vector of left endpoints \eqn{L_{sv}} for stratum \eqn{s}.
#' @param u List of length \eqn{S}. Each element \code{u[[s]]} is a numeric vector of right endpoints \eqn{U_{sv}} for stratum \eqn{s}.
#' @param u.star List of length \eqn{S}. Each element \code{u.star[[s]]} is a numeric vector of adjusted right endpoints \eqn{U^*_{sv}} for stratum \eqn{s}.
#' @param tau List of length \eqn{S}. Each element \code{tau[[s]]} is a numeric vector of sorted unique timepoints \eqn{\tau_{sr}} for stratum \eqn{s}.
#' @param x List of length \eqn{S}. Each element \code{x[[s]]} is a 3D array of covariates of dimension \eqn{n_s \times p \times \rho_s} for stratum \eqn{s}.
#' @param risk List of length \eqn{S}. Each element is a matrix of dimension \eqn{n_s \times \rho_s} indicating risk set membership at each timepoint in stratum \eqn{s}.
#' @param mapping Optional list. Required when clustering and stratification differ. Should contain:
#' \describe{
#'   \item{\code{stratum}}{List of length \eqn{M}. Maps each observation in cluster \eqn{i} to its stratum index \eqn{s}.}
#'   \item{\code{indiv}}{List of length \eqn{M}. Maps each observation in cluster \eqn{i} to its subject index \eqn{v} within stratum \eqn{s}.}
#' }
#' @param control List of control parameters for the EM algorithm. Should include:
#' \describe{
#'   \item{\code{tol}}{Stopping criterion for convergence, based on the change in log-likelihood.}
#'   \item{\code{maxit.var}}{Maximum number of EM iterations.}
#' }
#'
#' @return A list of length \eqn{S}, where each element is a numeric vector representing the estimated baseline hazard function for one stratum.
#'
#' @details
#' This function uses a composite EM algorithm to iteratively estimate the baseline hazard function while holding the fixed effect coefficients \eqn{\beta} constant. It stops when the change in the composite log-likelihood is below \code{control$tol} or after \code{control$maxit.var} iterations.
#'
#' @examples
#' \dontrun{
#' lambda_hat <- argmax_profile(beta, l, u, u.star, tau, x, risk, control = list(tol = 1e-6, maxit.var = 50))
#' }
#'
#' @export

argmax_profile <- function(beta, l, u, u.star, tau, x, risk, mapping=NULL, control){

  # Initial values for the baseline hazard and composite likelihood
  lambda <- lapply(tau, FUN=function(x)rep(1/length(x), length(x)))
  cl <- composite(beta, lambda, l, u, tau, x, mapping)

  # Obtaining profile maximum likelihood estimators for lambda given the specified value of beta
  l1.norm <- 1; n.iter <- 0
  while (l1.norm > control$tol & n.iter < control$maxit.var){
    w <- calc_w(l, u, tau, x, lambda, beta)
    lambda.iter <- update_lambda(beta, x, w, risk)
    cl.iter <- composite(beta, lambda.iter, l, u, tau, x, mapping)
    l1.norm <- abs(cl.iter - cl)
    lambda <- lambda.iter
    cl <- cl.iter
    n.iter <- n.iter + 1
  }

  return(lambda)
}
