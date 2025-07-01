#' Perform a Single Iteration of the Composite EM Algorithm
#'
#' Executes one iteration of the composite Expectation-Maximization (EM) algorithm
#' for estimating the regression parameters and baseline hazard functions in marginal
#' Cox regression models with clustered interval-censored data.
#'
#' @param beta Numeric vector. Current estimates of the fixed effect coefficients, of length \code{p}.
#' @param lambda List. Current estimates of the baseline hazard functions, one vector \code{lambda.s} for each stratum.
#' @param l List. Left endpoints \code{l.s} of the censoring intervals for each stratum.
#' @param u List. Right endpoints \code{u.s} of the censoring intervals for each stratum.
#' @param u.star List. Modified upper bounds \code{u.s.star} with \code{Inf} replaced for computational convenience.
#' @param tau List. Sorted unique time points \code{tau.s} for each stratum.
#' @param x List. Covariate arrays \code{x.s} for each stratum.
#' @param risk List. Binary matrices \code{risk.vr} indicating risk set membership for each observation and time point.
#' @param weights Optional list. Perturbation weights \code{weight.s} for each stratum, used in bootstrap variance estimation.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{beta}: Updated coefficient estimates after one EM iteration.
#'   \item \code{lambda}: Updated baseline hazard estimates after one EM iteration.
#' }
#'
#' @seealso \code{\link{composite_coxIC}}, \code{\link{em_init}}, \code{\link{update_beta}}, \code{\link{update_lambda}}
#'
#' @export

point_iter <- function(beta, lambda, l, u, u.star, tau, x, risk, weights=NULL){

  # Expectation step
  w <- calc_w(l, u, tau, x, lambda, beta)

  # Maximization step
  beta.iter <- update_beta(beta, u.star, tau, x, w, risk, weights)
  lambda.iter <- update_lambda(beta.iter, x, w, risk, weights)

  return(list(beta=beta.iter, lambda=lambda.iter))
}
