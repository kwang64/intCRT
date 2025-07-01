#' Calculate Conditional Expectations of Poisson Augmentation Variables
#'
#' Computes the conditional expectations of the latent Poisson augmentation variables \eqn{W_{svr}}
#' used in the E-step of the composite EM algorithm for marginal Cox regression.
#'
#' @param l List. Left endpoints \code{l.s} for each stratum.
#' @param u List. Right endpoints \code{u.s} for each stratum.
#' @param tau List. Unique time points \code{tau.s} for each stratum.
#' @param x List. Covariate arrays \code{x.s} for each stratum.
#' @param lambda List. Current estimates of the baseline hazards \code{lambda.s} for each stratum.
#' @param beta Numeric vector. Current fixed effect estimates.
#'
#' @return A list of matrices \code{w.s}, one for each stratum, where each matrix contains
#' the conditional expectations \eqn{E[W_{svr} | data]}.
#'
#' @seealso \code{\link{w_s}}, \code{\link{point_iter}}
#' @export
calc_w <- function(l, u, tau, x, lambda, beta){

  # Taking the projection of w.s for each stratum
  mapply(function(l.s, u.s, tau.s, x.s, lambda.s){
    w_s(l.s, u.s, tau.s, x.s, lambda.s, beta)
  }, l.s=l, u.s=u, tau.s=tau, x.s=x, lambda.s=lambda, SIMPLIFY=FALSE)

}

#' Compute Projected Poisson Expectations for One Stratum
#'
#' For a given stratum, computes a matrix of conditional expectations \eqn{w_{svr}}
#' for each subject and time point.
#'
#' @param l.s Numeric vector. Left endpoints for individuals in the stratum.
#' @param u.s Numeric vector. Right endpoints for individuals in the stratum.
#' @param tau.s Numeric vector. Time points in the stratum.
#' @param x.s Array. Covariate data \code{x.s} for the stratum; dimensions are individuals × covariates × time points.
#' @param lambda.s Numeric vector. Current baseline hazard estimates for the stratum.
#' @param beta Numeric vector. Current fixed effect estimates.
#'
#' @return A matrix with one row per individual and one column per time point \code{tau.s},
#' containing the expected values \eqn{E[W_{svr} | data]}.
#'
#' @seealso \code{\link{w_sv}}, \code{\link{calc_w}}
#' @export
w_s <- function(l.s, u.s, tau.s, x.s, lambda.s, beta){

  # Taking the projection of w.sv for each individual
  do.call(rbind, mapply(FUN=function(l.sv, u.sv, x.pr) w_sv(l.sv, u.sv, tau.s, x.pr, lambda.s, beta),
                        l.sv=l.s, u.sv=u.s, x.pr=lapply(seq_len(dim(x.s)[1]), FUN=function(x)matrix(x.s[x, , ], nrow=dim(x.s)[2])), SIMPLIFY=FALSE))

}

#' Compute Expected Poisson Variables for One Individual
#'
#' Computes the vector of expected latent Poisson variables \eqn{w_{svr}} for a single individual
#' in a stratum, used in the E-step of the EM algorithm.
#'
#' @param l.sv Numeric. Left endpoint of the interval for the individual.
#' @param u.sv Numeric. Right endpoint of the interval for the individual.
#' @param tau.s Numeric vector. Sorted time points for the stratum.
#' @param x.pr Matrix. Covariate matrix for the individual (covariates × time points).
#' @param lambda.s Numeric vector. Current baseline hazard estimates for the stratum.
#' @param beta Numeric vector. Current fixed effect estimates.
#'
#' @return A numeric vector of length \code{length(tau.s)} with conditional expectations \eqn{w_{svr}}.
#'
#' @seealso \code{\link{w_s}}, \code{\link{calc_w}}
#' @export
w_sv <- function(l.sv, u.sv, tau.s, x.pr, lambda.s, beta){

  # Initializing the vector of projections for all W_svr such that tau_sr <= L_sv
  w.sv <- rep(0, sum(tau.s <= l.sv))

  # Computing the conditional expectation for all W_svr such that tau_sr \in (L_sv, U_sv] for U_sv < infinity
  if (u.sv != Inf){
    r.ids <- (tau.s > l.sv & tau.s <= u.sv)
    lin.pred <- beta %*% x.pr[, r.ids]
    denom <- 1 - exp(-sum(lambda.s[r.ids] * exp(lin.pred)))
    w.sv <- c(w.sv, lambda.s[r.ids]*exp(lin.pred)/denom)
  }

  # Padding with additional zeroes for all W_svr such that tau_sr > U_sv* (to aid with computation and matrix algebra)
  w.sv <- c(w.sv, rep(0, length(tau.s) - length(w.sv)))

  return(w.sv)
}
