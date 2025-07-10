#' Solve the Profile Composite Likelihood Score Equation for \code{beta}
#'
#' Updates the fixed effects \code{beta} by solving the profile composite likelihood
#' score equation using a root-finding procedure, with optional perturbation weights.
#'
#' @param beta Numeric vector of length \code{p}, the current estimates of fixed effects.
#' @param u.star List of length \code{S}, each element a numeric vector of \code{u^*_sv} for stratum \code{s}.
#' @param tau List of length \code{S}, each element a numeric vector of unique timepoints \code{tau.s} for stratum \code{s}.
#' @param x List of length \code{S}, each element a 3D covariate array \code{x.s} for stratum \code{s}.
#' @param w List of length \code{S}, each element a matrix of expected values from the E-step \code{w.s} for stratum \code{s}.
#' @param risk List of length \code{S}, each element a matrix of risk set indicators \code{risk.vr} for stratum \code{s}.
#' @param weights Optional list of length \code{S}, each element a numeric vector of weights \code{weights.s} for stratum \code{s}.
#' If \code{NULL}, defaults to equal weights for all individuals.
#'
#' @return A numeric vector of updated \code{beta} values (length \code{p}).
#' @export
update_beta <- function(beta, u.star, tau, x, w, risk, weights=NULL){

  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }

  # Finding and returning zeroes of the profile composite score equation for beta
  suppressWarnings(rootSolve::multiroot(beta_score, start=beta, u.star=u.star, tau=tau, x=x, w=w, risk=risk, weights=weights, maxiter=1)$root, classes="warning")

}

#' Evaluate Profile Composite Likelihood Score for All Strata
#'
#' Computes the composite score function for the fixed effects \code{beta} across all strata,
#' by summing the stratum-level profile composite likelihood score contributions.
#'
#' @param beta Numeric vector of length \code{p}, the current estimates of fixed effects.
#' @param u.star List of length \code{S}, where each element is a numeric vector of \code{u^*_sv} values for stratum \code{s}.
#' @param tau List of length \code{S}, where each element is a numeric vector of unique timepoints for stratum \code{s}.
#' @param x List of length \code{S}, where each element is a 3D array of covariates \code{x.s} for stratum \code{s}.
#' @param w List of length \code{S}, where each element is a matrix \code{w.s} containing expected counts from the E-step.
#' @param risk List of length \code{S}, where each element is a matrix indicating risk set membership at each timepoint in stratum \code{s}.
#' @param weights List of length \code{S}, where each element is a numeric vector of perturbation weights for stratum \code{s}.
#'
#' @return A numeric vector of length \code{p}, representing the total score contributions for each covariate across all strata.
#'
#' @export
#'
#' @keywords internal

beta_score <- function(beta, u.star, tau, x, w, risk, weights){

  #' Evaluate Profile Composite Likelihood Score for a Single Stratum
  #'
  #' Computes the score function (gradient of the composite log-likelihood)
  #' for the profile composite likelihood in a single stratum.
  #'
  #' @param beta Numeric. A \code{p}-dimensional vector of current fixed effect estimates.
  #' @param u.s.star Numeric vector of length \code{n.s}, representing the right endpoints (\code{u^*_sv}) for individuals in stratum \code{s}.
  #' @param tau.s Numeric vector of length \code{rho.s}, giving the sorted unique left and right endpoints for stratum \code{s}.
  #' @param x.s 3D array. Covariate array of dimension \code{n.s x p x rho.s}, where:
  #'   \itemize{
  #'     \item \code{x.s[v, , ]} returns a \code{p x rho.s} matrix for individual \code{v}
  #'     \item \code{x.s[, p, ]} returns an \code{n.s x rho.s} matrix for covariate \code{p}
  #'     \item \code{x.s[, , r]} returns an \code{n.s x p} matrix for timepoint \code{r}
  #'   }
  #' @param w.s Matrix of dimension \code{n.s x rho.s}, containing the expected counts from the E-step (i.e., \code{w_{svr}}).
  #' @param risk.vr Matrix of dimension \code{n.s x rho.s}, indicating risk set membership at each timepoint.
  #' @param weights.s Numeric vector of length \code{n.s}, giving the perturbation weights for individuals in stratum \code{s}.
  #'
  #' @return A numeric vector of length \code{p}, containing the score contributions for each covariate in the stratum.
  #'
  #' @keywords internal
  #' @seealso \code{\link{x_bar_s}}, \code{\link{score}}, \code{\link{update_beta}}
  #' @export
  score_s <- function(beta, u.s.star, tau.s, x.s, w.s, risk.vr, weights.s){

    # Determining the weight matrix for each individual at each timepoint
    weight.vr <- risk.vr * w.s * weights.s

    # Calculating stratum-level score contribution for each covariate
    x.bar <- x_bar_s(beta, u.s.star, tau.s, x.s, risk.vr, weights.s)
    list.x.p <- lapply(seq_len(dim(x.s)[2]), function(x) x.s[, x, ])
    list.x.bar <- lapply(seq_len(ncol(x.bar)), function(x)x.bar[, x])
    score_out <- mapply(function(x.vr, x.bar){
      sum(rowSums(weight.vr * t(apply(x.vr, 1, FUN=function(x)x-x.bar))))
    }, x.vr=list.x.p, x.bar=list.x.bar)

    return(score_out)
  }

  # Calculating and aggregating the score contributions across strata
  rowSums(do.call(cbind, mapply(function(u.s.star, tau.s, x.s, w.s, risk.vr, weights.s) score_s(beta, u.s.star, tau.s, x.s, w.s, risk.vr, weights.s),
                                u.s.star=u.star, tau.s=tau, x.s=x, w.s=w, risk.vr=risk, weights.s=weights, SIMPLIFY=FALSE)))
}

#' Computes the mean covariate vector at each timepoint in stratum s
#'
#' This function calculates the weighted average of covariates at each timepoint \code{tau.s}
#' in a given stratum \code{s}, using the current estimate of the regression coefficients \code{beta}.
#'
#' @param beta Numeric vector of length \code{p}, the current fixed effect estimates.
#' @param u.s.star Numeric vector of length \code{n.s}, containing the right endpoints for each subject in stratum \code{s}.
#' @param tau.s Numeric vector of length \code{rho.s}, the sorted unique timepoints (left and right endpoints) for stratum \code{s}.
#' @param x.s Array of dimension \code{n.s x p x rho.s}, the covariate array for stratum \code{s}.
#' @param risk.vr Numeric matrix of dimension \code{n.s x rho.s}, indicating risk set membership for each subject at each timepoint.
#' @param weights.s Numeric vector of length \code{n.s}, containing perturbation weights for each individual in stratum \code{s}.
#'
#' @return A matrix of dimension \code{rho.s x p}, where each row gives the mean covariate vector at a specific timepoint.
#'
#' @keywords internal
#' @export
x_bar_s <- function(beta, u.s.star, tau.s, x.s, risk.vr, weights.s){

  # Returns a matrix with rows corresponding to individuals and columns to timepoints
  constant.vr <- do.call(rbind, lapply(seq_len(dim(x.s)[1]), FUN=function(x) exp(beta %*% x.s[x, , ])))

  # Returns a matrix with rows corresponding to timepoints and columns to covariates
  out <- matrix(NA, nrow=length(tau.s), ncol=length(beta))
  for (p in 1:length(beta)){
    x.vr <- x.s[ , p, ]
    out[, p] <- colSums(risk.vr * x.vr * constant.vr * weights.s)/colSums(risk.vr * constant.vr * weights.s)
  }

  return(out)
}

#' Update All Baseline Hazard Components Across Strata
#'
#' Applies the M-step update of the baseline hazard estimates \eqn{\lambda_{sr}} to all strata, using current estimates of the coefficients and expected values from the E-step.
#'
#' @param beta Numeric vector of length \eqn{p}, containing the current fixed-effect coefficient estimates.
#' @param x List of length \eqn{S}, where each element is a covariate array \code{x.s} for a stratum.
#' @param w List of length \eqn{S}, where each element is a matrix \code{w.s} of expected Poisson augmentation variables for a stratum.
#' @param risk List of length \eqn{S}, where each element is a matrix \code{risk.vr} indicating risk set membership for each stratum.
#' @param weights Optional list of length \eqn{S}, where each element is a numeric vector of perturbation weights \code{weights.s} for the corresponding stratum. If \code{NULL}, equal weights are assumed.
#'
#' @return A list of length \eqn{S}, where each element is a numeric vector of updated baseline hazard estimates \eqn{\lambda_{sr}} for stratum \eqn{s}.
#' @export
#' @keywords internal

update_lambda <- function(beta, x, w, risk, weights=NULL){

  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }

  # Updating the estimated baseline hazard function for each stratum and aggregating as a list
  mapply(function(x.s, w.s, risk.vr, weights.s){
    update_lambda_s(beta, x.s, w.s, risk.vr, weights.s)
  }, x.s=x, w.s=w, risk.vr=risk, weights.s=weights, SIMPLIFY=FALSE)

}

#' Update Baseline Hazard Components for a Stratum
#'
#' Computes the updated baseline hazard estimates \eqn{\lambda_{sr}} at each timepoint \eqn{r} for a given stratum \eqn{s}, as part of the M-step in the composite EM algorithm.
#'
#' @param beta Numeric vector of length \eqn{p}, representing the current estimates of the fixed-effect coefficients.
#' @param x.s 3-dimensional array of covariates, where \code{x.s[v, , ]} returns a \eqn{p \times \rho_s} matrix for individual \eqn{v}; \code{x.s[, p, ]} gives an \eqn{n_s \times \rho_s} matrix for covariate \eqn{p}; and \code{x.s[, , r]} gives an \eqn{n_s \times p} matrix for timepoint \eqn{r}.
#' @param w.s Matrix (\eqn{n_s \times \rho_s}) of expected Poisson augmentation variables from the E-step, for stratum \eqn{s}.
#' @param risk.vr Matrix (\eqn{n_s \times \rho_s}) indicating risk set membership for each individual at each timepoint in stratum \eqn{s}.
#' @param weights.s Numeric vector of length \eqn{n_s}, containing perturbation weights for each individual in stratum \eqn{s}.
#'
#' @return A numeric vector of length \eqn{\rho_s}, containing the updated baseline hazard estimates \eqn{\lambda_{sr}} for stratum \eqn{s}.
#' @export
#' @keywords internal

update_lambda_s <- function(beta, x.s, w.s, risk.vr, weights.s){

  # Calculating the denominator at each timepoint
  constant.vr <- do.call(rbind, lapply(seq_len(dim(x.s)[1]), FUN=function(x) exp(beta %*% x.s[x, , ])))
  denom <- colSums(risk.vr * constant.vr * weights.s)

  # Calculating the numerator at each timepoint
  num <- colSums(risk.vr * w.s * weights.s)

  return(num/denom)

}
