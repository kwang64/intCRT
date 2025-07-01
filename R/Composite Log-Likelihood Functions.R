#' Evaluates the composite log-likelihood function
#'
#' Computes the independence composite log-likelihood under a marginal Cox regression model for clustered interval-censored data.
#'
#' @param beta Numeric vector. Current estimates of the fixed effect coefficients (length \eqn{p}).
#' @param lambda List of length \eqn{S}. Each element is a numeric vector \eqn{\lambda_s} representing the estimated baseline hazard function for stratum \eqn{s}.
#' @param l List of length \eqn{S}. Each element \code{l[[s]]} is a numeric vector of left endpoints \eqn{L_{sv}} for stratum \eqn{s}.
#' @param u List of length \eqn{S}. Each element \code{u[[s]]} is a numeric vector of right endpoints \eqn{U_{sv}} for stratum \eqn{s}.
#' @param tau List of length \eqn{S}. Each element \code{tau[[s]]} is a numeric vector of sorted unique timepoints \eqn{\tau_{sr}} for stratum \eqn{s}.
#' @param x List of length \eqn{S}. Each element \code{x[[s]]} is a 3D array of covariates of dimension \eqn{n_s \times p \times \rho_s} for stratum \eqn{s}.
#' @param mapping Optional list used when clustering and stratification differ. Should contain:
#' \describe{
#'   \item{\code{stratum}}{List of length \eqn{M}. Each element maps subject \eqn{j} in cluster \eqn{i} to their corresponding stratum index \eqn{s}.}
#'   \item{\code{indiv}}{List of length \eqn{M}. Each element maps subject \eqn{j} in cluster \eqn{i} to their subject index \eqn{v} within stratum \eqn{s}.}
#' }
#'
#' @return Numeric. The overall composite log-likelihood across all strata or clusters.
#' @export
composite <- function(beta, lambda, l, u, tau, x, mapping=NULL){

  # If no mapping is provided, the clustering and stratification variable(s) is(are) assumed to be the same
  if (is.null(mapping)){
    do.call(sum, mapply(FUN=function(lambda.i, l.i, u.i, tau.i, x.i) composite_i_1(beta, lambda.i, l.i, u.i, tau.i, x.i),
                        lambda.i=lambda, l.i=l, u.i=u, tau.i=tau, x.i=x, SIMPLIFY=FALSE))

    # Otherwise, the mapping between cluster and stratum memberships are used to construct the composite log-likelihood
  } else{
    do.call(sum, mapply(FUN=function(stratum.i, indiv.i) composite_i_2(beta, lambda, l, u, tau, x, stratum.i, indiv.i),
                        stratum.i=mapping$stratum, indiv.i=mapping$indiv, SIMPLIFY=FALSE))
  }
}

#' Computes the composite log-likelihood contribution for a cluster (when clustering and stratification coincide)
#'
#' Calculates the independence composite log-likelihood contribution of cluster \eqn{i}, assuming the clustering and stratification variable(s) are the same,
#' under a marginal Cox regression model for interval-censored data.
#'
#' @param beta Numeric vector. Current estimates of the fixed effect coefficients (length \eqn{p}).
#' @param lambda.i Numeric vector of length \eqn{\rho_i}. Estimated baseline hazard function for cluster (stratum) \eqn{i}.
#' @param l.i Numeric vector of length \eqn{n_i}. Left endpoints \eqn{L_{ij}} of censoring intervals for individuals in cluster \eqn{i}.
#' @param u.i Numeric vector of length \eqn{n_i}. Right endpoints \eqn{U_{ij}} of censoring intervals for individuals in cluster \eqn{i}.
#' @param tau.i Numeric vector of length \eqn{\rho_i}. Sorted unique timepoints \eqn{\tau_{ir}} for cluster (stratum) \eqn{i}.
#' @param x.i Array of dimension \eqn{n_i \times p \times \rho_i}. Covariate array where \code{x.i[j, , ]} gives a \eqn{p \times \rho_i} matrix for individual \eqn{j} in cluster \eqn{i}.
#'
#' @return Numeric. The composite log-likelihood contribution from cluster \eqn{i}.
#'
#' @details
#' For each individual \eqn{j} in cluster \eqn{i}, this function extracts the individual covariate matrix, identifies the appropriate interval,
#' and evaluates the individual composite log-likelihood contribution using \code{\link{composite_ij}}. The total contribution from cluster \eqn{i}
#' is then calculated as the sum of these individual contributions.
#'
#' @keywords internal
#' @export
composite_i_1 <- function(beta, lambda.i, l.i, u.i, tau.i, x.i){

  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(l.ij, u.ij, x.pr) composite_ij(l.ij, u.ij, tau.i, x.pr, lambda.i, beta),
                      l.ij=l.i, u.ij=u.i, x.pr=lapply(seq_len(dim(x.i)[1]), FUN=function(x)matrix(x.i[x, , ], nrow=dim(x.i)[2])),
                      SIMPLIFY=FALSE))

}

#' Computes the composite log-likelihood contribution for a cluster (stratified and clustered data)
#'
#' Calculates the independence composite log-likelihood contribution of cluster \eqn{i}, when the clustering and stratification variables differ,
#' under a marginal Cox regression model with interval-censored data.
#'
#' @param beta Numeric vector. Current estimates of the fixed effect coefficients (length \eqn{p}).
#' @param lambda List of length \eqn{S}. Each element \eqn{lambda[[s]]} is the estimated baseline hazard vector for stratum \eqn{s}.
#' @param l List of length \eqn{S}. Each \eqn{l[[s]]} contains the left endpoints \eqn{L_{sv}} of the censoring intervals in stratum \eqn{s}.
#' @param u List of length \eqn{S}. Each \eqn{u[[s]]} contains the right endpoints \eqn{U_{sv}} of the censoring intervals in stratum \eqn{s}.
#' @param tau List of length \eqn{S}. Each \eqn{tau[[s]]} is the sorted vector of unique timepoints \eqn{\tau_{sr}} in stratum \eqn{s}.
#' @param x List of length \eqn{S}. Each \eqn{x[[s]]} is a 3D array of covariates for stratum \eqn{s}, where \eqn{x[[s]][v, , ]} returns a \eqn{p \times \rho_s} matrix.
#' @param stratum.i Numeric vector of length \eqn{n_i}. Maps subject \eqn{j} in cluster \eqn{i} to their stratum index \eqn{s}.
#' @param indiv.i Numeric vector of length \eqn{n_i}. Maps subject \eqn{j} in cluster \eqn{i} to their subject index \eqn{v} within stratum \eqn{s}.
#'
#' @return Numeric. The composite log-likelihood contribution from cluster \eqn{i}.
#'
#' @details
#' The function iterates over each individual \eqn{j} in cluster \eqn{i}, determines their associated stratum and individual index,
#' extracts the corresponding data, and evaluates the individual composite log-likelihood contribution using \code{\link{composite_ij}}.
#'
#' @keywords internal
#' @export
composite_i_2 <- function(beta, lambda, l, u, tau, x, stratum.i, indiv.i){

  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(stratum.ij, indiv.ij){
    tau.s <- tau[[stratum.ij]]; lambda.s <- lambda[[stratum.ij]]
    l.ij <- l[[stratum.ij]][indiv.ij]; u.ij <- u[[stratum.ij]][indiv.ij]
    x.pr <- adrop(x[[stratum.ij]][indiv.ij, , , drop=FALSE], drop=1)
    composite_ij(l.ij, u.ij, tau.s, x.pr, lambda.s, beta)
  }, stratum.ij=stratum.i, indiv.ij=indiv.i, SIMPLIFY=FALSE))

}
#' Computes the composite log-likelihood contribution for an individual
#'
#' Calculates the independence composite log-likelihood contribution for individual \eqn{j} in cluster \eqn{i}
#' under a marginal Cox regression model with interval-censored data.
#'
#' @param l.ij Numeric. Left endpoint of the censoring interval for individual \eqn{j} in cluster \eqn{i}.
#' @param u.ij Numeric. Right endpoint of the censoring interval for individual \eqn{j} in cluster \eqn{i}.
#' @param tau.s Numeric vector. Sorted vector of unique timepoints \eqn{\tau_{sr}} (length \eqn{\rho_s}) for stratum \eqn{s}.
#' @param x.pr Matrix (p × rho.s). Covariate values for individual \eqn{j} in cluster \eqn{i} across all timepoints.
#' @param lambda.s Numeric vector. Estimated baseline hazard components for stratum \eqn{s} at each timepoint.
#' @param beta Numeric vector. Current estimates of the fixed effect coefficients (length \eqn{p}).
#'
#' @return Numeric. The composite log-likelihood contribution of individual \eqn{j} in cluster \eqn{i}.
#' @export
#' @details The function computes two terms:
#' \enumerate{
#'   \item The cumulative hazard up to the left endpoint \eqn{L_{ij}}.
#'   \item The conditional probability of an event in the interval \eqn{(L_{ij}, U_{ij}]} if \eqn{U_{ij} < \infty}.
#' }
#'
composite_ij <- function(l.ij, u.ij, tau.s, x.pr, lambda.s, beta){

  # Computing the first term of the individual log-likelihood
  r.ids.1 <- (tau.s <= l.ij)
  term.1 <- -sum(lambda.s[r.ids.1] * exp(beta %*% x.pr[, r.ids.1]))

  # Computing the second term of the individual log-likelihood
  if (u.ij < Inf){
    r.ids.2 <- (tau.s > l.ij & tau.s <= u.ij)
    term.2 <- log(1-exp(-sum(lambda.s[r.ids.2] * exp(beta %*% x.pr[, r.ids.2]))))
  } else{
    term.2 <- 0
  }

  # Individual composite log-likelihood contribution
  term.1 + term.2
}
