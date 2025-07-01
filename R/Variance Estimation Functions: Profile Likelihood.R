#' Construct the k-th Canonical Basis Vector in \eqn{\mathbb{R}^p}
#'
#' Returns the \eqn{k}-th standard basis vector (canonical vector) of dimension \eqn{p},
#' which is a vector of zeros except a 1 in position \eqn{k}.
#'
#' @param k Integer. The index of the canonical basis vector to construct (1 ≤ k ≤ p).
#' @param p Integer. The dimension of the vector space.
#'
#' @return Numeric vector of length \eqn{p} with 1 at position \eqn{k} and 0 elsewhere.
#'
#' @examples
#' make_basis(2, 5)
#' # [1] 0 1 0 0 0
#'
#' @export
make_basis <- function(k, p) replace(numeric(p), k, 1)

#' Compute the Covariance Matrix for Fixed Effect Estimates
#'
#' Calculates the estimated covariance matrix of the fixed effect estimators \code{beta}
#' based on the profile composite likelihood. Returns either the model-based covariance
#' matrix or the robust sandwich covariance matrix, depending on the data structure
#' and the presence of clustering and stratification mappings.
#'
#' @param beta Numeric vector of length \eqn{p}. Final fixed effect estimates.
#' @param processed List. Processed data object containing inputs used for point estimation, including covariates, event intervals, and risk sets.
#' @param mapping List or \code{NULL}. Optional mapping from cluster indexing to stratum indexing for robust variance estimation.
#'                If \code{NULL}, clustering and stratification variables are assumed to be the same.
#' @param control List. Control parameters for the composite EM algorithm and variance estimation, including perturbation constants and tolerance.
#'
#' @return A \eqn{p \times p} numeric covariance matrix for the fixed effect estimates.
#'
#' @details
#' This function internally computes profile maximum likelihood estimators of the baseline hazard
#' under various perturbations of \code{beta} to approximate Hessian and variability matrices using
#' finite differences. It returns:
#' \itemize{
#'   \item The inverse Hessian (model-based covariance) when data clustering matches stratification or only one stratum is present.
#'   \item The robust sandwich covariance matrix (Godambe information) when clustering and stratification differ.
#' }
#'
#' The perturbation size \code{h.n} is adaptively scaled by the sample size.
#'
#' @seealso \code{\link{argmax_profile}}, \code{\link{H}}, \code{\link{J}}
#'
#' @export


variance_beta <- function(beta, processed, mapping, control){

  # Determining the number of parameters
  p <- length(beta)
  n <- do.call(sum, lapply(processed$x, FUN=function(z)dim(z)[1]))
  h.n <- control$c/sqrt(n)

  # Computing profile composite likelihood maximizers
  lambda.list <- list(CMLE=argmax_profile(beta, processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control))
  lambda.list$first.order <- vector("list", p)
  for (j in 1:p){
    lambda.list$first.order[[j]] <- argmax_profile(beta + h.n*make_basis(j, p), processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control)
    lambda.list$second.order[[j]] <- vector("list", j)
    for (k in 1:j){
      lambda.list$second.order[[j]][[k]] <- argmax_profile(beta + h.n*make_basis(j, p) + h.n*make_basis(k, p), processed$l, processed$u, processed$u.star, processed$tau, processed$x, processed$risk, control=control)
    }
  }

  if (is.null(mapping) & length(processed$l)==1){

    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Returning the model-based covariance matrix
    return(solve(-H.hat))

  } else if (is.null(mapping) & length(processed$l) > 1){

    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Computing the covariance matrix for the profile composite likelihood score
    J.hat <- J(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Computing the inverse robust sandwich (Godambe) information
    H.hat.inv <- solve(H.hat)
    G.hat.inv <- H.hat.inv %*% J.hat %*% H.hat.inv

    # Returning the robust sandwich covariance matrix
    return(G.hat.inv)

  } else if (length(mapping$stratum)==1){

    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Returning the model-based covariance matrix
    return(solve(-H.hat))

  } else {

    # Computing the Hessian matrix
    H.hat <- H(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Computing the covariance matrix for the profile composite likelihood score
    J.hat <- J(beta, lambda.list, processed$l, processed$u, processed$tau, processed$x, mapping, p, h.n)

    # Computing the inverse robust sandwich (Godambe) information
    H.hat.inv <- solve(H.hat)
    G.hat.inv <- H.hat.inv %*% J.hat %*% H.hat.inv

    # Returning the robust sandwich covariance matrix
    return(G.hat.inv)

  }

}
#' Computes the Full-Data Sensitivity Matrix
#'
#' Calculates the sensitivity matrix (Hessian matrix) for the profile composite log-likelihood
#' with respect to the fixed effect estimates \code{beta} using finite difference approximations.
#'
#' @param beta Numeric vector of length \eqn{p}. Final fixed effect estimates.
#' @param lambda.list List. Profile maximum likelihood estimators for \code{lambda} under various first- and second-order perturbations of \code{beta}.
#' @param l List of length \eqn{S}. Each element is a numeric vector of left endpoints for each stratum.
#' @param u List of length \eqn{S}. Each element is a numeric vector of right endpoints for each stratum.
#' @param tau List of length \eqn{S}. Each element is a numeric vector of sorted unique timepoints for each stratum.
#' @param x List of length \eqn{S}. Each element is a covariate array for each stratum.
#' @param mapping List or \code{NULL}. Mapping of cluster indexing \eqn{(i,j)} to stratum indexing \eqn{(s,v)}; used for robust variance estimation. If \code{NULL}, clustering and stratification variables are assumed the same.
#' @param p Integer. Dimension of the coefficient vector \code{beta}.
#' @param h.n Numeric. Perturbation constant for finite difference approximation.
#'
#' @return A \eqn{p \times p} numeric matrix representing the full-data sensitivity (Hessian) matrix.
#'
#' @export

H <- function(beta, lambda.list, l, u, tau, x, mapping, p, h.n){

  # Initializing output matrix
  out <- matrix(NA, nrow=p, ncol=p)

  # Computing each element of the Hessian matrix
  #   Can improve this coding to avoid multiple computations of the same composite likelihood value
  for (j in 1:p){
    for (k in 1:j){
      temp <- (composite(beta, lambda.list$CMLE, l, u, tau, x, mapping) - composite(beta + h.n*make_basis(j, p), lambda.list$first.order[[j]], l, u, tau, x, mapping) - composite(beta + h.n*make_basis(k, p), lambda.list$first.order[[k]], l, u, tau, x, mapping) + composite(beta + h.n*make_basis(j, p) + h.n*make_basis(k, p), lambda.list$second.order[[j]][[k]], l, u, tau, x, mapping))/h.n^2
      out[j, k] <- out[k, j] <- temp
    }
  }

  return(out)
}
#' Computes the Full-Data Variability Matrix
#'
#' Calculates the variability matrix using cluster-specific contributions.
#' Handles both cases where the clustering and stratification variables are the same or differ.
#'
#' @param beta Numeric vector of length \eqn{p}. Final fixed effect estimates.
#' @param lambda.list List. Profile maximum likelihood estimators for the baseline hazard function under various first- and second-order perturbations of \code{beta}.
#' @param l List of length \eqn{S}. Each element is a numeric vector of left endpoints of censoring intervals for each stratum.
#' @param u List of length \eqn{S}. Each element is a numeric vector of right endpoints of censoring intervals for each stratum.
#' @param tau List of length \eqn{S}. Each element is a numeric vector of sorted unique timepoints for each stratum.
#' @param x List of length \eqn{S}. Each element is a covariate array for each stratum.
#' @param mapping List or \code{NULL}. Mapping between cluster indexing \eqn{(i,j)} and stratum indexing \eqn{(s,v)} for robust variance estimation. If \code{NULL}, assumes clustering and stratification are identical.
#' @param p Integer. Dimension of the coefficient vector \code{beta}.
#' @param h.n Numeric. Perturbation constant used in finite difference approximations.
#'
#' @return A \eqn{p \times p} numeric matrix representing the full-data variability matrix.
#'
#' @details
#' This function sums cluster-specific outer products of score vectors, using either \code{J_i_1} or \code{J_i_2} depending on the relationship between clustering and stratification variables.
#' It restructures \code{lambda.list} for compatibility with cluster-specific functions if clustering and stratification coincide.
#'
#' @seealso \code{\link{J_i_1}}, \code{\link{J_i_2}}
#'
#' @export

J <- function(beta, lambda.list, l, u, tau, x, mapping, p, h.n){

  if (is.null(mapping)){

    # Restructuring lambda.list for use with mapply
    lambda <- vector("list", length(x))
    for (i in 1:length(x)){
      temp1 <- lambda.list$CMLE[[i]]
      temp2 <- lapply(seq_len(length(lambda.list$first.order)), FUN=function(z) lambda.list$first.order[[z]][[i]])
      lambda[[i]] <- list(CMLE=temp1, first.order=temp2)
    }

    # Summing across cluster-specific outer products (cluster-specific contributions to the full-data variability matrix)
    Reduce("+", mapply(FUN=function(lambda.list.i, l.i, u.i, tau.i, x.i){
      J_i_1(beta, lambda.list.i, l.i, u.i, tau.i, x.i, p, h.n)},
      lambda.list.i=lambda, l.i=l, u.i=u, tau.i=tau, x.i=x, SIMPLIFY=FALSE))

  } else {

    # Summing across cluster-specific outer products (cluster-specific contributions to the full-data variability matrix)
    Reduce("+", mapply(FUN=function(stratum.i, indiv.i){
      J_i_2(beta, lambda.list, l, u, tau, x, stratum.i, indiv.i, p, h.n)},
      stratum.i=mapping$stratum, indiv.i=mapping$indiv, SIMPLIFY=FALSE))

  }

}
#' Cluster-Specific Variability Matrix Contribution (Same Stratification and Clustering)
#'
#' Computes the contribution to the variability matrix from cluster \eqn{i}, assuming the stratification and clustering variables are the same.
#'
#' @param beta Numeric vector of length \eqn{p}. Final fixed effect estimates.
#' @param lambda.list.i List. Contains the profile maximum likelihood estimators for the baseline hazard function of cluster \eqn{i} under:
#' \describe{
#'   \item{\code{CMLE}}{Baseline hazard under the composite maximum likelihood estimator.}
#'   \item{\code{first.order}}{List of length \eqn{p}, each containing the baseline hazard under a first-order perturbation of \code{beta}.}
#' }
#' @param l.i Numeric vector of length \eqn{n_i}. Left endpoints of the censoring intervals for individuals in cluster \eqn{i}.
#' @param u.i Numeric vector of length \eqn{n_i}. Right endpoints of the censoring intervals for individuals in cluster \eqn{i}.
#' @param tau.i Numeric vector of length \eqn{\rho_i}. Sorted unique timepoints for cluster (stratum) \eqn{i}.
#' @param x.i 3D array. Covariate array for cluster \eqn{i}, where:
#' \itemize{
#'   \item \code{x.i[j, , ]} returns a \eqn{p \times \rho_i} matrix for subject \eqn{j}.
#'   \item \code{x.i[ , p, ]} returns an \eqn{n_i \times \rho_i} matrix for covariate \eqn{p}.
#'   \item \code{x.i[ , , r]} returns an \eqn{n_i \times p} matrix for timepoint \eqn{r}.
#' }
#' @param p Integer. Dimension of the coefficient vector \code{beta}.
#' @param h.n Numeric. Perturbation constant used for finite difference approximation of the score.
#'
#' @return A \eqn{p \times p} matrix representing the outer product of the gradient (score vector) for cluster \eqn{i}.
#'
#' @details
#' This function contributes to the robust sandwich variance estimator by calculating the outer product of the score vector for a cluster under the assumption that clustering and stratification are the same. The score is approximated using finite differences.
#'
#' @seealso \code{\link{composite_i_1}}, \code{\link{make_basis}}
#'
#' @export

J_i_1 <- function(beta, lambda.list.i, l.i, u.i, tau.i, x.i, p, h.n){

  # Profile likelihood contribution under the composite maximum likelihood estimator
  pl.cmle <- composite_i_1(beta, lambda.list.i$CMLE, l.i, u.i, tau.i, x.i)

  # Computing cluster-specific gradient
  score_i <- matrix(NA, nrow=p, ncol=1)
  for (j in 1:p){
    score_i[j,] <- (composite_i_1(beta + h.n*make_basis(j, p), lambda.list.i$first.order[[j]], l.i, u.i, tau.i, x.i) - pl.cmle)/h.n
  }

  score_i %*% t(score_i)
}
#' Cluster-Specific Variability Matrix Contribution (General Mapping)
#'
#' Computes the contribution to the variability matrix from cluster \eqn{i} when clustering and stratification variables differ.
#'
#' @param beta Numeric vector of length \eqn{p}. Final fixed effect estimates.
#' @param lambda.list List. Contains profile maximum likelihood estimators for the baseline hazard functions under various perturbations of \code{beta}. Must include:
#' \describe{
#'   \item{\code{CMLE}}{List of baseline hazard functions under the composite MLE.}
#'   \item{\code{first.order}}{List of length \eqn{p}, each containing the baseline hazards under first-order perturbations.}
#' }
#' @param l List of length \eqn{S}. Each element \code{l[[s]]} is a numeric vector of left endpoints for stratum \eqn{s}.
#' @param u List of length \eqn{S}. Each element \code{u[[s]]} is a numeric vector of right endpoints for stratum \eqn{s}.
#' @param tau List of length \eqn{S}. Each element \code{tau[[s]]} is a numeric vector of timepoints for stratum \eqn{s}.
#' @param x List of length \eqn{S}. Each element \code{x[[s]]} is a 3D array of covariates for stratum \eqn{s}.
#' @param stratum.i Numeric vector of length \eqn{n_i}. Maps each individual \eqn{j} in cluster \eqn{i} to their stratum index \eqn{s}.
#' @param indiv.i Numeric vector of length \eqn{n_i}. Maps each individual \eqn{j} in cluster \eqn{i} to their subject index \eqn{v} within stratum \eqn{s}.
#' @param p Integer. Dimension of the fixed effect coefficient vector \code{beta}.
#' @param h.n Numeric. Perturbation constant used for finite difference approximation of the score.
#'
#' @return A \eqn{p \times p} matrix representing the outer product of the gradient (score vector) for cluster \eqn{i}.
#'
#' @details
#' This function computes the contribution of cluster \eqn{i} to the empirical variability matrix using finite differences to approximate the gradient of the composite log-likelihood. It is designed for use in robust variance estimation when cluster and stratum identifiers differ.
#'
#' @seealso \code{\link{composite_i_2}}, \code{\link{make_basis}}
#'
#' @export

J_i_2 <- function(beta, lambda.list, l, u, tau, x, stratum.i, indiv.i, p, h.n){

  # Profile likelihood contribution under the composite maximum likelihood estimator
  pl.cmle <- composite_i_2(beta, lambda.list$CMLE, l, u, tau, x, stratum.i, indiv.i)

  # Computing cluster-specific gradient
  score.i <- matrix(NA, nrow=p, ncol=1)
  for (j in 1:p){
    score.i[j,] <- (composite_i_2(beta + h.n*make_basis(j, p), lambda.list$first.order[[j]], l, u, tau, x, stratum.i, indiv.i) - pl.cmle)/h.n
  }

  score.i %*% t(score.i)
}
