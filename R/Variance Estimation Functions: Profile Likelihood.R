### Constructs the kth canonical vector in R^p
make_basis <- function(k, p) replace(numeric(p), k, 1)

### Returns the estimated (either model or robust) profile composite likelihood covariance matrix for the fixed effects estimators
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
### Computes the full-data sensitivity matrix
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
### Computes the full-data variability matrix
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
### Computes the variability matrix contribution from cluster i (clustering and stratification variable(s) is(are) the same)
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
### Computes the variability matrix contribution from cluster i (clustering and stratification variable(s) differ)
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
