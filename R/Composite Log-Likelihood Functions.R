### Evaluates the independence composite log-likelihood function
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
### Calculates the independence composite log-likelihood contribution of cluster i (clustering and stratification variable(s) is(are) the same)
composite_i_1 <- function(beta, lambda.i, l.i, u.i, tau.i, x.i){

  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(l.ij, u.ij, x.pr) composite_ij(l.ij, u.ij, tau.i, x.pr, lambda.i, beta),
                      l.ij=l.i, u.ij=u.i, x.pr=lapply(seq_len(dim(x.i)[1]), FUN=function(x)matrix(x.i[x, , ], nrow=dim(x.i)[2])),
                      SIMPLIFY=FALSE))

}

### Calculates the independence composite log-likelihood contribution of cluster i (clustering and stratification variable(s) differ)
composite_i_2 <- function(beta, lambda, l, u, tau, x, stratum.i, indiv.i){

  # Computing and aggregating the individual contributions to the cluster-specific composite likelihood
  do.call(sum, mapply(FUN=function(stratum.ij, indiv.ij){
    tau.s <- tau[[stratum.ij]]; lambda.s <- lambda[[stratum.ij]]
    l.ij <- l[[stratum.ij]][indiv.ij]; u.ij <- u[[stratum.ij]][indiv.ij]
    x.pr <- adrop(x[[stratum.ij]][indiv.ij, , , drop=FALSE], drop=1)
    composite_ij(l.ij, u.ij, tau.s, x.pr, lambda.s, beta)
  }, stratum.ij=stratum.i, indiv.ij=indiv.i, SIMPLIFY=FALSE))

}
### Calculates the independence composite log-likelihood contribution of individual j in cluster i
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
