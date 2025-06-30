### Solves the profile composite likelihood score equation for beta
update_beta <- function(beta, u.star, tau, x, w, risk, weights=NULL){

  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }

  # Finding and returning zeroes of the profile composite score equation for beta
  suppressWarnings(multiroot(beta_score, start=beta, u.star=u.star, tau=tau, x=x, w=w, risk=risk, weights=weights, maxiter=1)$root, classes="warning")

}
### Evaluates the profile composite likelihood score equation for beta
beta_score <- function(beta, u.star, tau, x, w, risk, weights){

  ### Evaluates the profile composite likelihood score equation in stratum s
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

### Computes the mean covariate vector at each timepoint in stratum s
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

### Updates all baseline hazard components
update_lambda <- function(beta, x, w, risk, weights=NULL){

  if (is.null(weights)){
    weights <- lapply(x, FUN=function(x.s) rep(1, dim(x.s)[1]))
  }

  # Updating the estimated baseline hazard function for each stratum and aggregating as a list
  mapply(function(x.s, w.s, risk.vr, weights.s){
    update_lambda_s(beta, x.s, w.s, risk.vr, weights.s)
  }, x.s=x, w.s=w, risk.vr=risk, weights.s=weights, SIMPLIFY=FALSE)

}

### Updates the baseline hazard components for stratum s
update_lambda_s <- function(beta, x.s, w.s, risk.vr, weights.s){

  # Calculating the denominator at each timepoint
  constant.vr <- do.call(rbind, lapply(seq_len(dim(x.s)[1]), FUN=function(x) exp(beta %*% x.s[x, , ])))
  denom <- colSums(risk.vr * constant.vr * weights.s)

  # Calculating the numerator at each timepoint
  num <- colSums(risk.vr * w.s * weights.s)

  return(num/denom)

}
