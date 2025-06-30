### Returns the estimated baseline hazard function(s) that maximize(s) the profile composite log-likelihood
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
