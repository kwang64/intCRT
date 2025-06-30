### Performs a single iteration of the composite EM algorithm
point_iter <- function(beta, lambda, l, u, u.star, tau, x, risk, weights=NULL){

  # Expectation step
  w <- calc_w(l, u, tau, x, lambda, beta)

  # Maximization step
  beta.iter <- update_beta(beta, u.star, tau, x, w, risk, weights)
  lambda.iter <- update_lambda(beta.iter, x, w, risk, weights)

  return(list(beta=beta.iter, lambda=lambda.iter))
}
