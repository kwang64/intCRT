### Calculates the conditional expectation of the Poisson augmentation variables
calc_w <- function(l, u, tau, x, lambda, beta){

  # Taking the projection of w.s for each stratum
  mapply(function(l.s, u.s, tau.s, x.s, lambda.s){
    w_s(l.s, u.s, tau.s, x.s, lambda.s, beta)
  }, l.s=l, u.s=u, tau.s=tau, x.s=x, lambda.s=lambda, SIMPLIFY=FALSE)

}

### Returns the matrix of projected w.s for stratum s
w_s <- function(l.s, u.s, tau.s, x.s, lambda.s, beta){

  # Taking the projection of w.sv for each individual
  do.call(rbind, mapply(FUN=function(l.sv, u.sv, x.pr) w_sv(l.sv, u.sv, tau.s, x.pr, lambda.s, beta),
                        l.sv=l.s, u.sv=u.s, x.pr=lapply(seq_len(dim(x.s)[1]), FUN=function(x)matrix(x.s[x, , ], nrow=dim(x.s)[2])), SIMPLIFY=FALSE))

}

### Returns the vector w.sv for individual v in stratum s
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
