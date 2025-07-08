#' Compute Individual Survival Curve under UTT Policy Intervention
#'
#' Returns the survival curve for an individual in a given stratum, adjusting for the effect of a
#' Universal Test-and-Treat (UTT) policy and intervention assignment.
#'
#' @param baseline.i List. A list of baseline survival components for stratum \code{i}, including:
#' \itemize{
#'   \item \code{survival}: a data frame with columns \code{time} and \code{surv}.
#'   \item \code{density}: (not used in this function, but typically part of the baseline).
#'   \item \code{hazard}: (not used here, but included in the baseline list).
#' }
#' @param tau.i Numeric. Internal time point indicating the implementation of the UTT policy.
#' @param beta Numeric vector. Coefficient vector with elements corresponding to:
#' \itemize{
#'   \item \code{beta[1]}: pre-UTT intervention effect
#'   \item \code{beta[2]}: main effect of UTT
#'   \item \code{beta[3]}: interaction effect of UTT with intervention assignment
#' }
#' @param x.i Numeric. Indicator variable for individual intervention assignment (e.g., 0 = control, 1 = treated).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{time}{The time points from the baseline survival curve.}
#'   \item{surv}{The adjusted individual survival probabilities at each time point.}
#' }
#'
#' @details The survival function is modified based on the individual's intervention assignment and the
#' timing of the UTT policy. Prior to \code{tau.i}, survival probabilities are scaled by \code{exp(beta[1] * x.i)}.
#' Post-\code{tau.i}, an interaction model accounts for both the UTT main and interaction effects.
#'
#' @examples
#' # Simulated baseline curve
#' base <- return_exp_baseline(maxT = 10, lambda0 = 0.1, alpha = 0)
#' get_surv(base, tau.i = 5, beta = c(0.3, -0.2, 0.1), x.i = 1)
#'
#' @export

get_surv <- function(baseline.i, tau.i, beta, x.i){

  surv.i <- baseline.i$survival
  index <- sum(surv.i$time <= tau.i)

  # Calculating survival at times prior to the policy implementation
  s.prior <- surv.i$surv[1:index]^(exp(beta[1]*x.i))

  if (index == nrow(surv.i)){
    s.curve <- data.frame(time = surv.i$time, surv=s.prior)
  } else{
    # Calculating survival at times post the policy implementation
    s.tau.i <- surv.i$surv[index]
    s.post <- s.tau.i^(exp(beta[1]*x.i))*s.tau.i^(-exp(beta[1]*x.i + beta[2] + beta[3]*x.i))*surv.i$surv[(index+1):nrow(surv.i)]^(exp(beta[1]*x.i + beta[2] + beta[3]*x.i))

    s.curve <- data.frame(time=surv.i$time, surv=c(s.prior, s.post))
  }

  # Replacing missing/NaN survival times with zero
  if (sum(is.na(s.curve$surv)) > 0){
    na.index <- which(is.na(s.curve$surv))
    s.curve$surv[na.index] <- 0
  }

  return(s.curve)
}

#' Compute Failure Time Corresponding to a Given Survival Probability
#'
#' Returns the failure time (i.e., time-to-event) at which an individual's survival curve
#' reaches a specified survival probability threshold.
#'
#' @param s Numeric. Survival probability threshold (a value between 0 and 1).
#' @param tau.i Numeric. Time of Universal Test-and-Treat (UTT) policy implementation.
#' @param beta Numeric vector. Coefficient vector with elements corresponding to:
#' \itemize{
#'   \item \code{beta[1]}: pre-UTT intervention effect
#'   \item \code{beta[2]}: main effect of UTT
#'   \item \code{beta[3]}: interaction effect of UTT with intervention assignment
#' }
#' @param x.i Numeric. Indicator for individual intervention assignment (e.g., 0 = control, 1 = treated).
#' @param baseline.i List. Baseline survival data for stratum \code{i}, must include:
#' \itemize{
#'   \item \code{survival}: a data frame with columns \code{time} and \code{surv}.
#'   \item \code{density}, \code{hazard}: included for compatibility but not used directly.
#' }
#'
#' @return A numeric value representing the earliest failure time at which the individual's
#' survival probability drops below or equals the specified value \code{s}.
#'
#' @details The function first computes the individual's survival curve using \code{\link{get_surv}},
#' then finds the first time point where the survival probability falls below or equals \code{s}.
#' If \code{s} is greater than all survival probabilities, \code{NA} is returned.
#'
#' @examples
#' base <- return_exp_baseline(maxT = 10, lambda0 = 0.1, alpha = 0)
#' get_inv_surv(s = 0.5, tau.i = 5, beta = c(0.3, -0.2, 0.1), x.i = 1, baseline.i = base)
#'
#' @export

get_inv_surv <- function(s, tau.i, beta, x.i, baseline.i){

  # Obtaining piecewise-approximated individual survival function
  s.curve.i <- get_surv(baseline.i, tau.i, beta, x.i)

  # Returning the failure time corresponding to s
  index <- which(s.curve.i$surv <= s)[1]
  t <- s.curve.i$time[index]
  return(t)
}


#' Generate Interval-Censored Long-Form Data for an Individual with Time-Varying UTT Indicator
#'
#' Constructs the long-form data for a single individual in a survival study with time-varying
#' intervention (UTT policy) and potential interval censoring due to scheduled study visits and loss to follow-up.
#'
#' @param u.t Numeric. Uniform random variate used to simulate failure time.
#' @param tau.i Numeric. Time of the UTT policy implementation for the individual's stratum.
#' @param beta Numeric vector. Coefficient vector for the survival model (including pre-UTT effect, UTT effect, and interaction).
#' @param x.i Numeric. Indicator for intervention assignment (e.g., 0 = control, 1 = treated).
#' @param baseline.i List. Baseline survival components for the individual's stratum (from \code{return_exp_baseline} or similar), including:
#' \itemize{
#'   \item \code{survival}: Data frame of survival probabilities.
#'   \item \code{density}: Density estimates (not directly used).
#'   \item \code{hazard}: Hazard rates (not directly used).
#' }
#' @param visits.i Numeric vector. Scheduled HIV testing visit times for the individual.
#' @param cens Numeric. Censoring rate (0 if no censoring).
#' @param stratum Numeric. Stratum ID to which the individual belongs.
#' @param group Numeric. Community or group ID.
#' @param id Numeric. Unique identifier for the individual.
#'
#' @return A data frame with two rows (pre- and post-UTT policy) and the following columns:
#' \describe{
#'   \item{stratum, group, id}{Identifiers for stratum, group, and individual.}
#'   \item{treat}{Indicator for intervention assignment.}
#'   \item{utt}{Time-varying UTT policy indicator (0 = before policy, 1 = after).}
#'   \item{start, stop}{Start and stop times of each time interval.}
#'   \item{left, right}{Interval-censored survival information for the individual.}
#' }
#'
#' @details The function simulates a failure time from an inverse survival function, applies right censoring if specified, and
#' determines interval-censored observation based on scheduled visits. It returns two rows for modeling time-varying effects of UTT.
#'
#' @examples
#' base <- return_exp_baseline(10, lambda0 = 0.1, alpha = 0)
#' visits <- seq(1, 9, by = 2)
#' gen_individual_tvc(u.t = runif(1), tau.i = 5, beta = c(0.2, -0.1, 0.05),
#'                    x.i = 1, baseline.i = base, visits.i = visits,
#'                    cens = 0.01, stratum = 1, group = 1, id = 1001)
#'
#' @export

gen_individual_tvc <- function(u.t, tau.i, beta, x.i, baseline.i, visits.i, cens, stratum, group, id){

  # Generating survival and censoring times
  s.time <- get_inv_surv(u.t, tau.i, beta, x.i, baseline.i)
  if (cens==0){
    c.time <- visits.i[length(visits.i)] + 1
  } else{
    u.c <- runif(1)
    c.time <- -log(u.c)/cens
  }

  # Deriving final interval censored observation and return (left, right)
  visits.i[visits.i > c.time] <- NA
  event <- as.numeric(visits.i > s.time)
  if (sum(is.na(event))==length(event)){
    out <- data.frame(stratum=stratum, group=group, treat=x.i, left=0, right=Inf)
    return(out)
  }
  event.status <- max(event, na.rm=TRUE)
  if (event.status==1){
    first.event.time <- min(which(event==1))
    l <- ifelse(first.event.time==1, 0,
                ifelse(sum(!is.na(visits.i[1:(first.event.time-1)]))==0, 0,
                       max(visits.i[1:(first.event.time-1)], na.rm=TRUE)))
    r <- visits.i[first.event.time]
  } else{
    l <- max(visits.i, na.rm=TRUE)
    r <- Inf
  }

  # Returning the final long-form data for individual i
  out <- data.frame(stratum=rep(stratum, 2), group=rep(group, 2), id=rep(id, 2), treat=rep(x.i, 2), utt=c(0, 1),
                    start=c(0, tau.i), stop=c(tau.i, Inf), left=rep(l, 2), right=rep(r, 2))
  return(out)
}

#' Generate Stratified Dataset with Time-Varying Covariates and Interval-Censored Survival Times
#'
#' Constructs a stratified, clustered dataset for a survival study with a binary treatment assignment
#' (constant within each cluster) and a time-varying intervention (UTT policy). The outcome is subject
#' to interval censoring and optional within-cluster dependence.
#'
#' @param M Numeric. Number of clusters or communities.
#' @param ni Numeric vector of length 2. Range (min, max) of individuals per community, sampled uniformly.
#' @param beta Numeric vector. Coefficients for: (1) the time-constant treatment variable, (2) the time-varying UTT indicator, and (3) their interaction.
#' @param gen_visits Function. Function to generate internal visit times (e.g., HIV testing dates).
#' @param cens Numeric. Loss-to-follow-up (censoring) rate. 0 implies no censoring.
#' @param rho Numeric. Kendall’s tau indicating within-cluster dependence. Set to 0 for independence.
#' @param baseline List. List of baseline distribution functions (generated from `return_baseline()` or `return_exp_baseline()`) for each stratum.
#' @param S Either a function that returns stratum assignment (e.g., sampling function), or a vector of length `M` giving fixed stratum assignments for each cluster.
#'
#' @return A data frame in long format with two rows per individual (before and after UTT policy) and the following columns:
#' \describe{
#'   \item{id}{Individual ID.}
#'   \item{stratum, group}{Stratum and cluster identifiers.}
#'   \item{treat}{Cluster-assigned binary treatment indicator.}
#'   \item{utt}{Time-varying UTT indicator (0 before policy, 1 after).}
#'   \item{start, stop}{Start and stop times of the interval.}
#'   \item{left, right}{Interval-censored survival times.}
#' }
#'
#' @details This function simulates data for survival models with interval censoring and time-varying covariates.
#' The time of UTT implementation varies across individuals. Dependence within clusters can be introduced via a
#' Clayton copula if `rho > 0`.
#'
#' @examples
#' set.seed(1)
#' baseline <- list(return_exp_baseline(10, 0.2, 0))
#' gen_visits <- function(max.time=10) sort(runif(5, 0, max.time))
#' simdat <- gen_time_varying_cov_tx(M = 4, ni = c(4, 6),
#'                                   beta = c(0.1, -0.2, 0.05),
#'                                   gen_visits = gen_visits,
#'                                   cens = 0.01,
#'                                   rho = 0,
#'                                   baseline = baseline,
#'                                   S = rep(1, 4))
#'
#' @export

gen_time_varying_cov_tx <- function(M, ni, beta, gen_visits, cens, rho=0, baseline, S){

  # Generating the design matrix if observations are independent within each cluster
  if (rho==0){

    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }

    # Generating the design matrix if observations belong to an exchangeable subnetwork within each cluster
  } else {

    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }

  }

  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }

  data <- do.call(rbind, data.list.2)
  n <- nrow(data)

  # Generating (internal) inspection times
  visits <- lapply(data[,4], FUN=gen_visits)

  # Generating (internal) time of UTT implementation
  tau <- runif(n, 0, max(do.call(c, visits)))

  data <- cbind(seq(1, n, by=1), data, tau)
  data.list.3 <- split(data, seq(nrow(data)))

  # Generating long-form dataset
  out <- do.call(rbind, mapply(FUN=function(y, v){
    stratum <- y[2]; group <- y[3]; id <- y[1]
    u.t <- y[4]; x.i <- y[5]; tau.i <- y[6]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvc(u.t, tau.i, beta, x.i, baseline.i, v, cens, stratum, group, id)
  },
  y=data.list.3, v=visits, SIMPLIFY=FALSE))

  # Returning the final dataset
  return(out)
}

#' Generate Stratified Dataset with Clustered Binary Covariate and Time-Varying UTT Indicator
#'
#' Simulates a stratified dataset for interval-censored survival analysis where each cluster (community)
#' is assigned a binary treatment indicator (constant within cluster), and individuals are affected by a
#' time-varying intervention (UTT). Within-cluster dependence can be introduced via a copula.
#'
#' @param M Numeric. Number of communities (clusters).
#' @param ni Numeric vector of length 2. The minimum and maximum number of individuals per community. Sampled uniformly.
#' @param beta Numeric vector. Coefficients for: (1) the cluster-level treatment variable, (2) the time-varying UTT indicator, and (3) their interaction.
#' @param gen_visits Function. A function to generate internal visit times (e.g., HIV testing dates) for each individual.
#' @param cens Numeric. Loss-to-follow-up rate. Set to 0 for no censoring.
#' @param rho Numeric. Kendall’s tau indicating the strength of within-cluster dependence. Set to 0 for independence.
#' @param baseline List. A list of baseline distribution functions (e.g., from \code{return_baseline()} or \code{return_exp_baseline()}) for each stratum.
#' @param S Either a function that generates a vector of stratum memberships (of length equal to cluster size) or a vector of length \code{M} indicating fixed stratum assignments for each cluster.
#'
#' @return A long-format \code{data.frame} with two rows per individual (pre- and post-UTT), including:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{stratum, group}{Stratum and cluster IDs.}
#'   \item{treat}{Binary treatment indicator (constant within cluster).}
#'   \item{utt}{Time-varying UTT policy indicator (0 before, 1 after).}
#'   \item{start, stop}{Time interval corresponding to each policy stage.}
#'   \item{left, right}{Interval-censored survival time for each individual.}
#' }
#'
#' @details If \code{rho > 0}, the function introduces dependence among individuals in a cluster using a Clayton copula. Otherwise, individual failure times are generated independently. Each individual is subject to internal visit times and possible loss to follow-up.
#'
#' @examples
#' set.seed(123)
#' baseline <- list(return_exp_baseline(10, lambda0 = 0.2, alpha = 0))
#' gen_visits <- function(n) t(replicate(n, sort(runif(4, 0, 10))))
#' simdat <- gen_time_varying_cov(M = 4, ni = c(5, 6),
#'                                 beta = c(0.1, -0.2, 0.05),
#'                                 gen_visits = gen_visits,
#'                                 cens = 0.01,
#'                                 rho = 0,
#'                                 baseline = baseline,
#'                                 S = rep(1, 4))
#'
#' @export

gen_time_varying_cov <- function(M, ni, beta, gen_visits, cens, rho=0, baseline, S){

  # Generating the design matrix if observations are independent within each cluster
  if (rho==0){

    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      n.i <- sample(ni[1]:ni[2], 1)
      u <- runif(n.i, 0, 1)
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }

    # Generating the design matrix if observations belong to an exchangeable subnetwork within each cluster
  } else {

    data.list <- list()
    tx.sequence <- rep(c(1, 0), M/2)
    for (i in 1:M){
      networks <- rpois(ni[2]/2, 2) + 2
      networks.cm <- cumsum(networks)
      n.i <- sample(ni[1]:ni[2], 1)
      sub.ni <- networks[networks.cm < n.i]
      copula.fn.list <- lapply(sub.ni, FUN=function(x) claytonCopula(-2*rho/(rho-1), dim = x))
      u <- do.call(c, lapply(copula.fn.list, FUN=function(x) as.numeric(rCopula(1, x))))
      temp <- matrix(NA, nrow=length(u), ncol=3)
      temp[, 1] <- i; temp[, 2] <- u; temp[, 3] <- tx.sequence[i]
      data.list[[length(data.list) + 1]] <- temp
    }

  }

  # Determining stratum membership
  if (is.function(S)) {
    data.list.2 <- lapply(data.list, FUN=function(dat){
      n.i <- nrow(dat)
      s <- S(n.i)
      cbind(s, dat)
    })
  } else{
    data.list.2 <- mapply(FUN=function(s, dat){
      cbind(rep(s, nrow(dat)), dat)
    }, s=S, dat=data.list, SIMPLIFY=FALSE)
  }

  data <- do.call(rbind, data.list.2)
  n <- nrow(data)

  # Generating (internal) inspection times
  visits <- gen_visits(n)

  # Generating (internal) time of UTT implementation
  tau <- runif(n, 0, max(visits[, ncol(visits)]))

  data <- cbind(seq(1, n, by=1), data, tau, visits)

  # Generating long-form dataset
  out <- do.call(rbind, apply(data, 1, FUN=function(y){
    stratum <- y[2]; group <- y[3]; id <- y[1]
    u.t <- y[4]; x.i <- y[5]; tau.i <- y[6]; visits.i <- y[7:length(y)]
    baseline.i <- baseline[[stratum]]
    gen_individual_tvc(u.t, tau.i, beta, x.i, baseline.i, visits.i, cens, stratum, group, id)
  }))

  # Returning the final dataset
  return(out)
}
