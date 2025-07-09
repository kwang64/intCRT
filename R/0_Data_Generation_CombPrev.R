#' Time-Dependent Covariate Effect Functions
#'
#' Provides a list of possible functional forms for modeling the time-varying effect
#' of a covariate in a Cox model or other survival models.
#'
#' Each element in the list corresponds to a specific form—logarithmic, linear, or piecewise—
#' and contains:
#' \describe{
#'   \item{name}{A string identifying the type of function.}
#'   \item{args}{A character vector naming the arguments used by the function.}
#'   \item{beta.t}{A function specifying how the time-varying coefficient \eqn{\beta(t)} is computed at time \eqn{u}, given coefficient parameters \eqn{\beta} (and \eqn{\tau} if applicable).}
#' }
#'
#' ## Functional Forms:
#' \describe{
#'   \item{log}{\eqn{\beta(t) = \beta_1 + \beta_2 \log(t)}}
#'   \item{linear}{\eqn{\beta(t) = \beta_1 + \beta_2 t}}
#'   \item{piecewise}{\eqn{\beta(t) = \beta_1 + \beta_2 \cdot I(t > \tau)} where \eqn{I} is the indicator function.}
#' }
#'
#' @format A named list of three time-varying coefficient function specifications.
#' @examples
#' # Evaluate beta(t) for t = 5 using the log form
#' tvbeta_fns$log$beta.t(u = 5, beta = c(1, 0.5))
#'
#' # Use the piecewise form with a threshold tau = 3
#' tvbeta_fns$piecewise$beta.t(u = 5, beta = c(1, 0.5), tau = 3)
#'
#' @export

tvbeta_fns <- list(
  log = list(
    name = "log",
    args = c("u", "beta"),
    beta.t = function(u, beta) beta[1] + beta[2]*log(u)
  ),
  linear = list(
    name = "linear",
    args = c("u", "beta"),
    beta.t = function(u, beta) beta[1] + beta[2]*u
  ),
  piecewise = list(
    name = "piecewise",
    args = c("u", "beta", "tau"),
    beta.t = function(u, beta, tau) beta[1] + ifelse(u > tau, 1, 0)*beta[2]
  )
)

#' Compute Individual Survival Curve with Time-Dependent Covariate Effects
#'
#' Computes the survival curve \eqn{S_{ij}(t)} for an individual given a stratum-specific baseline hazard
#' and a time-dependent covariate effect \eqn{\beta(t)}.
#'
#' The function evaluates the linear predictor over time using a specified time-varying effect form,
#' adjusts the hazard accordingly, and approximates the survival function by numerically integrating
#' the cumulative hazard.
#'
#' @param baseline.i List. Contains stratum-specific baseline functions:
#' \code{survival}, \code{density}, and \code{hazard}, each as data frames with time.
#' @param bfns List. A list defining the functional form for the time-varying effect \eqn{\beta(t)}.
#' Must contain a function \code{beta.t()}.
#' @param aux List. Additional arguments required by \code{beta.t()}, such as parameters \code{beta} and possibly \code{tau}.
#' @param x.i Numeric. Indicator of intervention assignment (e.g., 0 or 1).
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{time}{Time points at which the survival curve is evaluated.}
#'   \item{surv}{Estimated survival probability at each time point.}
#' }
#'
#' @examples
#' # Example assuming you already have baseline.i, bfns, aux, and x.i defined:
#' # result <- get_surv_tvb(baseline.i, bfns, aux, x.i = 1)
#' # plot(result$time, result$surv, type = "l")
#'
#' @export

get_surv_tvb <- function(baseline.i, bfns, aux, x.i){

  h.i <- baseline.i$hazard

  # Obtaining the linear predictor at each timepoint
  aux$u <- h.i$time
  lp <- do.call(bfns$beta.t, aux)*x.i
  if (sum(is.na(lp)) > 0){
    na.index <- which(is.na(lp))
    if (length(na.index)==1){
      if (na.index != length(lp)){
        lp[na.index] <- lp[na.index + 1]
      } else if (na.index == length(lp)){
        lp[na.index] <- lp[na.index -1]
      }
    } else{
      for (i in na.index){
        if (i==length(lp)){
          lp[i] <- lp[tail(seq(1, i-1, by=1)[which(!is.na(lp[1:(i-1)]))], 1)]
        } else if (length(seq(i+1, length(lp), by=1)[which(!is.na(lp[(i+1):length(lp)]))]) > 0){
          lp[i] <- lp[head(seq(i+1, length(lp), by=1)[which(!is.na(lp[(i+1):length(lp)]))], 1)]
        } else{
          lp[i] <- lp[tail(seq(1, i-1, by=1)[which(!is.na(lp[1:(i-1)]))], 1)]
        }
      }
    }
  }
  if (sum(lp==Inf) > 0){
    inf.index <- which(lp==Inf)
    if ((length(inf.index)==1) & (inf.index != length(lp))){
      lp[inf.index] <- lp[inf.index + 1]
    } else if ((length(inf.index)==1) & (inf.index == length(lp))){
      lp[inf.index] <- lp[inf.index -1]
    } else {
      for (i in inf.index){
        lp[i] <- ifelse(i==length(lp), lp[tail(seq(1, i-1, by=1)[which(!(lp[1:(i-1)]==Inf))], 1)],
                        lp[head(seq(i+1, length(lp), by=1)[which(!(lp[(i+1):length(lp)]==Inf))], 1)])
      }
    }
  }
  h.ij <-  h.i$hazard*exp(lp)

  # Approximating the cumulative hazard
  H.ij <- cumsum(h.ij[-length(h.ij)]*diff(h.i$time))

  # Approximating the corresponding survival
  S.ij <- exp(-H.ij)
  # Replacing missing/NA/infinite survival times with zero
  if (sum(is.na(S.ij)) > 0){
    na.index <- which(is.na(S.ij))
    S.ij[na.index] <- 0
  }
  if (sum(S.ij==Inf) > 0){
    inf.index <- which(S.ij==Inf)
    S.ij[inf.index] <- 0
  }

  s.curve <- data.frame(time=h.i$time, surv=c(S.ij, 0))
  return(s.curve)
}

#' Invert Individual Survival Curve to Obtain Failure Time
#'
#' Returns the failure time \eqn{t} corresponding to a specified survival probability \eqn{S(t) = s}
#' for an individual, based on a time-dependent covariate effect and stratum-specific baseline hazard.
#'
#' Internally, this function computes the individual survival curve using \code{\link{get_surv_tvb}}, then finds
#' the earliest time at which the survival curve falls below or equals the target survival probability.
#'
#' @param s Numeric. Target survival probability (e.g., 0.5 for median survival).
#' @param baseline.i List. Contains stratum-specific baseline functions: \code{survival}, \code{density}, and \code{hazard}, each as data frames with time.
#' @param bfns List. A list defining the functional form for the time-varying effect \eqn{\beta(t)}.
#' Must contain a function \code{beta.t()}.
#' @param aux List. Additional arguments required by \code{beta.t()}, such as parameters \code{beta} and possibly \code{tau}.
#' @param x.i Numeric. Indicator of intervention assignment (e.g., 0 or 1).
#'
#' @return Numeric. The estimated failure time \eqn{t} such that \eqn{S(t) = s}.
#'
#' @examples
#' # Example assuming you already have baseline.i, bfns, aux, and x.i defined:
#' # median_time <- get_inv_surv_tvb(s = 0.5, baseline.i, bfns, aux, x.i = 1)
#'
#' @seealso \code{\link{get_surv_tvb}}
#'
#' @export

get_inv_surv_tvb <- function(s, baseline.i, bfns, aux, x.i){

  # Obtaining piecewise-approximated individual survival function
  s.curve.i <- get_surv_tvb(baseline.i, bfns, aux, x.i)

  # Returning the failure time corresponding to s
  index <- which(s.curve.i$surv <= s)[1]
  t <- s.curve.i$time[index]
  return(t)
}

#' Generate Interval-Censored Individual Data with Time-Dependent Effects
#'
#' Simulates an individual's interval-censored time-to-event data under a time-dependent covariate effect
#' using a specified baseline hazard and a given intervention assignment.
#'
#' The function draws a survival time based on a uniform random variate and inverts the individual survival curve.
#' It then draws a censoring time (if applicable), and determines the interval \code{[left, right)} within which the event occurred,
#' based on the timing of scheduled study visits.
#'
#' @param u.t Numeric. Uniform random variate used to simulate failure time via inverse survival function.
#' @param baseline.i List. Contains stratum-specific baseline functions: \code{survival}, \code{density}, and \code{hazard}, each as data frames with time.
#' @param bfns List. List containing a function \code{beta.t} defining the time-dependent covariate effect.
#' @param aux List. Additional arguments required by \code{beta.t}, such as \code{beta} and optionally \code{tau}.
#' @param x.i Numeric. Indicator of intervention assignment (e.g., 0 or 1).
#' @param visits.i Numeric vector. Times of all study-mandated follow-up visits for HIV testing.
#' @param cens Numeric. Loss-to-follow-up rate; 0 indicates no censoring.
#' @param stratum Numeric. ID number of the individual's associated stratum.
#' @param group Numeric. ID number of the individual's associated community or group.
#'
#' @return A one-row data frame with columns:
#' \describe{
#'   \item{stratum}{Stratum ID.}
#'   \item{group}{Group/community ID.}
#'   \item{treat}{Treatment assignment (0 or 1).}
#'   \item{left}{Left endpoint of interval-censored event time.}
#'   \item{right}{Right endpoint of interval-censored event time.}
#' }
#'
#' @examples
#' # Example assumes availability of: baseline.i, bfns, aux, visits.i
#' # gen_individual_tvbeta(u.t = runif(1), baseline.i, bfns, aux, x.i = 1,
#' #                       visits.i = c(1, 2, 3, 4, 5), cens = 0.1, stratum = 1, group = 2)
#'
#' @seealso \code{\link{get_inv_surv_tvb}}, \code{\link{get_surv_tvb}}
#'
#' @export

gen_individual_tvbeta <- function(u.t, baseline.i, bfns, aux, x.i, visits.i, cens, stratum, group){

  # Generating survival and censoring times
  s.time <- get_inv_surv_tvb(u.t, baseline.i, bfns, aux, x.i)
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

  # Returning the final observation for individual i
  out <- data.frame(stratum=stratum, group=group, treat=x.i, left=l, right=r)
  return(out)
}
#' Generate Stratified Time-to-Event Data with Time-Dependent Covariate Effects
#'
#' Simulates individual-level, interval-censored survival data across multiple communities and strata,
#' incorporating a time-dependent covariate effect (e.g., intervention assignment) and optional
#' within-community correlation through positive stable frailty or copula-based dependence.
#'
#' This function allows flexible specification of the baseline hazard and time-varying covariate effects,
#' enabling simulations of complex, realistic trial or cohort structures.
#'
#' @param M Integer. Number of communities (groups) to simulate.
#' @param ni Integer vector of length 2. Range of individual counts per community (e.g., \code{c(10, 20)}).
#' @param beta List or character. Either a list containing the function \code{beta.t()} for time-dependent effect,
#' or a string name matching one of the predefined options in \code{\link{tvbeta_fns}} (e.g., "log", "linear").
#' @param beta.args List. Arguments (e.g., \code{beta}, \code{tau}) required by the time-dependent effect function \code{beta.t()}.
#' @param gen_visits Function. A function that takes an individual's latent time and returns internal visit times.
#' @param cens Numeric. Loss-to-follow-up rate (0 means no censoring).
#' @param rho Numeric. Kendall's tau controlling dependence within sub-networks in a community. Set to 0 for independence.
#' @param alpha Numeric. Index parameter of the positive stable distribution for community-level frailty (0 means no frailty).
#' @param baseline List. A list of baseline survival/density/hazard distributions, one per stratum.
#' @param S Either a function or vector. If a function, it should generate stratum membership for each cluster.
#' If a vector, it assigns fixed stratum IDs to each cluster.
#'
#' @return A data frame with interval-censored survival data. Each row represents an individual and includes:
#' \describe{
#'   \item{stratum}{Stratum ID.}
#'   \item{group}{Community/group ID.}
#'   \item{treat}{Treatment assignment (0 or 1).}
#'   \item{left}{Left endpoint of the interval containing the failure time.}
#'   \item{right}{Right endpoint of the interval containing the failure time.}
#' }
#'
#' @examples
#' # Example setup (pseudo-code):
#' # baseline_list <- list(baseline1, baseline2)
#' # visits_fn <- function(u) seq(1, 10, by = 2)
#' # S_fn <- function(n) sample(1:2, n, replace = TRUE)
#' # df <- gen_time_dependent_beta_tx(M = 4, ni = c(10, 20),
#' #                                  beta = "log", beta.args = list(beta = c(0.5, -0.2)),
#' #                                  gen_visits = visits_fn, cens = 0.05,
#' #                                  rho = 0.2, alpha = 0.3, baseline = baseline_list, S = S_fn)
#'
#' @seealso \code{\link{gen_individual_tvbeta}}, \code{\link{tvbeta_fns}}, \code{\link{g}}
#'
#' @export

gen_time_dependent_beta_tx <- function(M, ni, beta, beta.args, gen_visits, cens, rho=0, alpha = 0, baseline, S){

  # Determining beta(t) parametrization
  if (missing(beta)){
    stop("Time-dependent covariate effect \"beta\" not specified")
  } else if (is.list(beta)){
    bfns <- beta
  } else {
    beta <- match.arg(tolower(beta), tolower(names(tvbeta_fns)))
    bfns <- tvbeta_fns[[beta]]
  }

  # Generating cluster-specific positive stable frailty terms
  if (alpha==0){
    omega <- rep(1, M)
  } else {
    omega <- rstable(M, alpha=alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
  }

  # Generating the design matrix if observations are independent within each stratum
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

    # Generating the design matrix if observations belong to an exchangeable subnetwork within each stratum
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
  data.list.3 <- split(data, seq(nrow(data)))

  # Generating dataset
  out <- do.call(rbind, mapply(FUN=function(y, v){
    stratum <- y[1]; group <- y[2]; u.t <- y[3]; x.i <- y[4]
    omega.i <- omega[group]
    baseline.i <- baseline[[stratum]]
    baseline.i$hazard$hazard <- baseline.i$hazard$hazard*omega.i
    gen_individual_tvbeta(u.t, baseline.i, bfns, beta.args, x.i, v, cens, stratum, group)
  },
  y=data.list.3, v=visits, SIMPLIFY=FALSE))

  # Returning the final dataset
  return(out)
}

#' Generate Stratified Survival Data with Time-Dependent Covariate Effects
#'
#' Simulates a stratified, community-level survival dataset with interval-censored failure times
#' and a time-dependent covariate effect (e.g., treatment or intervention). Individuals may be clustered
#' within communities and optionally connected through copula-based dependence or positive stable frailty.
#'
#' The time-dependent covariate effect may be specified as a named option (e.g., "log", "linear", "piecewise")
#' or supplied as a custom list with function \code{beta.t()}. When frailty (\code{alpha > 0}) is present,
#' the function applies a cluster-specific effect multiplier to the baseline hazard.
#'
#' @param M Integer. Number of communities (clusters).
#' @param ni Integer vector of length 2. Range of individual counts per community, e.g., \code{c(10, 20)}.
#' @param beta List or character. Either a list with a function \code{beta.t()} describing the time-dependent effect,
#' or a string matching one of \code{\link{tvbeta_fns}} options (e.g., "log", "linear").
#' If \code{alpha > 0}, this determines the cluster-conditional effect, and the marginal effect is \eqn{\beta(t) \cdot \alpha}.
#' @param beta.args List. Arguments required by the time-dependent effect function (e.g., \code{beta}, \code{tau}).
#' @param gen_visits Function. A function that takes the total number of individuals and returns a list of vectors
#' containing each individual's internal visit times.
#' @param cens Numeric. Loss-to-follow-up rate (0 indicates no censoring).
#' @param rho Numeric. Kendall's tau controlling intra-cluster dependence using Clayton copulas; set to 0 for independence.
#' @param alpha Numeric. Index parameter for the positive stable distribution; controls cluster-level frailty. Set to 0 for no frailty.
#' @param baseline List. A list of baseline hazard/survival/density functions (one per stratum).
#' @param S Function or vector. A function to generate stratum membership per cluster, or a fixed vector of stratum IDs (one per cluster).
#'
#' @return A data frame with interval-censored observations. Columns include:
#' \describe{
#'   \item{stratum}{Stratum ID.}
#'   \item{group}{Community/group ID.}
#'   \item{treat}{Treatment assignment (0 or 1).}
#'   \item{left}{Left endpoint of the event time interval.}
#'   \item{right}{Right endpoint of the event time interval.}
#' }
#'
#' @details
#' The function supports both frailty-based dependence using positive stable distributions (via \code{alpha})
#' and subnetwork dependence within clusters using Clayton copulas (via \code{rho}).
#' If both are zero, individuals are treated as conditionally independent.
#'
#' @examples
#' # Minimal example (assuming supporting functions and objects defined elsewhere):
#' # baseline_list <- list(baseline1, baseline2)
#' # gen_visits_fn <- function(n) replicate(n, sort(runif(5, 0, 10)), simplify = FALSE)
#' # S_fn <- function(n) sample(1:2, n, replace = TRUE)
#' # df <- gen_time_dependent_beta(M = 4, ni = c(10, 20),
#' #     beta = "linear", beta.args = list(beta = c(0.3, -0.1)),
#' #     gen_visits = gen_visits_fn, cens = 0.05,
#' #     rho = 0.2, alpha = 0.5,
#' #     baseline = baseline_list, S = S_fn)
#'
#' @seealso \code{\link{gen_individual_tvbeta}}, \code{\link{tvbeta_fns}}, \code{\link{get_surv_tvb}}
#'
#' @export

gen_time_dependent_beta <- function(M, ni, beta, beta.args, gen_visits, cens, rho=0, alpha=0, baseline, S){

  # Determining beta(t) parametrization
  if (missing(beta)){
    stop("Time-dependent covariate effect \"beta\" not specified")
  } else if (is.list(beta)){
    bfns <- beta
  } else {
    beta <- match.arg(tolower(beta), tolower(names(tvbeta_fns)))
    bfns <- tvbeta_fns[[beta]]
  }

  # Generating cluster-specific positive stable frailty terms
  if (alpha==0){
    omega <- rep(1, M)
  } else {
    omega <- rstable(M, alpha=alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
  }

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
  data <- cbind(data, visits)

  # Generating dataset
  out <- do.call(rbind, apply(data, 1, FUN=function(y){
    stratum <- y[1]; group <- y[2]; u.t <- y[3]; x.i <- y[4]
    visits.i <- y[5:length(y)]
    omega.i <- omega[group]
    baseline.i <- baseline[[stratum]]
    baseline.i$hazard$hazard <- baseline.i$hazard$hazard*omega.i
    gen_individual_tvbeta(u.t, baseline.i, bfns, beta.args, x.i, visits.i, cens, stratum, group)
  }))

  # Returning the final dataset
  return(out)
}
