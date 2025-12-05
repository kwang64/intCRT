#' Composite Maximum Likelihood Estimation for Marginal Cox Models with Clustered Interval-Censored Data
#'
#' Obtains composite maximum likelihood estimators for marginal Cox regression models
#' with clustered interval-censored data. Supports optional stratification, clustering,
#' time-varying covariates, and multiple methods for variance estimation.
#'
#' @importFrom abind adrop

#' @param formula A formula object specifying the marginal Cox model (e.g., `Surv(start, stop, event) ~ covariates`). Required if \code{data} is unprocessed.
#' @param data A \code{data.frame} or a pre-processed \code{list} with two elements:
#'   \describe{
#'     \item{\code{pointEstimation}}{A list containing: \code{l}, \code{u}, \code{u.star}, \code{tau}, \code{risk}, \code{x}, each a list of length S (number of strata).}
#'     \item{\code{varMapping}}{Optional list of cluster-to-stratum mapping. Contains:
#'       \code{stratum} (list of stratum IDs per cluster) and
#'       \code{indiv} (list of individual IDs per cluster). May be \code{NULL} if not needed.}
#'   }
#' @param strata Character vector specifying stratification variable(s) in \code{data}.
#' @param clus Character vector specifying clustering variable(s) in \code{data}.
#' @param id Character name of the subject identifier variable, used if covariates are time-varying.
#' @param variance Logical; whether to compute and return the variance-covariance matrix of regression coefficients.
#' @param method Character string; either \code{"profile"} (default) or \code{"bootstrap"}, specifying the method for variance estimation.
#' @param control A named list of EM algorithm control parameters:
#'   \describe{
#'     \item{\code{init}}{Initial values for coefficients and baseline hazard(s).}
#'     \item{\code{tol}}{Convergence tolerance for the L1 norm of coefficient updates.}
#'     \item{\code{maxit}}{Maximum number of EM iterations.}
#'     \item{\code{c}}{Perturbation constant for numerical derivatives.}
#'     \item{\code{trace}}{Logical; if \code{TRUE}, stores the full sequence of EM iterates.}
#'   }
#' @param boot.control A named list of bootstrap control parameters:
#'   \describe{
#'     \item{\code{perturb}}{Logical; whether to use perturbation bootstrap (\code{TRUE}) or nonparametric bootstrap (\code{FALSE}).}
#'     \item{\code{resample.clusters}}{Logical; if \code{TRUE}, resample clusters instead of individuals. Used only if \code{clus} is specified.}
#'     \item{\code{n.boot}}{Number of bootstrap replicates.}
#'   }
#' @param boot.data Optional list of pre-processed bootstrap datasets (only used if \code{data} is already processed and \code{perturb = FALSE}).
#'
#' @return A list of class \code{"compCoxIC"} with the following components:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{data.attr}}{List of data attributes: number of observations, strata, events.}
#'   \item{\code{coefficients}}{Estimated regression coefficients.}
#'   \item{\code{var}}{Variance-covariance matrix of coefficients, if requested.}
#'   \item{\code{var.method}}{Method used for variance estimation ("profile" or "bootstrap").}
#'   \item{\code{var.args}}{Arguments used for variance estimation.}
#'   \item{\code{baseline.hazard}}{List of estimated baseline hazards, one per stratum.}
#'   \item{\code{n.iter}}{Number of EM iterations completed.}
#'   \item{\code{l1.norm}}{Final L1 norm of coefficient updates (used for convergence check).}
#'   \item{\code{trace}}{Trace of EM iterations if \code{trace = TRUE}.}
#'   \item{\code{fit.times}}{Time used for model fitting and variance estimation.}
#' }
#'
#' @export

composite_coxIC <- function(formula=NULL, data, strata=NULL, clus=NULL, id=NULL, variance=TRUE, method="profile", control=NULL, boot.control=NULL, boot.data=NULL){

  call <- match.call()

  args <- list(init=NULL, tol=1e-4, maxit.fit=1000, maxit.var=1000, c=1, trace=FALSE)
  if (!is.null(control)){
    if ("maxit" %in% names(control)){
      args$maxit.fit <- control$maxit
      args$maxit.var <- control$maxit
      args[setdiff(names(control), "maxit")] <- control[-which(names(control)=="maxit")]
    } else {
      args[names(control)] <- control
    }
  }

  boot.args <- list(perturb=FALSE, resample.clusters=FALSE, n.boot=100)
  if (!is.null(boot.control)){
    boot.args[names(boot.control)] <- boot.control
  }
  if (!is.data.frame(data) & method=="bootstrap" & is.null(boot.data) & boot.args$perturb==FALSE){
    stop("To conduct bootstrap inference for analyses with pre-processed data, provide a list of pre-processed bootstrapped datasets using the boot.data option.")
  }

  # Processing the input data
  if (is.data.frame(data)){
    temp <- process_data(formula, data, strata, id, clus)
    processed <- temp$pointEstimation
    var.mapping <- temp$varMapping
  } else{
    processed <- data$pointEstimation
    var.mapping <- data$varMapping
  }
  l <- processed$l; u <- processed$u; u.star <- processed$u.star; tau <- processed$tau
  risk <- processed$risk
  x <- processed$x

  # Initializing the point estimates
  if (is.null(args$init)){
    init <- em_init(dim(x[[1]])[2], tau)
  } else {
    init <- args$init
  }
  beta.current <- init$beta; lambda.current <- init$lambda

  # Initializing lists to store point estimates at each iteration
  if (args$trace){
    trace.beta <- list(beta.current)
    trace.lambda <- list(lambda.current)
  }

  # Iterating until convergence
  fit.start <- Sys.time()
  l1.norm <- 1; n.iter <- 0
  while (l1.norm > args$tol & n.iter < args$maxit.fit){
    update <- point_iter(beta.current, lambda.current, l, u, u.star, tau, x, risk)
    l1.norm <- sum(abs(update$beta - beta.current))
    beta.current <- update$beta; lambda.current <- update$lambda
    if (args$trace){
      trace.beta[[length(trace.beta)+1]] <- beta.current
      trace.lambda[[length(trace.lambda)+1]] <- lambda.current
    }
    n.iter <- n.iter + 1
  }
  fit.end <- Sys.time()

  # Estimating the covariance matrix for beta
  var.start <- Sys.time()
  if (variance){
    if (method=="profile"){
      beta.cov <- variance_beta(beta.current, processed, var.mapping, args)
      colnames(beta.cov) <- dimnames(x[[1]])[[2]]
    } else if (method=="bootstrap"){
      # Creating the bootstrap distribution from user-provided processed datasets
      if (!is.null(boot.data)){
        boot.processed <- lapply(boot.data, FUN=function(b) b$pointEstimation)
        boot.weights <- lapply(boot.processed, FUN=function(b){
          lapply(b$l, FUN=function(b.l.i) rep(1, length(b.l.i)))
        })
        # Creating the bootstrap distribution by perturbing individuals
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l)==1)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          boot.weights[[length(boot.weights) + 1]] <- list(rexp(sum(do.call(c, lapply(processed$x, FUN=function(x.s)dim(x.s)[1]))), 1))
        }
        # Creating the bootstrap distribution by perturbing clusters
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l) > 1)){
        warning("Use of the perturbation bootstrap is discouraged when the clustering and stratification factors coincide.")
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          boot.weights[[length(boot.weights) + 1]] <- lapply(processed$l, FUN=function(b.l.s){
            b.weight.s <- rexp(1, 1)
            rep(b.weight.s, length(b.l.s))
          })
        }
        # Creating the bootstrap distribution by perturbing individuals
      } else if (boot.args$perturb==TRUE & length(var.mapping$stratum)==1){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          b.weights.vec <- rexp(sum(do.call(c, lapply(processed$x, FUN=function(x.s)dim(x.s)[1]))), 1)
          b.weights <- list(); n.s.vec <- c(0, do.call(c, lapply(processed$l, FUN=function(b.l.s)length(b.l.s))))
          for (s in 1:length(processed$l)){
            b.weights[[length(b.weights) + 1]] <- b.weights.vec[(n.s.vec[s] + 1):n.s.vec[s+1]]
          }
          boot.weights[[length(boot.weights) + 1]] <- b.weights
        }
        # Creating the bootstrap distribution by perturbing clusters
      } else if (boot.args$perturb==TRUE & (length(var.mapping$stratum) > 1)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          boot.processed[[length(boot.processed) + 1]] <- processed
          b.weights <- lapply(processed$l, FUN=function(l.s) rep(NA, length(l.s)))
          b.weights.i <- rexp(length(var.mapping$stratum), 1)
          for (i in 1:length(var.mapping$stratum)){
            for (j in 1:length(var.mapping$stratum[[i]])){
              b.weights[[var.mapping$stratum[[i]][j]]][var.mapping$indiv[[i]][j]] <- b.weights.i[i]
            }
          }
          boot.weights[[length(boot.weights) + 1]] <- b.weights
        }
        # Creating the bootstrap distribution by resampling individuals
      } else if (boot.args$perturb==FALSE & is.null(clus)){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          b.index <- sample(nrow(data), nrow(data), replace=TRUE)
          b.data <- data[b.index, ]
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
        # Creating the bootstrap distribution by resampling individuals within clusters: perturb is FALSE, there are clusters, resample.clusters=FALSE
      } else if (boot.args$perturb==FALSE & !is.null(clus) & boot.args$resample.clusters==FALSE){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          data.list <- split(data, data[, as.character(clus)])
          b.data.list <- lapply(data.list, FUN=function(b.d){
            b.index <- sample(nrow(b.d), nrow(b.d), replace=TRUE)
            b.d[b.index, ]
          })
          b.data <- do.call(rbind, b.data.list)
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
        # Creating the bootstrap distribution by resampling clusters: perturb is FALSE, there are clusters, resample.clusters=TRUE
      } else if (boot.args$perturb==FALSE & !is.null(clus) & boot.args$resample.clusters==TRUE){
        boot.processed <- list(); boot.weights <- list()
        for (b in 1:boot.args$n.boot){
          data.list <- split(data, data[, as.character(clus)])
          b.index <- sample(length(data.list), length(data.list), replace=TRUE)
          b.data.list <- data.list[b.index]
          b.data <- do.call(rbind, b.data.list)
          b.temp <- process_data(formula, b.data, strata, id, clus)
          boot.processed[[length(boot.processed) + 1]] <- b.temp$pointEstimation
          boot.weights[[length(boot.weights) + 1]] <- lapply(b.temp$pointEstimation$l, FUN=function(b.l.s) rep(1, length(b.l.s)))
        }
      }
      beta.boot <- list()
      for (b in 1:length(boot.processed)){
        b.l <- boot.processed[[b]]$l; b.u <- boot.processed[[b]]$u; b.u.star <- boot.processed[[b]]$u.star; b.tau <- boot.processed[[b]]$tau
        b.risk <- boot.processed[[b]]$risk
        b.x <- boot.processed[[b]]$x
        b.init <- em_init(dim(b.x[[1]])[2], b.tau)
        b.beta.current <- b.init$beta; b.lambda.current <- b.init$lambda
        b.l1.norm <- 1; b.n.iter <- 0
        while (b.l1.norm > args$tol & b.n.iter < args$maxit.var){
          b.update <- point_iter(b.beta.current, b.lambda.current, b.l, b.u, b.u.star, b.tau, b.x, b.risk, weights=boot.weights[[b]])
          b.l1.norm <- sum(abs(b.update$beta - b.beta.current))
          b.beta.current <- b.update$beta; b.lambda.current <- b.update$lambda
          b.n.iter <- b.n.iter + 1
        }
        beta.boot[[length(beta.boot)+1]] <- b.beta.current
      }
      beta.cov <- var(do.call(rbind, beta.boot))
      colnames(beta.cov) <- dimnames(x[[1]])[[2]]
    }
  } else {
    beta.cov <- NULL
  }
  var.end <- Sys.time()

  # Formatting the final results
  lambda.mat <- lapply(seq_len(length(lambda.current)), FUN=function(x){
    temp <- cbind(tau[[x]], lambda.current[[x]])
    colnames(temp) <- c("tau.sr", "lambda.sr")
    return(temp)
  })
  names(lambda.mat) <- names(lambda.current)
  names(beta.current) <- dimnames(x[[1]])[[2]]
  trace.out <- NULL
  if (args$trace){
    trace.lambda.2 <- lapply(seq_len(length(trace.lambda[[1]])), FUN=function(x)lapply(trace.lambda, `[[`, x))
    names(trace.lambda.2) <- names(trace.lambda[[1]])
    trace.out <- list('beta'=trace.beta, 'lambda'=trace.lambda.2)
  }
  if (variance){
    if (method=="profile"){
      if (is.null(var.mapping) & length(processed$l)==1){
        M <- 1
      } else if (is.null(var.mapping) & length(processed$l) > 1){
        M <- length(processed$l)
      } else if (length(var.mapping$stratum)==1){
        M <- 1
      } else {
        M <- length(var.mapping$stratum)
      }
      var.args <- list('c'=args$c)
    } else if (method=="bootstrap"){
      if (!is.null(boot.data)){
        M <- NA
        var.args <- list('n.boot'=length(boot.data))
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l)==1)){
        M <- 1
        var.args <- list('n.boot'=boot.args$n.boot)
      } else if (boot.args$perturb==TRUE & (is.null(var.mapping) & length(processed$l) > 1)){
        M <- length(processed$l)
        var.args <- list('n.boot'=boot.args$n.boot)
      } else if (boot.args$perturb==TRUE & (length(var.mapping$stratum)==1)){
        M <- 1
        var.args <- list('n.boot'=boot.args$n.boot)
      } else if (boot.args$perturb==TRUE & (length(var.mapping$stratum) > 1)){
        M <- length(var.mapping$stratum)
        var.args <- list('n.boot'=boot.args$n.boot)
      } else if (boot.args$perturb==FALSE & is.null(clus)){
        M <- 1
        var.args <- list('n.boot'=boot.args$n.boot)
      } else if (boot.args$perturb==FALSE & !is.null(clus)){
        M <- length(unique(data[, as.character(clus)]))
        var.args <- list('n.boot'=boot.args$n.boot)
      }
    }
  } else{
    var.args <- M <- NA
  }
  fit.time <- difftime(fit.end, fit.start, units="min")
  var.time <- ifelse(variance==TRUE, difftime(var.end, var.start, units="min"), NA)

  out <- list(call=call, data.attr=list('n'=sum(do.call(c, lapply(x, FUN=function(dat)dim(dat)[1]))),
                                        'M'= M, 'S'=length(lambda.mat),
                                        'events'=sum(do.call(c, lapply(u, FUN=function(right)right != Inf)))),
              coefficients=beta.current, var=beta.cov, var.method=method, var.args=var.args,
              baseline.hazard=lambda.mat, n.iter=n.iter, l1.norm=l1.norm, trace=trace.out,
              fit.times=list('model'=fit.time, 'variance'=var.time))
  class(out) <- 'compCoxIC'
  return(out)
}



