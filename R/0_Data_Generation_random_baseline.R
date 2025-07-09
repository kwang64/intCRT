#' Generate Random Baseline Survival and Hazard Functions
#'
#' Randomly generates a non-increasing baseline survival function \eqn{S(t)}
#' on the interval \eqn{[0, T]} and returns the corresponding hazard and density
#' functions.
#'
#' The survival function is constructed via a monotonic cubic spline interpolation
#' through randomly selected anchor points.
#'
#' @param maxT Numeric. Upper bound of the support for the failure time distribution.
#' @param knots Integer. Number of internal knots used for spline representation of the survival function.
#' @param max.grid.width Numeric. Maximum size of intervals in the time grid (default = 0.001).
#' @param min.grid.number Integer. Minimum number of intervals in the time grid (default = 16000).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{survival}{A data frame with time points and corresponding survival probabilities \eqn{S(t)}.}
#'   \item{density}{A data frame with time points and corresponding density values \eqn{f(t) = \lambda(t) S(t)}.}
#'   \item{hazard}{A data frame with time points and corresponding hazard values \eqn{\lambda(t)}.}
#' }
#'
#' @examples
#' result <- return_baseline(maxT = 10, knots = 5)
#' head(result$survival)
#'
#' @export

return_baseline <- function(maxT, knots, max.grid.width=0.001, min.grid.number=16000){

  # Selecting random anchor points for the construction of the baseline functions
  partition <- c(0, sort(runif(knots, 0, maxT)), maxT)
  surv <- c(1, sort(runif(knots), decreasing = TRUE), 0)

  # Fitting a monotonic cubic spline model to the anchor points
  survSpline <- stats::splinefun(x=partition, y=surv, method="hyman")

  # Piecewise-constant approximations to the baseline survival, density, and hazard funcitons
  timeGrid <- seq(0, maxT, by=min(max.grid.width, maxT/min.grid.number))
  survGrid <- survSpline(timeGrid)
  pdfGrid <- -survSpline(timeGrid, 1)
  hazGrid <- pdfGrid/survGrid

  # Quality assurance on the estimated components
  if (sum(survGrid > 1) > 0){
    index <- which(survGrid > 1)
    survGrid[index] <- 1
  }
  if (sum(survGrid < 0) > 0){
    index <- which(survGrid < 0)
    survGrid[index] <- 0
  }
  if (sum(is.na(hazGrid)) > 0){
    index  <- which(is.na(hazGrid))
    hazGrid[index] <- 0
  }
  if (sum(hazGrid < 0) > 0){
    index <- which(hazGrid < 0)
    hazGrid[index] <- 0
  }

  return(list(survival=data.frame(time=timeGrid, surv=survGrid),
              density=data.frame(time=timeGrid, density=pdfGrid),
              hazard=data.frame(time=timeGrid, hazard=hazGrid)))
}

#' Generate Exponential Baseline Survival and Hazard Functions
#'
#' Generates baseline survival, hazard, and density functions under the
#' exponential model with stratum-specific hazard variation.
#'
#' The hazard is modeled as \eqn{\lambda(t) = \lambda_0 e^b}, where \eqn{b \sim N(0, \alpha)}.
#'
#' @param maxT Numeric. Upper bound of the support for the failure time distribution.
#' @param lambda0 Numeric. Baseline hazard rate.
#' @param alpha Numeric. Variance of the random effect added to the log-hazard.
#' @param max.grid.width Numeric. Maximum size of intervals in the time grid (default = 0.001).
#' @param min.grid.number Integer. Minimum number of intervals in the time grid (default = 16000).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{survival}{A data frame with time points and corresponding survival probabilities \eqn{S(t)}.}
#'   \item{density}{A data frame with time points and corresponding density values \eqn{f(t) = \lambda(t) S(t)}.}
#'   \item{hazard}{A data frame with time points and corresponding hazard values \eqn{\lambda(t)}.}
#' }
#'
#' @examples
#' result <- return_exp_baseline(maxT = 10, lambda0 = 0.1, alpha = 0.5)
#' head(result$hazard)
#'
#' @export

return_exp_baseline <- function(maxT, lambda0, alpha, max.grid.width=0.001, min.grid.number=16000){

  timeGrid <- seq(0, maxT, by=min(max.grid.width, maxT/min.grid.number))
  b <- rnorm(1, 0, sqrt(alpha))
  hazGrid <- rep(lambda0*exp(b), length(timeGrid)-1)
  survGrid <- c(1, exp(-cumsum(hazGrid[1:(length(hazGrid)-1)]*diff(timeGrid[-length(hazGrid)]))))
  pdfGrid <- hazGrid*survGrid

  return(list(survival=data.frame(time=timeGrid, surv=c(survGrid, 0)),
              density=data.frame(time=timeGrid[-length(timeGrid)], density=pdfGrid),
              hazard=data.frame(time=timeGrid[-length(timeGrid)], hazard=hazGrid)))
}


#' Positive Stable Distribution Constant
#'
#' Computes the normalization constant \eqn{\gamma} in the positive stable distribution
#' with index parameter \eqn{\alpha}.
#'
#' @param a Numeric. Index parameter \eqn{\alpha} of the desired positive stable distribution (0 < a < 1).
#'
#' @return Numeric. Value of the constant \eqn{\gamma}.
#'
#' @examples
#' g(0.5)
#'
#' @export

g <- function(a) {
  iu <- complex(real=0, imaginary=1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}
