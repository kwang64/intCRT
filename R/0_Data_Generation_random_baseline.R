#' Generate Baseline Survival, Density, and Hazard Functions
#'
#' Constructs piecewise-constant approximations of the baseline survival, density, and hazard functions
#' using a monotonic cubic spline interpolation of randomly selected anchor points.
#'
#' @param maxT Numeric. The maximum observed time (end of follow-up).
#' @param knots Integer. Number of internal knots to use when generating random anchor points.
#' @param max.grid.width Numeric. Maximum spacing between grid points for time discretization. Default is 0.001.
#' @param min.grid.number Integer. Minimum number of grid points. Default is 16000.
#'
#' @return A list with three data frames:
#' \describe{
#'   \item{survival}{Data frame with columns \code{time} and \code{surv}, representing the baseline survival function.}
#'   \item{density}{Data frame with columns \code{time} and \code{density}, representing the baseline density function.}
#'   \item{hazard}{Data frame with columns \code{time} and \code{hazard}, representing the baseline hazard function.}
#' }
#'
#' @examples
#' result <- return_baseline(maxT = 5, knots = 3)
#' head(result$survival)
#' head(result$hazard)
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

#' Generate Exponential Baseline Survival, Density, and Hazard Functions
#'
#' Constructs piecewise-constant approximations of the baseline survival, density, and hazard functions
#' under an exponential hazard model with a Gaussian frailty term.
#'
#' @param maxT Numeric. Upper bound of the support for the failure time distribution.
#' @param lambda0 Numeric. Baseline hazard rate.
#' @param alpha Numeric. Variance of the Gaussian random effect (frailty) on the log-hazard scale.
#' @param max.grid.width Numeric. Maximum width between grid points in the time grid. Default is 0.001.
#' @param min.grid.number Integer. Minimum number of intervals in the time grid. Default is 16000.
#'
#' @return A list with three data frames:
#' \describe{
#'   \item{survival}{Data frame with columns \code{time} and \code{surv}, representing the baseline survival function.}
#'   \item{density}{Data frame with columns \code{time} and \code{density}, representing the baseline density function.}
#'   \item{hazard}{Data frame with columns \code{time} and \code{hazard}, representing the baseline hazard function.}
#' }
#'
#' @details The function simulates a frailty-adjusted constant hazard using a single draw \eqn{b \sim N(0, \alpha)},
#' leading to a hazard function of the form \eqn{\lambda_0 \exp(b)}.
#'
#' @examples
#' set.seed(1)
#' result <- return_exp_baseline(maxT = 10, lambda0 = 0.3, alpha = 0.2)
#' head(result$survival)
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


#' Compute Gamma Term for Positive Stable Distribution
#'
#' Derives the expression for the scaling constant \eqn{\gamma(\alpha)} in the characteristic function
#' of a positive stable distribution, based on the index parameter \eqn{\alpha}.
#'
#' @param a Numeric. The index parameter \eqn{\alpha} (0 < \eqn{\alpha} < 1) of the desired positive stable distribution.
#'
#' @return Numeric. The value of \eqn{\gamma(\alpha)}, computed as \eqn{|1 - i \tan(\pi \alpha / 2)|^{-1/\alpha}}.
#'
#' @details This function is used in the characteristic function formulation of a positive \eqn{\alpha}-stable distribution.
#' It relies on complex arithmetic involving the imaginary unit \eqn{i = \sqrt{-1}}.
#'
#' @examples
#' g(0.5)
#' g(0.8)
#'
#' @export

g <- function(a) {
  iu <- complex(real=0, imaginary=1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}
