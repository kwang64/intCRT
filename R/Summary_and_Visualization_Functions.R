#' Summarize a Fitted Composite Cox Proportional Hazards Model
#'
#' Prints a summary of the fitted Cox proportional hazards model returned by \code{composite_coxIC()},
#' including data attributes, regression coefficients, hazard ratios, standard errors, z-statistics,
#' p-values, and significance stars.
#'
#' @param x List. Model output from \code{composite_coxIC()} containing estimated coefficients, variance,
#'   variance method, data attributes, and call information.
#'
#' @details
#' The summary prints the total number of subjects, clusters, strata, and events used in the model.
#' Regression coefficients and their standard errors are displayed alongside hazard ratios
#' (\eqn{\exp(\beta)}) with associated z-scores and p-values. Significance levels are annotated as
#' \code{***} (p < 0.001), \code{**} (p < 0.01), \code{*} (p < 0.05), \code{.} (p < 0.1), or blank otherwise.
#' If variance estimates are unavailable, related statistics are shown as NA.
#'
#' Additionally, the method used for variance estimation is described (profile composite likelihood or bootstrap).
#'
#' @return Invisibly returns \code{x} after printing the summary.
#'
#' @export
summary.compCoxIC <- function(x){

  # Summarizing data attributes
  data.attr <- paste0("Total subjects = ", x$data.attr$n, ", clusters = ", x$data.attr$M, ", strata = ", x$data.attr$S, ", events = ", x$data.attr$events)

  # Summarizing estimated regression coefficients
  log.hr <- matrix(NA, nrow=length(x$coefficients), ncol=5)
  signif.ind <- c()
  for (i in 1:nrow(log.hr)){
    log.hr[i, 1] <- format(round(x$coefficients[i], 3), nsmall=3)
    log.hr[i, 2] <- format(round(exp(x$coefficients[i]), 3), nsmall=3)
    if (is.null(x$var)){
      log.hr[i, 3] <- log.hr[i, 4] <- log.hr[i, 5] <- NA
      signif.ind <- c(signif.ind, "")
    } else {
      log.hr[i, 3] <- format(round(sqrt(x$var[i, i]), 3), nsmall=3)
      log.hr[i, 4] <- format(round(x$coefficients[i]/sqrt(x$var[i, i]), 3), nsmall=3)
      log.hr[i, 5] <- format(round(1-pchisq(x$coefficients[i]^2/x$var[i, i], 1), 6), nsmall=6)
      if (as.numeric(log.hr[i, 5]) < 0.001){
        signif.ind <- c(signif.ind, "***")
      } else if (as.numeric(log.hr[i, 5]) < 0.01){
        signif.ind <- c(signif.ind, "**")
      } else if (as.numeric(log.hr[i, 5]) < 0.05){
        signif.ind <- c(signif.ind, "*")
      } else if (as.numeric(log.hr[i, 5]) < 0.1){
        signif.ind <- c(signif.ind, ".")
      } else{
        signif.ind <- c(signif.ind, "")
      }
    }
  }
  log.hr <- as.data.frame(log.hr)
  log.hr$signif <- signif.ind
  rownames(log.hr) <- names(x$coefficients)
  colnames(log.hr) <- c("Coef", "exp(Coef)", "se(Coef)", "z", "Pr(>|z|)", "")

  # Determining variance type
  if (!is.null(x$var)){
    var.statement <- "Note: Variance computed using "
    if (x$var.method=="profile"){
      var.statement <- paste0(var.statement, "the profile composite log-likelihood function with c = ", x$var.args$c, ".")
    } else{
      var.statement <- paste0(var.statement, x$var.args$n.boot, " bootstrap resamples.")
    }
  }

  # Printing out results
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Data Attributes:\n")
  cat(data.attr)
  cat("\n\n")
  cat("Coefficients:\n")
  print(noquote(log.hr))
  if (!is.null(x$var)){
    cat("\n")
    cat(var.statement)
  }
}

#' Plot Estimated Baseline Cumulative Hazard, Survival, or Cumulative Incidence Functions
#'
#' Visualizes the estimated baseline distribution functions from the output of \code{composite_coxIC()}.
#' Users can plot cumulative hazard, survival, or cumulative incidence functions, either combined on one plot
#' or separately by stratum.
#'
#' @param x List. Model output from \code{composite_coxIC()}, containing baseline hazard estimates by stratum.
#' @param type Character vector. Which baseline distribution function(s) to plot. Options include
#'   \code{"cumulative hazard"}, \code{"survival"}, and \code{"cumulative incidence"}. Can specify multiple.
#' @param combine Logical. If \code{TRUE}, plots all strata on a single combined plot. If \code{FALSE},
#'   plots separate plots for each stratum.
#' @param x.title Character. Label for the x-axis (default \code{"Time"}).
#' @param legend.title Character. Title for the legend indicating strata (default \code{"Stratum"}).
#' @param y.title Optional character vector of length matching \code{type} to customize y-axis titles.
#'   If \code{NULL}, default titles will be used.
#'
#' @return NULL. This function prints ggplot2 plots to the current graphics device.
#'
#' @details
#' The function computes cumulative sums of the estimated baseline hazard to derive cumulative hazard,
#' survival (via \eqn{\exp(-\text{cumulative hazard})}), or cumulative incidence (via \eqn{1-\exp(-\text{cumulative hazard})}).
#' Multiple types can be plotted by specifying them in \code{type}.
#'
#' @import ggplot2
#' @export
plot.compCoxIC <- function(x, type=c("cumulative hazard", "survival", "cumulative incidence"), combine=TRUE, x.title="Time", legend.title="Stratum", y.title=NULL){

  if (combine){
    if ("cumulative hazard" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="cumulative hazard")])){
        y.title.2 <- "Cumulative Baseline Hazard"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="cumulative hazard")]
      }
      temp <- do.call(rbind, lapply(seq_len(length(x$baseline.hazard)), FUN=function(z){
        data.frame(tau.ir=x$baseline.hazard[[z]][,1], ch=cumsum(x$baseline.hazard[[z]][,2]), stratum=names(x$baseline.hazard)[z])

      }))
      p <- ggplot2::ggplot(temp, aes(x=tau.ir, y=ch, group=stratum)) + geom_line(aes(color=stratum))
      p <- p + xlab(x.title) + ylab(y.title.2) + guides(color=guide_legend(title=legend.title))
      print(p)
    }
    if ("survival" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="survival")])){
        y.title.2 <- "Baseline Survival"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="survival")]
      }
      temp <- do.call(rbind, lapply(seq_len(length(x$baseline.hazard)), FUN=function(z){
        data.frame(tau.ir=x$baseline.hazard[[z]][,1], surv=exp(-cumsum(x$baseline.hazard[[z]][,2])), stratum=names(x$baseline.hazard)[z])

      }))
      p <- ggplot2::ggplot(temp, aes(x=tau.ir, y=surv, group=stratum)) + geom_line(aes(color=stratum))
      p <- p + xlab(x.title) + ylab(y.title.2) + guides(color=guide_legend(title=legend.title))
      print(p)
    }
    if ("cumulative incidence" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="cumulative incidence")])){
        y.title.2 <- "Baseline Cumulative Incidence"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="cumulative incidence")]
      }
      temp <- do.call(rbind, lapply(seq_len(length(x$baseline.hazard)), FUN=function(z){
        data.frame(tau.ir=x$baseline.hazard[[z]][,1], surv=1-exp(-cumsum(x$baseline.hazard[[z]][,2])), stratum=names(x$baseline.hazard)[z])

      }))
      p <- ggplot2::ggplot(temp, aes(x=tau.ir, y=surv, group=stratum)) + geom_line(aes(color=stratum))
      p <- p + xlab(x.title) + ylab(y.title.2) + guides(color=guide_legend(title=legend.title))
      print(p)
    }
  } else{
    if ("cumulative hazard" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="cumulative hazard")])){
        y.title.2 <- "Cumulative Baseline Hazard"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="cumulative hazard")]
      }
      colors <- rainbow(length(x$baseline.hazard))
      for (i in 1:length(x$baseline.hazard)){
        stratum.name <- names(x$baseline.hazard)[i]
        temp <- data.frame(tau.ir=x$baseline.hazard[[i]][,1], ch=cumsum(x$baseline.hazard[[i]][,2]))
        p <- ggplot2::ggplot(temp, aes(x=tau.ir, y=ch)) + geom_line(color=colors[i])
        p <- p + xlab(x.title) + ylab(y.title.2) + title(paste0(legend.title," ", stratum.name))
        print(p)
      }
    }
    if ("survival" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="survival")])){
        y.title.2 <- "Baseline Survival"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="survival")]
      }
      colors <- rainbow(length(x$baseline.hazard))
      for (i in 1:length(x$baseline.hazard)){
        stratum.name <- names(x$baseline.hazard)[i]
        temp <- data.frame(tau.ir=x$baseline.hazard[[i]][,1], surv=exp(-cumsum(x$baseline.hazard[[i]][,2])))
        p <- ggplot2::ggplot(temp, aes(x=tau.ir, y=surv)) + geom_line(color=colors[i])
        p <- p + xlab(x.title) + ylab(y.title.2) + ggtitle(paste0(legend.title," ", stratum.name))
        print(p)
      }
    }
    if ("cumulative incidence" %in% tolower(type)){
      if (is.null(y.title[which(tolower(type)=="cumulative incidence")])){
        y.title.2 <- "Baseline Cumulative Incidence"
      } else{
        y.title.2 <- y.title[which(tolower(type)=="cumulative incidence")]
      }
      colors <- rainbow(length(x$baseline.hazard))
      for (i in 1:length(x$baseline.hazard)){
        stratum.name <- names(x$baseline.hazard)[i]
        temp <- data.frame(tau.ir=x$baseline.hazard[[i]][,1], surv=1-exp(-cumsum(x$baseline.hazard[[i]][,2])))
        p <-  ggplot2::ggplot(temp, aes(x=tau.ir, y=surv)) + geom_line(color=colors[i])
        p <- p + xlab(x.title) + ylab(y.title.2) + ggtitle(paste0(legend.title," ", stratum.name))
        print(p)
      }
    }
  }

}
