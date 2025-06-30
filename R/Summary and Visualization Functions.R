### Summarizes the fitted Cox proportional hazards model
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

### Plots the estimated baseline cumulative hazard (survival) function(s)
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
      p <- ggplot(temp, aes(x=tau.ir, y=ch, group=stratum)) + geom_line(aes(color=stratum))
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
      p <- ggplot(temp, aes(x=tau.ir, y=surv, group=stratum)) + geom_line(aes(color=stratum))
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
      p <- ggplot(temp, aes(x=tau.ir, y=surv, group=stratum)) + geom_line(aes(color=stratum))
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
        p <- ggplot(temp, aes(x=tau.ir, y=ch)) + geom_line(color=colors[i])
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
        p <- ggplot(temp, aes(x=tau.ir, y=surv)) + geom_line(color=colors[i])
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
        p <- ggplot(temp, aes(x=tau.ir, y=surv)) + geom_line(color=colors[i])
        p <- p + xlab(x.title) + ylab(y.title.2) + ggtitle(paste0(legend.title," ", stratum.name))
        print(p)
      }
    }
  }

}
