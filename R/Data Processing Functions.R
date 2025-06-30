### Processes the input dataset

#' Process Stratified Survival Data
#' Main wrapper that processes a dataset with optional stratification and clustering.
#'
#' @param formula formula object for the stratified COX model
#' @param data full set of observations
#' @param strata vector indicating the stratification variables
#' @param id individual subject identification variable
#' @param clus vector indicating the clustering variable
#'
#' @return A list with two components
#' \describe{
#'   \item{pointEstimation}{A list containing processed data objects (like time intervals, risk sets, and covariates) ready for survival model estimation, organized by strata.}
#'   \item{varMapping}{A list describing how clusters map to strata and individuals for variance estimation; or NULL if no clustering or stratification is used.}
#' }
#' @export
#'
process_data <- function(formula, data, strata, id, clus){

  row.names(data) <- 1:nrow(data)

  # Creating list of datasets for each stratum
  if (!is.null(strata)){
    data.est <- split(data, data[, as.character(strata)])
  } else{
    data.est <- list(data)
  }

  # Creating processed data structures for each stratum
  processed.list <- lapply(data.est, FUN=function(d)process_data_s(formula, d, id))

  # Forming processed data object for point estimation
  processed <- list(l=lapply(processed.list, '[[', 'l.s'),
                    u=lapply(processed.list, '[[', 'u.s'),
                    u.star=lapply(processed.list, '[[', 'u.s.star'),
                    tau=lapply(processed.list, '[[', 'tau.s'),
                    risk=lapply(processed.list, '[[', 'risk.vr'),
                    x=lapply(processed.list, '[[', 'x.s'))

  if (is.null(clus) & is.null(strata)){

    # Unstratified model with model-based variance estimation
    mapping <- NULL

  } else if (is.null(clus)){

    # Creating list of datasets for each cluster
    data.clus <- list(data)

    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))

    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))

  } else if (is.null(strata)){

    # Creating list of datasets for each cluster
    data.clus <- split(data, data[, as.character(clus)])

    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))

    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))

  } else if (clus == strata){

    # Stratification and clustering variable(s) coincide
    mapping <- NULL

  } else {

    # Creating list of datasets for each cluster
    data.clus <- split(data, data[, as.character(clus)])

    # Creating the cluster-to-stratum mappings for each cluster
    mapping.list <- lapply(data.clus, FUN=function(d)cluster_to_stratum_mapping_i(d, data.est, strata, id))

    # Forming the final mapping object
    mapping <- list(stratum=lapply(mapping.list, '[[', 'stratum.id'),
                    indiv=lapply(mapping.list, '[[', 'indiv.id'))

  }

  return(list(pointEstimation=processed, varMapping=mapping))
}

#' Process Data for a Single Stratum in Survival Analysis
#'
#' This function processes the input data for one stratum, extracting relevant
#' time intervals, risk sets, and covariate matrices needed for stratified Cox modeling.
#'
#' @param formula A formula object specifying the survival model, typically
#'   with start and stop times on the left side.
#' @param data A data frame containing the observations for the given stratum.
#' @param id Optional character string specifying the subject ID variable
#'   for handling time-varying covariates. If NULL, assumes no time-varying covariates.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{l.s}{Vector of start times for each individual in the stratum.}
#'   \item{u.s}{Vector of stop times (event or censoring times) for each individual.}
#'   \item{u.s.star}{Modified stop times where infinite times are replaced by start times.}
#'   \item{tau.s}{Sorted unique event or censoring times within the stratum.}
#'   \item{risk.vr}{A matrix indicating risk sets for each time point (rows: individuals, columns: time points).}
#'   \item{x.s}{A 3-dimensional array of covariate values for each individual and time point.}
#' }
#'
#' @export

process_data_s <- function(formula, data, id){

  # Determining relevant column names and indices
  vars <- all.vars(formula); l <- vars[1]; u <- vars[2]

  # Identifying one unique row for each individual
  if (!is.null(id)){
    id.s <- data[, as.character(id)]
    row.id.s <- match(unique(id.s), id.s)
  } else{
    row.id.s <- 1:nrow(data)
  }

  # Creating outcome, time, and riskset variables
  l.s <- data[row.id.s, l]; u.s <- data[row.id.s, u]
  u.s.star <- ifelse(u.s==Inf, l.s, u.s)
  tau.s <- sort(unique(c(l.s[l.s > 0], u.s[u.s < Inf])))
  risk.vr <- do.call(rbind, lapply(u.s.star, FUN=function(x)as.numeric(x >= tau.s)))

  # Creating covariate array
  mat <- model.matrix(formula[c(1, 3)], data=data)
  mat <- mat[, -which(colnames(mat)=="(Intercept)"), drop=FALSE]
  if (is.null(id)){
    x.s <- array(rep(mat, length(tau.s)), dim=c(nrow(mat), ncol(mat), length(tau.s)))
  } else{
    temp <- matrix(NA, nrow=length(unique(id.s))*ncol(mat), ncol=length(tau.s))
    temp.index <- seq(1, length(unique(id.s))*(ncol(mat)), by=length(unique(id.s)))
    for (v in unique(id.s)){
      index <- which(id.s==v)
      temp.time <- data[index, c("start", "stop")]
      temp.mat <- mat[index, ]
      temp[temp.index, ] <- do.call(cbind, lapply(tau.s, FUN=function(r){
        if (length(index)==1){
          return(temp.mat)
        } else{
          cov.index <- sum(temp.time[, "start"] < r)
          if (!is.matrix(temp.mat)){
            return(temp.mat[cov.index])
          } else{
            return(temp.mat[cov.index, ])
          }
        }
      }))
      temp.index <- temp.index + 1
    }
    x.s <- array(temp, dim=c(length(unique(id.s)), ncol(mat), length(tau.s)))
  }

  dimnames(x.s)[[2]] <- as.list(colnames(mat))

  # Returning the processed data structures for stratum s
  return(list(l.s=l.s, u.s=u.s, u.s.star=u.s.star, tau.s=tau.s, risk.vr=risk.vr, x.s=x.s))

}

#' Map Cluster Members to Their Corresponding Strata and Individuals
#'
#' This function determines how individuals within a given cluster correspond to
#' strata and their positions within those strata. It creates mappings needed
#' for variance estimation in stratified and clustered survival analysis.
#'
#' @param data A data frame containing observations for one cluster.
#' @param stratum.list A list of data frames, each representing a stratum with processed data.
#' @param strata Optional character string specifying the stratification variable.
#' @param id Optional character string specifying the individual subject ID variable.
#'
#' @return A list with two components:
#' \describe{
#'   \item{stratum.id}{An integer vector indicating the stratum membership of each individual in the cluster.}
#'   \item{indiv.id}{An integer vector indicating each individual's row index within their corresponding stratum data frame.}
#' }
#'
#' @export
cluster_to_stratum_mapping_i <- function(data, stratum.list, strata, id){

  # Identifying one unique row for each individual
  if (!is.null(id)){
    id.i <- data[, as.character(id)]
    row.id.i <- match(unique(id.i), id.i)
  } else{
    row.id.i <- 1:nrow(data)
  }

  # Mapping each individual in cluster i to their corresponding stratum
  if (is.null(strata)){
    stratum.id <- rep(1, nrow(data[row.id.i, ]))
  } else {
    stratum.id <- apply(data[row.id.i, as.character(strata), drop=FALSE], 1, FUN=function(r){
      which(names(stratum.list) == paste0(as.character(r), collapse="."))
    })
  }

  # Mapping each individual in cluster i to their individual index in the corresponding stratum
  indiv.id <- mapply(function(row.name.ij, stratum.ij){
    strat <- stratum.list[[stratum.ij]]
    if (!is.null(id)){
      strat.id <- strat[, as.character(id)]
      strat.row.id <- match(unique(strat.id), strat.id)
    } else {
      strat.row.id <- 1:nrow(strat)
    }
    which(row.names(strat[strat.row.id, ]) == row.name.ij)
  }, row.name.ij = row.names(data[row.id.i, ]), stratum.ij = stratum.id)


  return(list(stratum.id=stratum.id, indiv.id=indiv.id))

}
