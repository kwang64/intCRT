### Processes the input dataset
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

### Processes the data for stratum s
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

### Determines the mapping between cluster membership and stratum membership
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
