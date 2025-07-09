#'
#'
gen_visits_tx <- function(x){
  if (x==1){
    temp <- cumsum(runif(15, 12, 24))
    temp <- c(temp, sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4)))
    temp <- sort(temp)
    return(temp)
  } else if (x==0){
    return(sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4)))
  }
}

#'
#'
gen_visits_bcpp <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4))))
  return(visits)
}

#'
#'
gen_visits_freq <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)cumsum(runif(20, 0, 16))))
  return(visits)
}
