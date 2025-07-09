#' Generate Covariate-Dependent Monitoring Schedule
#'
#' Simulates a vector of internal study visit times based on treatment assignment.
#' Participants in the treatment group receive more frequent monitoring, while
#' control group participants receive only scheduled quarterly (annual) visits.
#'
#' @param x Numeric. Treatment assignment indicator (1 = treatment, 0 = control).
#'
#' @return A numeric vector of visit times (in the same time unit as your survival model, e.g., weeks).
#' Treatment group participants receive ~15 random visits spaced between 12 and 24 units,
#' plus 4 regularly scheduled visits (around week 52, 104, 156, 208).
#' Control group participants only receive the 4 regularly scheduled visits.
#'
#' @examples
#' # Simulate visit times for a treated and untreated participant:
#' gen_visits_tx(1)  # More intensive monitoring
#' gen_visits_tx(0)  # Standard quarterly visits
#'
#' @seealso \code{\link{gen_time_dependent_beta}}, \code{\link{gen_individual_tvbeta}}
#'
#' @export

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

#' Generate Yearly Monitoring Schedule with Visit Time Jitter
#'
#' Simulates annual visit schedules for each participant in a study, with small random variation (±4 weeks)
#' around each scheduled yearly visit (e.g., week 52, 104, 156, 208). This function is useful for generating
#' internal inspection times in survival simulations with interval censoring.
#'
#' @param n Numeric. Total number of individuals to simulate visit schedules for.
#'
#' @return A numeric matrix with \code{n} rows and 4 columns, where each row corresponds to an individual and
#' each column represents a jittered annual visit time (in weeks).
#'
#' @examples
#' # Generate visit times for 5 participants
#' gen_visits_bcpp(5)
#'
#' @seealso \code{\link{gen_time_dependent_beta}}, \code{\link{gen_individual_tvbeta}}, \code{\link{gen_visits_tx}}
#'
#' @export
gen_visits_bcpp <- function(n){ ... }

gen_visits_bcpp <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)sapply(1:4, FUN=function(y)runif(1, 52*y-4, 52*y+4))))
  return(visits)
}

#' Generate Frequent Monitoring Visit Schedule
#'
#' Simulates frequent, irregular internal visit times for each individual by generating cumulative
#' sums of random intervals drawn uniformly between 0 and 16 time units (e.g., weeks). This function
#' is useful for modeling intensive follow-up or high-frequency monitoring designs.
#'
#' @param n Numeric. Total number of individuals to simulate visit schedules for.
#'
#' @return A numeric matrix with \code{n} rows and 20 columns. Each row corresponds to one individual,
#' and each column represents a visit time (cumulative over 20 visits).
#'
#' @examples
#' # Generate visit times for 3 individuals under frequent monitoring
#' gen_visits_freq(3)
#'
#' @seealso \code{\link{gen_visits_bcpp}}, \code{\link{gen_visits_tx}}, \code{\link{gen_time_dependent_beta}}
#'
#' @export

gen_visits_freq <- function(n){
  visits <- do.call(rbind, lapply(seq_len(n), FUN=function(x)cumsum(runif(20, 0, 16))))
  return(visits)
}
