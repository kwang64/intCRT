# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(intCRT)

test_check("intCRT")


test_that("simulation generates valid data structure", {
  set.seed(123)

  baseline <- list(return_baseline(100, 3))

  sim_data <- gen_time_dependent_cov(
    M = 2,
    ni = c(3, 4),
    beta = c(-1, 0.5, 0),
    gen_visits = gen_visits_bcpp,
    cens = 0,
    rho = 0.3,
    baseline = baseline,
    S = function(n) rep(1, n)
  )

  expect_true(is.data.frame(sim_data))
  expect_true(all(c("left", "right", "id") %in% names(sim_data)))
  expect_true(nrow(sim_data) > 0)
})

test_that("baseline generator returns list with valid components", {
  b <- return_baseline(200, 4)

  expect_true(is.list(b))
  expect_true(all(c("time", "hazard") %in% names(b)))
  expect_true(length(b$time) == length(b$hazard))
})

test_that("composite_coxIC fits on tiny simulated dataset", {
  set.seed(123)

  baseline <- list(return_baseline(100, 3))

  data <- gen_time_dependent_cov(
    M = 2,
    ni = c(3, 4),
    beta = c(-1, 0.5, 0),
    gen_visits = gen_visits_bcpp,
    cens = 0,
    rho = 0.3,
    baseline = baseline,
    S = function(n) rep(1, n)
  )

  fit <- composite_coxIC(
    Surv(left, right, type = "interval2") ~ treat + utt,
    data = data,
    strata = "stratum",
    clus = "group",
    id = "id",
    variance = FALSE   # 加快速度
  )

  expect_s3_class(fit, "compCoxIC")
  expect_true(is.numeric(fit$coefficients))
})

test_that("summary method works", {
  set.seed(123)

  baseline <- list(return_baseline(100, 3))

  data <- gen_time_dependent_cov(
    M = 2,
    ni = c(3, 4),
    beta = c(-1, 0.5, 0),
    gen_visits = gen_visits_bcpp,
    cens = 0,
    rho = 0.3,
    baseline = baseline,
    S = function(n) rep(1, n)
  )

  fit <- composite_coxIC(
    Surv(left, right, type = "interval2") ~ treat + utt,
    data = data,
    strata = "stratum",
    clus = "group",
    id = "id",
    variance = FALSE
  )

  expect_no_error(summary(fit))
})

test_that("plot method works", {
  set.seed(123)

  baseline <- list(return_baseline(100, 3))

  data <- gen_time_dependent_cov(
    M = 2,
    ni = c(3, 4),
    beta = c(-1, 0.5, 0),
    gen_visits = gen_visits_bcpp,
    cens = 0,
    rho = 0.3,
    baseline = baseline,
    S = function(n) rep(1, n)
  )

  fit <- composite_coxIC(
    Surv(left, right, type = "interval2") ~ treat + utt,
    data = data,
    strata = "stratum",
    clus = "group",
    id = "id",
    variance = FALSE
  )

  expect_no_error(
    plot(fit, type = "cumulative hazard", combine = TRUE)
  )
})
