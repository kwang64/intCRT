test_that("w_sv return the vector successfully with correct length and structure", {
  l.sv <- 1
  u.sv <- 3
  tau.s <- c(0.5, 1, 2, 3,4)
  # 2 covariates and 5 time points
  x.pr <- matrix(c(1,0,
                   1,1,
                   0,1,
                   1,-1,
                   0,0),nrow=2)
  lambda.s <- rep(0.2, length(tau.s))
  beta <- c(0.1, -0.2)

  result <- w_sv(l.sv, u.sv, tau.s, x.pr, lambda.s,beta)

  #check length
  expect_length(result, length(tau.s))

  #check vector initialization (values before l.sv and after u.sv are zeros)
  expect_equal(result[1],0)
  expect_equal(result[2],0)
  expect_equal(result[5],0)

  #check the value in (1,3] are positive and correctly normalized
  in_window <- 3
  r.ids <- c(FALSE, FALSE, TRUE, TRUE, FALSE)
  lin.pred <- beta %*% x.pr[, r.ids]
  denom <- 1 - exp(-sum(lambda.s[r.ids] * exp(lin.pred)))
  expected <- lambda.s[r.ids] * exp(lin.pred) / denom
  expect_equal(result[3:4], as.numeric(expected))
})
