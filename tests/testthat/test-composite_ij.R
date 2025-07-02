test_that("composite_ij returns correct individual log-likelihood", {
  tau.s <- c(1,2,3,4)
  x.pr <- matrix(c(1,2,1,2,
                   2,1,2,1), nrow = 2, byrow = TRUE)
  lambda.s <- rep(0.1, length(tau.s))
  beta <- c(0.5, -0.2)
  l.ij <- 1

  #for finite u.ij
  u.ij <- 3
  result1 <- composite_ij(l.ij, u.ij, tau.s, x.pr, lambda.s, beta)

  expect_type(result1, "double")
  expect_length(result1, 1)

  #compute the expected value manually
  r1<- (tau.s <= l.ij) #tau is 1
  r2 <- (tau.s > l.ij & tau.s <= u.ij) #tau is 2,3
  term1 <- -sum(lambda.s[r1] * exp(beta %*% x.pr[, r1]))
  term2 <- log(1-exp(-sum(lambda.s[r2] * exp(beta %*% x.pr[, r2]))))
  expected1 <- term1 + term2
  expect_equal(result1, expected1)

  #for infinite u.ij (right-censored)
  u.ij.inf <- Inf
  result2 <- composite_ij(l.ij, u.ij.inf, tau.s, x.pr, lambda.s, beta)
  r1_inf <- (tau.s <= l.ij)
  expected2 <- -sum(lambda.s[r1_inf] * exp(beta %*% x.pr[, r1_inf]))
  expect_equal(result2,expected2)
})
