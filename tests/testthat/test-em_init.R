test_that("em_init initialize coefficients and terms successfully", {
  p <- 3
  tau <- list(
    c(1,2,3),
    c(3,4)
  )
  result <- em_init(p,tau)

  expect_type(result$beta, "double")
  expect_equal(length(result$beta),p)
  expect_equal(result$beta, rep(0,p))

  expect_type(result$lambda, "list")
  expect_equal(length(result$lambda), length(tau))

  #unit test to verify the lambda output is correctly initialized
  for (i in seq_along(tau)) {
    expect_equal(length(result$lambda[[i]]), length(tau[[i]]))
    expect_equal(result$lambda[[i]], rep(1 / length(tau[[i]]), length(tau[[i]])))
  }
})
