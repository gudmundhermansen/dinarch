test_that("dinarch_fit(method = 'ml') matches dinarch_fit_ml() directly", {
  n <- 300
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 42)
  dat <- data.frame(y = sim$y, index = sim$index)

  fit_direct <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1)
  fit_dispatch <- dinarch_fit(dat, y = "y", index = "index", n_lags = 1, method = "ml")

  expect_s3_class(fit_dispatch, "dinarch_fit")
  expect_equal(fit_dispatch$method, "ml")
  expect_equal(fit_dispatch$coefficients, fit_direct$coefficients)
})

test_that("dinarch_fit() defaults to method = 'ml'", {
  dat <- data.frame(y = rpois(100, 3), index = 1:100)
  fit <- dinarch_fit(dat, y = "y", index = "index", n_lags = 1, vcov = FALSE)
  expect_equal(fit$method, "ml")
})

test_that("dinarch_fit(method = 'bayes') dispatches to dinarch_fit_bayes()", {
  testthat::skip_if_not_installed("rstan")

  dat <- data.frame(y = rpois(100, 3), index = 1:100)
  fit <- suppressWarnings(dinarch_fit(
    dat, y = "y", index = "index", n_lags = 1,
    method = "bayes", iter = 500, chains = 2, seed = 1
  ))

  expect_s3_class(fit, "dinarch_fit")
  expect_equal(fit$method, "bayes")
  expect_s4_class(fit$stanfit, "stanfit")
})

test_that("dinarch_fit() errors on an unknown method", {
  dat <- data.frame(y = rpois(100, 3), index = 1:100)
  expect_error(dinarch_fit(dat, method = "mle"))
})
