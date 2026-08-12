test_that("dinarch_fit_bayes recovers roughly plausible parameters", {
  testthat::skip_if_not_installed("rstan")

  n <- 200
  true_a   <- log(4)
  true_b   <- 0.3
  true_phi <- 8

  sim <- dinarch_simulate(n = n, a = true_a, b = true_b, phi = true_phi, seed = 55)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_bayes(
    dat, y = "y", index = "index", k = 1,
    iter = 1000, chains = 2, seed = 1
  )

  expect_s3_class(fit, "dinarch_fit")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$coefficients$b, true_b, tolerance = 0.2)
  expect_equal(exp(fit$coefficients$a), exp(true_a), tolerance = 2)
})

test_that("dinarch_fit_bayes returns the expected structure", {
  testthat::skip_if_not_installed("rstan")

  dat <- data.frame(y = rpois(100, 3), index = 1:100)
  # Deliberately few iterations here since this test only checks structure,
  # not estimation quality - low-ESS warnings are expected and harmless.
  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", k = 1,
    iter = 500, chains = 2, seed = 2
  ))

  expect_true(all(c("a", "b", "phi", "beta") %in% names(fit$coefficients)))
  expect_length(fit$coefficients$b, 1)
  expect_s4_class(fit$stanfit, "stanfit")
  expect_true(is.list(fit$posterior))
})

test_that("dinarch_fit_bayes pools multiple groups without crossing group boundaries", {
  testthat::skip_if_not_installed("rstan")

  # group = a series identity (e.g. country), index = order within that
  # group (e.g. year). Lags must never be taken across a group boundary,
  # so one observation per group is lost to the k = 1 lag, regardless of
  # how many groups there are.
  n <- 60
  sim1 <- dinarch_simulate(n = n, a = log(3), b = 0.3, phi = 6, seed = 10)
  sim2 <- dinarch_simulate(n = n, a = log(3), b = 0.3, phi = 6, seed = 11)

  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n), group = "B")
  )

  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", group = "group", k = 1,
    iter = 500, chains = 2, seed = 3
  ))

  expect_equal(fit$n_obs, 2 * (n - 1))
})
