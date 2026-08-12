test_that("dinarch_fit_ml recovers approximately correct parameters on simulated data", {
  n <- 500
  true_a   <- log(4)
  true_b   <- 0.4
  true_phi <- 8

  sim <- dinarch_simulate(n = n, a = true_a, b = true_b, phi = true_phi, seed = 99)
  dat <- data.frame(y = sim$y, index = sim$index)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", k = 1,
                         optimizer = "nelder-mead+bfgs")

  expect_s3_class(fit, "dinarch_fit")
  expect_equal(fit$coefficients$b, true_b, tolerance = 0.15)
  expect_equal(exp(fit$coefficients$a), exp(true_a), tolerance = 1)
  expect_gt(fit$coefficients$phi, 1)
})

test_that("dinarch_fit_ml works with a population covariate", {
  n <- 400
  pop <- rep(1000, n)
  sim <- dinarch_simulate(n = n, a = log(0.01), b = 0.3, phi = 6,
                           population = pop, seed = 21)
  dat <- data.frame(y = sim$y, index = seq_len(n), population = pop)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", population = "population",
                         k = 1, optimizer = "nelder-mead+bfgs")

  expect_true(fit$population_used)
  expect_false(is.na(fit$coefficients$a))
})

test_that("optim_restarts is reported and control/max_restarts are accepted", {
  n <- 150
  sim <- dinarch_simulate(n = n, a = log(3), b = 0.3, phi = 6, seed = 5)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_ml(
    dat, y = "y", index = "index", k = 1,
    control = list(maxit = 500), max_restarts = 1
  )

  expect_true(is.numeric(fit$optim_restarts))
  expect_gte(fit$optim_restarts, 0)
  expect_lte(fit$optim_restarts, 1)
})

test_that("multiple groups are pooled correctly (lags never cross group boundaries)", {
  n <- 150
  sim1 <- dinarch_simulate(n = n, a = log(3), b = 0.3, phi = 6, seed = 1)
  sim2 <- dinarch_simulate(n = n, a = log(3), b = 0.3, phi = 6, seed = 2)

  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n), group = "B")
  )

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", group = "group", k = 1)

  expect_equal(fit$n_obs, 2 * (n - 1))  # one lag lost per group
})

test_that("missing required columns raise an informative error", {
  dat <- data.frame(y = 1:10, index = 1:10)
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", population = "pop"),
    "missing required column"
  )
})

test_that("a mis-sized `start` raises an informative error", {
  dat <- data.frame(y = rpois(20, 3), index = 1:20)
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", k = 1, start = c(0, 0)),
    "start"
  )
})

test_that("print and coef methods work on a dinarch_fit object", {
  dat <- data.frame(y = rpois(50, 3), index = 1:50)
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", k = 1)

  expect_output(print(fit), "dinarch_fit")
  expect_identical(coef(fit), fit$coefficients)
})
