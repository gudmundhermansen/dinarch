test_that("dinarch_fit_ml recovers approximately correct parameters on simulated data", {
  n <- 500
  true_beta <- log(4)
  true_b    <- 0.4
  true_phi  <- 8

  sim <- dinarch_simulate(n = n, b = true_b, phi = true_phi, beta = true_beta, seed = 99)
  dat <- data.frame(y = sim$y, index = sim$index)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1,
                         optim_method = "Nelder-Mead")

  expect_s3_class(fit, "dinarch_fit")
  expect_equal(fit$coefficients$b, true_b, tolerance = 0.15)
  expect_equal(exp(unname(fit$coefficients$beta["(Intercept)"])), exp(true_beta), tolerance = 1)
  expect_gt(fit$coefficients$phi, 1)
})

test_that("dinarch_fit_ml works with a formula covariate that can take negative values", {
  n <- 400
  gdp <- rnorm(n)  # standardized-style covariate, negative values allowed
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 6,
                           formula = ~gdp, covariates = data.frame(gdp = gdp),
                           beta = c(log(3), 0.5), seed = 21)
  dat <- data.frame(y = sim$y, index = seq_len(n), gdp = gdp)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", formula = ~gdp,
                         n_lags = 1, optim_method = "Nelder-Mead")

  expect_setequal(names(fit$coefficients$beta), c("(Intercept)", "gdp"))
  expect_false(anyNA(fit$coefficients$beta))
  expect_equal(unname(fit$coefficients$beta["gdp"]), 0.5, tolerance = 0.3)
})

test_that("multiple covariates (including interactions) fit correctly and in order", {
  n <- 400
  cov_a <- runif(n, 1, 3)
  cov_b <- runif(n, 1, 3)
  true_beta <- c(log(2), log(0.5), 0.2)
  sim <- dinarch_simulate(
    n = n, b = 0.1, phi = 8,
    formula = ~ cov_a + cov_b, covariates = data.frame(cov_a = cov_a, cov_b = cov_b),
    beta = true_beta, seed = 3
  )
  dat <- data.frame(y = sim$y, index = sim$index, cov_a = cov_a, cov_b = cov_b)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", formula = ~ cov_a + cov_b, n_lags = 1)

  expect_identical(names(fit$coefficients$beta), c("(Intercept)", "cov_a", "cov_b"))
  expect_equal(unname(fit$coefficients$beta), true_beta, tolerance = 0.4)
})

test_that("formula = ~0 + ... drops the intercept", {
  n <- 200
  cov_a <- runif(n, 1, 3)
  dat <- data.frame(y = rpois(n, 3), index = seq_len(n), cov_a = cov_a)

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", formula = ~ 0 + cov_a, n_lags = 1)

  expect_identical(names(fit$coefficients$beta), "cov_a")
})

test_that("optim_restarts is reported and control/max_restarts are accepted", {
  n <- 150
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 5)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_ml(
    dat, y = "y", index = "index", n_lags = 1,
    control = list(maxit = 500), max_restarts = 1
  )

  expect_true(is.numeric(fit$optim_restarts))
  expect_gte(fit$optim_restarts, 0)
  expect_lte(fit$optim_restarts, 1)
})

test_that("n_starts > 1 does multi-start optimization and reports the best result", {
  n <- 150
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 5)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, n_starts = 3, seed = 1)

  expect_identical(fit$n_starts, 3)
  expect_false(is.na(fit$loglik))
})

test_that("dinarch_fit_ml converges for high n_lags instead of getting stuck at the start", {
  # Regression test: the old default start (rep(-1, n_lags) on the raw b logit
  # scale) put sum(b) >= 1 at the very first point once n_lags >= 4, landing in
  # the infeasible flat region with no gradient signal - the fit would not
  # move at all from its starting values. The stick-breaking
  # reparameterization removes that infeasible region entirely.
  n <- 400
  n_lags <- 6
  true_b <- c(0.3, rep(0, n_lags - 1))
  sim <- dinarch_simulate(n = n, b = true_b, phi = 6, beta = log(3), seed = 7)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = n_lags)

  expect_lt(sum(fit$coefficients$b), 1)
  # sum(b) is much better identified than the individual b_j with 6
  # collinear lags on modest n, so check the total rather than b[1] alone.
  expect_equal(sum(fit$coefficients$b), sum(true_b), tolerance = 0.3)
  expect_false(is.na(fit$loglik))
})

test_that("multiple groups are pooled correctly (lags never cross group boundaries)", {
  n <- 150
  sim1 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 1)
  sim2 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 2)

  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n), group = "B")
  )

  fit <- dinarch_fit_ml(dat, y = "y", index = "index", group = "group", n_lags = 1)

  expect_equal(fit$n_obs, 2 * (n - 1))  # one lag lost per group
})

test_that("missing required columns raise an informative error", {
  dat <- data.frame(yy = 1:10, index = 1:10)
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1),
    "missing required column"
  )
})

test_that("a formula referencing a nonexistent covariate errors", {
  dat <- data.frame(y = 1:10, index = 1:10)
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", formula = ~cov_a, n_lags = 1),
    "cov_a"
  )
})

test_that("a mis-sized `start` raises an informative error", {
  dat <- data.frame(y = rpois(20, 3), index = 1:20)
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, start = c(0, 0)),
    "start"
  )
})

test_that("print and coef methods work on a dinarch_fit object", {
  dat <- data.frame(y = rpois(50, 3), index = 1:50)
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1)

  expect_output(print(fit), "dinarch_fit")
  expect_identical(coef(fit), fit$coefficients)
})
