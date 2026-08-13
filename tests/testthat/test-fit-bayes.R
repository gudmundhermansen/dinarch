test_that("dinarch_fit_bayes recovers roughly plausible parameters", {
  testthat::skip_if_not_installed("rstan")

  n <- 200
  true_beta <- log(4)
  true_b    <- 0.3
  true_phi  <- 8

  sim <- dinarch_simulate(n = n, b = true_b, phi = true_phi, beta = true_beta, seed = 55)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  fit <- dinarch_fit_bayes(
    dat, y = "y", index = "index", n_lags = 1,
    iter = 1000, chains = 2, seed = 1
  )

  expect_s3_class(fit, "dinarch_fit")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$coefficients$b, true_b, tolerance = 0.2)
  expect_equal(exp(unname(fit$coefficients$beta["(Intercept)"])), exp(true_beta), tolerance = 2)
})

test_that("dinarch_fit_bayes enforces sum(b) < 1 for every posterior draw when n_lags > 1", {
  testthat::skip_if_not_installed("rstan")

  n <- 150
  sim <- dinarch_simulate(n = n, b = c(0.3, 0.3), phi = 6, beta = log(3), seed = 44)
  dat <- data.frame(y = sim$y, index = seq_len(n))

  # Deliberately vague per-lag shape1 so that, under the old independent
  # per-lag beta(1, 3) priors, sum(b) >= 1 would be readily reachable for
  # n_lags = 2 - this is the regression case for the Dirichlet reparameterization.
  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", n_lags = 2,
    iter = 500, chains = 2, seed = 6
  ))

  expect_true(all(rowSums(fit$posterior$b) < 1))
})

test_that("dinarch_fit_bayes works with a formula covariate that can take negative values", {
  testthat::skip_if_not_installed("rstan")

  n <- 200
  gdp <- rnorm(n)
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 6,
                           formula = ~gdp, covariates = data.frame(gdp = gdp),
                           beta = c(log(3), 0.4), seed = 8)
  dat <- data.frame(y = sim$y, index = seq_len(n), gdp = gdp)

  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", formula = ~gdp, n_lags = 1,
    iter = 800, chains = 2, seed = 9
  ))

  expect_setequal(names(fit$coefficients$beta), c("(Intercept)", "gdp"))
  expect_false(anyNA(fit$coefficients$beta))
})

test_that("dinarch_fit_bayes does not fail to initialize with large-scale raw covariates", {
  # Regression test: Stan's default random Uniform(-2, 2) init, combined
  # with a covariate on a very different raw scale (e.g. population in
  # the tens of millions), makes exp(eta) overflow to Inf/0 at every
  # attempted initialization ("Initialization failed"), and the resulting
  # stanfit with zero samples used to crash confusingly downstream in
  # colMeans(). The data-driven default init (non-intercept beta always
  # starts at exactly 0, so eta_init never depends on covariate scale)
  # fixes this.
  testthat::skip_if_not_installed("rstan")

  n <- 150
  population <- stats::runif(n, 1e6, 5e7)
  gdppc <- stats::runif(n, 500, 40000)
  dat <- data.frame(y = rpois(n, 5), index = seq_len(n), population = population, gdppc = gdppc)

  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", n_lags = 1,
    formula = ~ 1 + population + gdppc,
    iter = 300, chains = 1, seed = 1
  ))

  expect_s3_class(fit, "dinarch_fit")
  expect_false(anyNA(fit$coefficients$beta))
})

test_that("dinarch_fit_bayes with a single covariate (n_beta == 1) returns the expected structure", {
  # Regression test: rep(0, 1) collapses to a bare R scalar, which Stan
  # rejects for a declared vector[1] parameter ("dims declared=(1); dims
  # found=()") unless explicitly wrapped with array(..., dim = 1) - this
  # would error on init before ever reaching the structure checks below.
  testthat::skip_if_not_installed("rstan")

  dat <- data.frame(y = rpois(100, 3), index = 1:100)
  # Deliberately few iterations here since this test only checks structure,
  # not estimation quality - low-ESS warnings are expected and harmless.
  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", n_lags = 1, formula = ~1,
    iter = 500, chains = 2, seed = 2
  ))

  expect_true(all(c("b", "phi", "beta") %in% names(fit$coefficients)))
  expect_length(fit$coefficients$b, 1)
  expect_identical(names(fit$coefficients$beta), "(Intercept)")
  expect_s4_class(fit$stanfit, "stanfit")
  expect_true(is.list(fit$posterior))
})

test_that("dinarch_fit_bayes pools multiple groups without crossing group boundaries", {
  testthat::skip_if_not_installed("rstan")

  # group = a series identity (e.g. country), index = order within that
  # group (e.g. year). Lags must never be taken across a group boundary,
  # so one observation per group is lost to the n_lags = 1 lag, regardless of
  # how many groups there are.
  n <- 60
  sim1 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 10)
  sim2 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 11)

  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n), group = "B")
  )

  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", group = "group", n_lags = 1,
    iter = 500, chains = 2, seed = 3
  ))

  expect_equal(fit$n_obs, 2 * (n - 1))
})
