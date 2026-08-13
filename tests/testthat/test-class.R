# Shared fixtures ------------------------------------------------------------
# Fit once per file, reused across tests below, to avoid repeated (slow)
# optimization/MCMC calls.

.test_ml_fit <- {
  sim <- dinarch_simulate(n = 200, b = 0.3, phi = 8, beta = log(4), seed = 1)
  dinarch_fit_ml(
    data.frame(y = sim$y, index = seq_len(200)),
    y = "y", index = "index", n_lags = 1
  )
}

if (requireNamespace("rstan", quietly = TRUE)) {
  sim_bayes <- dinarch_simulate(n = 150, b = 0.3, phi = 6, beta = log(3), seed = 77)
  .test_bayes_fit <- suppressWarnings(dinarch_fit_bayes(
    data.frame(y = sim_bayes$y, index = seq_len(150)),
    y = "y", index = "index", n_lags = 1, iter = 800, chains = 2, seed = 8
  ))
}

# summary.dinarch_fit ---------------------------------------------------------

test_that("summary.dinarch_fit works for an ML fit", {
  s <- summary(.test_ml_fit)
  expect_s3_class(s, "summary.dinarch_fit")
  expect_equal(s$coefficients$term, c("b[1]", "phi", "(Intercept)"))
  expect_true(is.finite(s$aic) && is.finite(s$bic))
  expect_output(print(s), "DINARCH fit")
})

test_that("summary.dinarch_fit includes credible intervals for a Bayes fit", {
  skip_if_not_installed("rstan")
  s <- summary(.test_bayes_fit)
  expect_true(all(c("lower", "upper") %in% names(s$coefficients)))
  expect_true(all(s$coefficients$lower <= s$coefficients$estimate))
  expect_true(all(s$coefficients$estimate <= s$coefficients$upper))
})

# simulate.dinarch_fit ---------------------------------------------------------

test_that("simulate.dinarch_fit matches a manual dinarch_simulate() call with the same parameters", {
  sim1 <- simulate(.test_ml_fit, nsim = 20, seed = 99)
  sim2 <- dinarch_simulate(
    n = 20,
    b = .test_ml_fit$coefficients$b,
    phi = .test_ml_fit$coefficients$phi,
    beta = .test_ml_fit$coefficients$beta,
    seed = 99
  )
  expect_identical(sim1$y, sim2$y)
})

test_that("simulate.dinarch_fit warns (then errors informatively) if the fit's formula referenced covariates not supplied for simulation", {
  n <- 200
  gdp <- rnorm(n)
  sim_gdp <- dinarch_simulate(
    n = n, b = 0.3, phi = 6,
    formula = ~gdp, covariates = data.frame(gdp = gdp), beta = c(log(3), 0.2), seed = 3
  )
  fit_gdp <- dinarch_fit_ml(
    data.frame(y = sim_gdp$y, index = seq_len(n), gdp = gdp),
    y = "y", index = "index", formula = ~gdp, n_lags = 1
  )
  # covariates isn't supplied, so this both warns up front (formula
  # references gdp) and then errors from model.matrix() itself (gdp isn't
  # resolvable - formula's environment is stripped to baseenv() precisely
  # so this fails loudly rather than silently reusing gdp from wherever
  # the fit's formula happened to be written).
  expect_warning(
    expect_error(simulate(fit_gdp, nsim = 10, seed = 1), "gdp"),
    "gdp"
  )
})

test_that("simulate.dinarch_fit's use_draw errors sensibly", {
  expect_error(simulate(.test_ml_fit, nsim = 5, use_draw = 1), "Bayesian")

  skip_if_not_installed("rstan")
  n_draws <- length(.test_bayes_fit$posterior$phi)
  expect_error(simulate(.test_bayes_fit, nsim = 5, use_draw = n_draws + 1), "use_draw")
})

test_that("simulate.dinarch_fit's use_draw runs using a specific posterior draw", {
  skip_if_not_installed("rstan")
  sim_mean <- simulate(.test_bayes_fit, nsim = 30, seed = 5)
  sim_draw <- simulate(.test_bayes_fit, nsim = 30, seed = 5, use_draw = 1)
  expect_equal(nrow(sim_mean), 30)
  expect_equal(nrow(sim_draw), 30)
})
