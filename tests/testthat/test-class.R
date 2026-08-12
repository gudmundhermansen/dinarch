# Shared fixtures ------------------------------------------------------------
# Fit once per file, reused across tests below, to avoid repeated (slow)
# optimization/MCMC calls.

.test_ml_fit <- {
  sim <- dinarch_simulate(n = 200, a = log(4), b = 0.3, phi = 8, seed = 1)
  dinarch_fit_ml(
    data.frame(y = sim$y, index = seq_len(200)),
    y = "y", index = "index", k = 1
  )
}

if (requireNamespace("rstan", quietly = TRUE)) {
  sim_bayes <- dinarch_simulate(n = 150, a = log(3), b = 0.3, phi = 6, seed = 77)
  .test_bayes_fit <- suppressWarnings(dinarch_fit_bayes(
    data.frame(y = sim_bayes$y, index = seq_len(150)),
    y = "y", index = "index", k = 1, iter = 800, chains = 2, seed = 8
  ))
}

# summary.dinarch_fit ---------------------------------------------------------

test_that("summary.dinarch_fit works for an ML fit", {
  s <- summary(.test_ml_fit)
  expect_s3_class(s, "summary.dinarch_fit")
  expect_equal(s$coefficients$term, c("a", "b[1]", "phi"))
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
    a = .test_ml_fit$coefficients$a,
    b = .test_ml_fit$coefficients$b,
    phi = .test_ml_fit$coefficients$phi,
    seed = 99
  )
  expect_identical(sim1$y, sim2$y)
})

test_that("simulate.dinarch_fit warns if population was used in fitting but not supplied for simulation", {
  sim_pop <- dinarch_simulate(
    n = 100, a = log(0.01), b = 0.3, phi = 6,
    population = rep(1000, 100), seed = 3
  )
  fit_pop <- dinarch_fit_ml(
    data.frame(y = sim_pop$y, index = seq_len(100), population = rep(1000, 100)),
    y = "y", index = "index", population = "population", k = 1
  )
  expect_warning(simulate(fit_pop, nsim = 10, seed = 1), "population")
})

test_that("simulate.dinarch_fit's use_draw errors sensibly", {
  expect_error(simulate(.test_ml_fit, nsim = 5, use_draw = 1), "Bayesian")

  skip_if_not_installed("rstan")
  n_draws <- length(.test_bayes_fit$posterior$a)
  expect_error(simulate(.test_bayes_fit, nsim = 5, use_draw = n_draws + 1), "use_draw")
})

test_that("simulate.dinarch_fit's use_draw runs using a specific posterior draw", {
  skip_if_not_installed("rstan")
  sim_mean <- simulate(.test_bayes_fit, nsim = 30, seed = 5)
  sim_draw <- simulate(.test_bayes_fit, nsim = 30, seed = 5, use_draw = 1)
  expect_equal(nrow(sim_mean), 30)
  expect_equal(nrow(sim_draw), 30)
})
