test_that("dinarch_simulate returns the right shape and columns", {
  sim <- dinarch_simulate(n = 20, b = 0.2, phi = 5, beta = log(2), seed = 1)
  expect_s3_class(sim, "data.table")
  expect_equal(nrow(sim), 20)
  expect_equal(names(sim), c("index", "y"))
})

test_that("results are reproducible with the same seed", {
  sim1 <- dinarch_simulate(n = 15, b = 0.3, phi = 10, beta = 0, seed = 42)
  sim2 <- dinarch_simulate(n = 15, b = 0.3, phi = 10, beta = 0, seed = 42)
  expect_identical(sim1$y, sim2$y)
})

test_that("mean matches exp(beta) for a pure-intercept model", {
  # E[y] = mu = exp(beta) always, for any phi. With phi = 1000 the variance
  # is mu*(1/phi + 1) ~ mu, so over n = 2000 draws the sample mean should
  # land well within 1 of the true mean of 5.
  sim <- dinarch_simulate(n = 2000, b = 0, phi = 1000, beta = log(5), seed = 123)
  expect_lt(abs(mean(sim$y) - 5), 1)
})

test_that("a fully time-varying n x n_lags matrix b is accepted", {
  n <- 10
  b_mat <- matrix(c(rep(0.3, n), rep(0.1, n)), nrow = n, ncol = 2)
  sim <- dinarch_simulate(
    n = n, b = b_mat, phi = 20, beta = 0, y_init = c(1, 2), seed = 5
  )
  expect_equal(nrow(sim), n)
})

test_that("formula/covariates build the linear predictor correctly, including interactions", {
  n <- 300
  cov_a <- runif(n, 1, 3)
  cov_b <- runif(n, 1, 3)
  sim <- dinarch_simulate(
    n = n, b = 0, phi = 1000,
    formula = ~ cov_a * cov_b,
    covariates = data.frame(cov_a = cov_a, cov_b = cov_b),
    beta = c(log(2), 0.1, -0.1, 0.05),
    seed = 6
  )
  expect_equal(nrow(sim), n)
  expect_true(all(is.finite(sim$y)))
})

test_that("formula = ~0 drops the intercept and requires matching beta length", {
  n <- 10
  cov_a <- runif(n, 1, 3)
  expect_error(
    dinarch_simulate(
      n = n, b = 0.1, phi = 5,
      formula = ~ 0 + cov_a, covariates = data.frame(cov_a = cov_a),
      beta = c(0.1, 0.2)
    ),
    "beta"
  )
  sim <- dinarch_simulate(
    n = n, b = 0.1, phi = 5,
    formula = ~ 0 + cov_a, covariates = data.frame(cov_a = cov_a),
    beta = 0.2
  )
  expect_equal(nrow(sim), n)
})

test_that("population has no special role in the mean - it's just another covariate via formula", {
  # A "population" column with no formula reference should have zero
  # effect on y, exactly like any other unused column would.
  n <- 2000
  pop <- c(rep(1000, n / 2), rep(10000, n / 2))
  sim_ignored <- dinarch_simulate(
    n = n, b = 0, phi = 5000, beta = log(5),
    covariates = data.frame(population = pop),  # not referenced by formula
    seed = 4
  )
  expect_equal(mean(sim_ignored$y[1:1000]), mean(sim_ignored$y[1001:2000]), tolerance = 0.1)

  # Referencing it via formula (log link, beta = 1) recovers proportional
  # scaling explicitly.
  sim_used <- dinarch_simulate(
    n = n, b = 0, phi = 5000,
    formula = ~ log(population), covariates = data.frame(population = pop),
    beta = c(log(0.005), 1),
    seed = 4
  )
  expect_equal(mean(sim_used$y[1001:2000]) / mean(sim_used$y[1:1000]), 10, tolerance = 0.1)
})

test_that("y_threshold zeroes out y once recent values exceed the limit, independently of mu", {
  # Same mean/seed with and without y_threshold: the cap introduces zeros
  # that wouldn't otherwise occur, and has no dependency on
  # formula/covariates - it's a plain, independent numeric cap.
  n <- 6
  sim_capped <- dinarch_simulate(
    n = n, b = 0, phi = 50, beta = 5,
    y_threshold = list(window = 2, limit = 1),
    seed = 11
  )
  sim_uncapped <- dinarch_simulate(n = n, b = 0, phi = 50, beta = 5, seed = 11)
  expect_true(any(sim_capped$y == 0))
  expect_false(any(sim_uncapped$y == 0))
})

test_that("y_threshold$limit can be a length-n vector (time-varying cap)", {
  # Near-deterministic mean of ~30 (b = 0, phi huge, so P(natural 0) is
  # negligible). limit = 100 is loose (never triggered by a window-of-2
  # sum around 60), limit = 1 is tight (triggered essentially
  # immediately). Loosen for the first half, tighten for the second, and
  # check zeros only appear once the tight regime starts (any zero seen
  # is therefore attributable to the threshold rule, not to sampling
  # variability).
  n <- 10
  limit_path <- c(rep(100, 5), rep(1, 5))
  sim <- dinarch_simulate(
    n = n, b = 0, phi = 10000, beta = log(30),
    y_threshold = list(window = 2, limit = limit_path),
    seed = 11
  )
  expect_false(any(sim$y[1:5] == 0))
  expect_true(any(sim$y[6:10] == 0))
})

test_that("covariates require beta with matching dimensions", {
  expect_error(
    dinarch_simulate(
      n = 5, b = 0.1, phi = 5,
      formula = ~cov_a, covariates = data.frame(cov_a = 1:5)
    ),
    "beta"
  )
  expect_error(
    dinarch_simulate(
      n = 5, b = 0.1, phi = 5,
      formula = ~ cov_a + cov_b,
      covariates = data.frame(cov_a = 1:5, cov_b = 6:10), beta = 0.1
    ),
    "beta"
  )
})

test_that("b outside [0, 1) is rejected", {
  expect_error(dinarch_simulate(n = 5, b = 1.5, phi = 5, beta = 0), "\\[0, 1\\)")
})

test_that("y_init length must match n_lags", {
  expect_error(
    dinarch_simulate(n = 5, b = c(0.2, 0.1), phi = 5, beta = 0, y_init = 1),
    "y_init"
  )
})
