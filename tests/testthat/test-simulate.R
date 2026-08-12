test_that("dinarch_simulate returns the right shape and columns", {
  sim <- dinarch_simulate(n = 20, a = log(2), b = 0.2, phi = 5, seed = 1)
  expect_s3_class(sim, "data.table")
  expect_equal(nrow(sim), 20)
  expect_equal(names(sim), c("index", "y", "population"))
  expect_equal(sim$population, rep(1, 20))
})

test_that("results are reproducible with the same seed", {
  sim1 <- dinarch_simulate(n = 15, a = 0, b = 0.3, phi = 10, seed = 42)
  sim2 <- dinarch_simulate(n = 15, a = 0, b = 0.3, phi = 10, seed = 42)
  expect_identical(sim1$y, sim2$y)
})

test_that("mean matches exp(a) for a pure-intercept model", {
  # E[y] = mu = exp(a) always, for any phi. With phi = 1000 the variance
  # is mu*(1/phi + 1) ~ mu, so over n = 2000 draws the sample mean should
  # land well within 1 of the true mean of 5.
  sim <- dinarch_simulate(n = 2000, a = log(5), b = 0, phi = 1000, seed = 123)
  expect_lt(abs(mean(sim$y) - 5), 1)
})

test_that("a fully time-varying n x k matrix b is accepted", {
  n <- 10
  b_mat <- matrix(c(rep(0.3, n), rep(0.1, n)), nrow = n, ncol = 2)
  sim <- dinarch_simulate(
    n = n, a = 0, b = b_mat, phi = 20, y_init = c(1, 2), seed = 5
  )
  expect_equal(nrow(sim), n)
})

test_that("dynamic_population never lets population increase", {
  sim <- dinarch_simulate(
    n = 15, a = log(3), b = 0.2, phi = 5,
    population = 100, dynamic_population = TRUE, seed = 7
  )
  expect_true(all(diff(sim$population) <= 0))
})

test_that("population fully depleted does not error, and can reach exactly 0", {
  # Small population plus a high death rate forces depletion toward 0
  # well before n periods are up. This should never error - it's a
  # regression test for the mu = 0 boundary case (population depleted,
  # no lag momentum) - and population should be able to reach exactly 0.
  sim <- dinarch_simulate(
    n = 30, a = log(50), b = 0.1, phi = 5,
    population = 5, dynamic_population = TRUE, seed = 3
  )
  expect_true(any(sim$population == 0))
})

test_that("threshold zeroes out y once recent deaths exceed the limit", {
  sim <- dinarch_simulate(
    n = 6, a = 5, b = 0, phi = 50, population = 10,
    threshold = list(window = 2, fraction = 0.001, pop_unit = 1),
    seed = 11
  )
  expect_true(any(sim$y[-1] == 0))
})

test_that("covariates require beta with matching dimensions", {
  expect_error(
    dinarch_simulate(
      n = 5, a = 0, b = 0.1, phi = 5, covariates = matrix(1:5, ncol = 1)
    ),
    "beta"
  )
  expect_error(
    dinarch_simulate(
      n = 5, a = 0, b = 0.1, phi = 5,
      covariates = matrix(1:10, ncol = 2), beta = 0.1
    ),
    "beta"
  )
})

test_that("dynamic_population and threshold require population", {
  expect_error(
    dinarch_simulate(n = 5, a = 0, b = 0.1, phi = 5, dynamic_population = TRUE),
    "population"
  )
  expect_error(
    dinarch_simulate(
      n = 5, a = 0, b = 0.1, phi = 5,
      threshold = list(window = 2, fraction = 0.01, pop_unit = 1)
    ),
    "population"
  )
})

test_that("b outside [0, 1) is rejected", {
  expect_error(dinarch_simulate(n = 5, a = 0, b = 1.5, phi = 5), "\\[0, 1\\)")
})

test_that("y_init length must match k", {
  expect_error(
    dinarch_simulate(n = 5, a = 0, b = c(0.2, 0.1), phi = 5, y_init = 1),
    "y_init"
  )
})
