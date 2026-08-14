test_that("dinarch_fit_ml() stores $data with full per-group history, including pre-lag rows and raw covariates", {
  n <- 50
  gdp <- rnorm(n)
  dat <- data.frame(y = rpois(n, 3), index = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 2, formula = ~gdp, vcov = FALSE)

  expect_s3_class(fit$data, "data.table")
  expect_equal(nrow(fit$data), n)  # all rows, including the first n_lags
  expect_identical(names(fit$data), c("group", "index", "y", "gdp"))
  expect_equal(fit$data$gdp, gdp)
})

test_that("dinarch_project() out-of-sample: default holds last covariate row fixed", {
  n <- 100
  gdp <- rnorm(n)
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, formula = ~gdp,
                           covariates = data.frame(gdp = gdp), beta = c(log(3), 0.3), seed = 1)
  dat <- data.frame(y = sim$y, index = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, formula = ~gdp)

  proj <- dinarch_project(fit, horizon = 5, nsim = 2, seed = 1)

  expect_s3_class(proj, "data.table")
  expect_identical(names(proj)[1:4], c("group", "index", "sim", "y"))
  expect_equal(nrow(proj), 5 * 2)
  expect_equal(proj$index, rep(rep((n + 1):(n + 5), times = 2)))
  # gdp held fixed at the last observed value for every forecast period
  expect_true(all(proj$gdp == tail(gdp, 1)))
})

test_that("dinarch_project() in-sample: unconditional replay covers n - n_lags periods per group, starting after the seed", {
  n <- 60
  dat <- data.frame(y = rpois(n, 4), index = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 2, vcov = FALSE)

  proj <- dinarch_project(fit, horizon = NULL, nsim = 1, seed = 1)

  expect_equal(nrow(proj), n - 2)
  expect_equal(proj$index, seq(3, n))
})

test_that("dinarch_project() is group-aware: each group forecasts from its own last index/lags", {
  n <- 40
  sim1 <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(3), seed = 1)
  sim2 <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(3), seed = 2)
  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n) + 100, group = "B")  # different index range
  )
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", group = "group", n_lags = 1)

  proj <- dinarch_project(fit, horizon = 3, nsim = 1, seed = 1)

  expect_setequal(unique(proj$group), c("A", "B"))
  expect_equal(sort(proj$index[proj$group == "A"]), (n + 1):(n + 3))
  expect_equal(sort(proj$index[proj$group == "B"]), (n + 100 + 1):(n + 100 + 3))
})

test_that("dinarch_project() new_data overrides covariates only for matching (group, index)", {
  n <- 50
  gdp <- rnorm(n)
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, formula = ~gdp,
                           covariates = data.frame(gdp = gdp), beta = c(log(3), 0.3), seed = 1)
  dat <- data.frame(y = sim$y, index = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, formula = ~gdp)

  new_gdp <- c(1.5, 2.5, 3.5)
  nd <- data.frame(index = (n + 1):(n + 3), gdp = new_gdp)
  proj <- dinarch_project(fit, horizon = 3, nsim = 1, new_data = nd, seed = 1)

  expect_equal(proj$gdp, new_gdp)
})

test_that("dinarch_project() new_para overrides b/phi only for matching (group, index)", {
  n <- 50
  dat <- data.frame(y = rpois(n, 4), index = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, vcov = FALSE)

  # phi = 10000 (near-deterministic) forces mu very close to y for periods
  # covered by new_para, and mu = exp(beta_0) + b1 * y_lag with b1 = 0
  # collapses to mu = exp(beta_0) exactly - a directly checkable value.
  np <- data.frame(index = (n + 1):(n + 3), b1 = 0, phi = 1e6)
  proj <- dinarch_project(fit, horizon = 3, nsim = 1, new_para = np, seed = 1)

  expected_mu <- exp(unname(fit$coefficients$beta["(Intercept)"]))
  expect_equal(proj$y, rep(round(expected_mu), 3), tolerance = 1)
})

test_that("dinarch_project() nsim replicates vary stochastically", {
  n <- 50
  dat <- data.frame(y = rpois(n, 5), index = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, vcov = FALSE)

  proj <- dinarch_project(fit, horizon = 10, nsim = 5, seed = 1)
  by_sim <- split(proj$y, proj$sim)
  expect_length(unique(by_sim), 5)  # exceedingly unlikely all 5 replicates coincide exactly
})

test_that("dinarch_project() validates group/index columns in new_data and new_para", {
  n <- 40
  sim1 <- dinarch_simulate(n = n, b = 0.2, phi = 6, beta = log(3), seed = 1)
  sim2 <- dinarch_simulate(n = n, b = 0.2, phi = 6, beta = log(3), seed = 2)
  dat <- rbind(
    data.frame(y = sim1$y, index = seq_len(n), group = "A"),
    data.frame(y = sim2$y, index = seq_len(n), group = "B")
  )
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", group = "group", n_lags = 1)

  expect_error(
    dinarch_project(fit, horizon = 3, new_data = data.frame(group = "Z", index = (n + 1):(n + 3))),
    "not present in the fitted data"
  )
  expect_error(
    dinarch_project(fit, horizon = 3, new_para = data.frame(index = (n + 1):(n + 3), phi = 5)),
    "must have a `group` column"
  )
  expect_error(
    dinarch_project(fit, horizon = 3, new_data = data.frame(group = "A")),
    "must have an `index` column"
  )
})

test_that("dinarch_project() warns when new_data/new_para don't match the simulated window", {
  n <- 40
  dat <- data.frame(y = rpois(n, 4), index = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, vcov = FALSE)

  expect_warning(
    dinarch_project(fit, horizon = 3, new_data = data.frame(index = 9999, gdp = 1)),
    "did not match"
  )
  expect_warning(
    dinarch_project(fit, horizon = 3, new_para = data.frame(index = 9999, phi = 1)),
    "did not match"
  )
})

test_that("dinarch_project() Bayesian use_draw = 'random' varies parameters across replicates", {
  testthat::skip_if_not_installed("rstan")

  n <- 100
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 1)
  dat <- data.frame(y = sim$y, index = seq_len(n))
  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", index = "index", n_lags = 1, iter = 500, chains = 2, seed = 1
  ))

  proj <- dinarch_project(fit, horizon = 3, nsim = 8, use_draw = "random", seed = 2)
  expect_equal(nrow(proj), 24)
  expect_false(anyNA(proj$y))
})

test_that("dinarch_project() errors informatively on a fit without stored $data", {
  n <- 30
  dat <- data.frame(y = rpois(n, 3), index = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1, vcov = FALSE)
  fit$data <- NULL

  expect_error(dinarch_project(fit, horizon = 3), "no stored")
})
