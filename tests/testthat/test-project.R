test_that("dinarch_fit_ml() does not store `data` on the fit object", {
  n <- 50
  gdp <- rnorm(n)
  dat <- data.frame(y = rpois(n, 3), time = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 2, formula = ~gdp, vcov = FALSE)

  expect_null(fit$data)
})

test_that("dinarch_project() out-of-sample: default holds last covariate row fixed", {
  n <- 100
  gdp <- rnorm(n)
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, formula = ~gdp,
                           covariates = data.frame(gdp = gdp), beta = c(log(3), 0.3), seed = 1)
  dat <- data.frame(y = sim$y, time = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, formula = ~gdp)

  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 2, seed = 1)

  expect_s3_class(proj, "data.table")
  expect_identical(names(proj)[1:5], c("id", "time", "sim", "historical", "y"))
  expect_equal(nrow(proj), (n + 5) * 2)  # n historical rows + 5 projected rows, per sim

  future <- proj[!proj$historical]
  expect_equal(nrow(future), 5 * 2)
  expect_equal(future$time, rep((n + 1):(n + 5), times = 2))
  # gdp held fixed at the last observed value for every forecast period
  expect_true(all(future$gdp == tail(gdp, 1)))
})

test_that("dinarch_project() in-sample: unconditional replay covers n - n_lags periods per id, starting after the seed", {
  n <- 60
  dat <- data.frame(y = rpois(n, 4), time = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 2, vcov = FALSE)

  proj <- dinarch_project(fit, dat, horizon = NULL, nsim = 1, seed = 1)
  future <- proj[!proj$historical]

  expect_equal(nrow(future), n - 2)
  expect_equal(future$time, seq(3, n))
  expect_equal(nrow(proj), n + (n - 2))  # historical + simulated, nsim = 1
})

test_that("dinarch_project() is id-aware: each id forecasts from its own last time/lags", {
  n <- 40
  sim1 <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(3), seed = 1)
  sim2 <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(3), seed = 2)
  dat <- rbind(
    data.frame(y = sim1$y, time = seq_len(n), id = "A"),
    data.frame(y = sim2$y, time = seq_len(n) + 100, id = "B")  # different time range
  )
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", id = "id", n_lags = 1)

  proj <- dinarch_project(fit, dat, horizon = 3, nsim = 1, seed = 1, id = "id")
  future <- proj[!proj$historical]

  expect_setequal(unique(future$id), c("A", "B"))
  expect_equal(sort(future$time[future$id == "A"]), (n + 1):(n + 3))
  expect_equal(sort(future$time[future$id == "B"]), (n + 100 + 1):(n + 100 + 3))
})

test_that("dinarch_project() new_data overrides covariates only for matching (id, time)", {
  n <- 50
  gdp <- rnorm(n)
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, formula = ~gdp,
                           covariates = data.frame(gdp = gdp), beta = c(log(3), 0.3), seed = 1)
  dat <- data.frame(y = sim$y, time = seq_len(n), gdp = gdp)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, formula = ~gdp)

  new_gdp <- c(1.5, 2.5, 3.5)
  nd <- data.frame(time = (n + 1):(n + 3), gdp = new_gdp)
  proj <- dinarch_project(fit, dat, horizon = 3, nsim = 1, new_data = nd, seed = 1)
  future <- proj[!proj$historical]

  expect_equal(future$gdp, new_gdp)
})

test_that("dinarch_project() new_para overrides b/phi only for matching (id, time)", {
  n <- 50
  dat <- data.frame(y = rpois(n, 4), time = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)

  # phi = 10000 (near-deterministic) forces mu very close to y for periods
  # covered by new_para, and mu = exp(beta_0) + b1 * y_lag with b1 = 0
  # collapses to mu = exp(beta_0) exactly - a directly checkable value.
  np <- data.frame(time = (n + 1):(n + 3), b1 = 0, phi = 1e6)
  proj <- dinarch_project(fit, dat, horizon = 3, nsim = 1, new_para = np, seed = 1)
  future <- proj[!proj$historical]

  expected_mu <- exp(unname(fit$coefficients$beta["(Intercept)"]))
  expect_equal(future$y, rep(round(expected_mu), 3), tolerance = 1)
})

test_that("dinarch_project() nsim replicates vary stochastically", {
  n <- 50
  dat <- data.frame(y = rpois(n, 5), time = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)

  proj <- dinarch_project(fit, dat, horizon = 10, nsim = 5, seed = 1)
  future <- proj[!proj$historical]
  by_sim <- split(future$y, future$sim)
  expect_length(unique(by_sim), 5)  # exceedingly unlikely all 5 replicates coincide exactly
})

test_that("dinarch_project() validates id/time columns in new_data and new_para", {
  n <- 40
  sim1 <- dinarch_simulate(n = n, b = 0.2, phi = 6, beta = log(3), seed = 1)
  sim2 <- dinarch_simulate(n = n, b = 0.2, phi = 6, beta = log(3), seed = 2)
  dat <- rbind(
    data.frame(y = sim1$y, time = seq_len(n), id = "A"),
    data.frame(y = sim2$y, time = seq_len(n), id = "B")
  )
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", id = "id", n_lags = 1)

  expect_error(
    dinarch_project(fit, dat, horizon = 3, id = "id",
                     new_data = data.frame(id = "Z", time = (n + 1):(n + 3))),
    "not present in"
  )
  expect_error(
    dinarch_project(fit, dat, horizon = 3, id = "id",
                     new_para = data.frame(time = (n + 1):(n + 3), phi = 5)),
    "must have an `id` column"
  )
  expect_error(
    dinarch_project(fit, dat, horizon = 3, id = "id",
                     new_data = data.frame(id = "A")),
    "must have a `time` column"
  )
})

test_that("dinarch_project() warns when new_data/new_para don't match the simulated window", {
  n <- 40
  dat <- data.frame(y = rpois(n, 4), time = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)

  expect_warning(
    dinarch_project(fit, dat, horizon = 3, new_data = data.frame(time = 9999, gdp = 1)),
    "did not match"
  )
  expect_warning(
    dinarch_project(fit, dat, horizon = 3, new_para = data.frame(time = 9999, phi = 1)),
    "did not match"
  )
})

test_that("dinarch_project() Bayesian use_draw = 'random' varies parameters across replicates", {
  testthat::skip_if_not_installed("rstan")

  n <- 100
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 1)
  dat <- data.frame(y = sim$y, time = seq_len(n))
  fit <- suppressWarnings(dinarch_fit_bayes(
    dat, y = "y", time = "time", n_lags = 1, iter = 500, chains = 2, seed = 1
  ))

  proj <- dinarch_project(fit, dat, horizon = 3, nsim = 8, use_draw = "random", seed = 2)
  future <- proj[!proj$historical]
  expect_equal(nrow(future), 24)
  expect_false(anyNA(future$y))
})

test_that("dinarch_project() requires `data` and validates it like the fit functions do", {
  n <- 30
  dat <- data.frame(y = rpois(n, 3), time = seq_len(n))
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)

  expect_error(dinarch_project(fit, horizon = 3), "data")

  bad_dat <- dat
  bad_dat$time[5] <- bad_dat$time[4]  # duplicate time
  expect_error(dinarch_project(fit, bad_dat, horizon = 3), "strictly increasing")

  missing_col_dat <- data.frame(yy = dat$y, time = dat$time)
  expect_error(dinarch_project(fit, missing_col_dat, horizon = 3), "missing required column")
})

test_that("plot.dinarch_project() returns a ggplot for a single id, both types", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 10, seed = 2)

  expect_s3_class(proj, "dinarch_project")
  expect_s3_class(plot(proj), "ggplot")
  expect_s3_class(plot(proj, type = "spaghetti"), "ggplot")
})

test_that("plot.dinarch_project() facets by id only when there is more than one", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 40
  sim1 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 10)
  sim2 <- dinarch_simulate(n = n, b = 0.3, phi = 6, beta = log(3), seed = 11)
  dat <- rbind(
    data.frame(y = sim1$y, time = seq_len(n), id = "A"),
    data.frame(y = sim2$y, time = seq_len(n), id = "B")
  )
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", id = "id", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 10, seed = 4, id = "id")

  p_multi <- plot(proj)
  expect_true(length(p_multi$facet$params$facets) > 0)

  dat1 <- data.frame(y = sim1$y, time = seq_len(n))
  fit1 <- dinarch_fit_ml(dat1, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj1 <- dinarch_project(fit1, dat1, horizon = 5, nsim = 10, seed = 4)
  p_single <- plot(proj1)
  expect_length(p_single$facet$params$facets, 0)
})

test_that("plot.dinarch_project() works for an in-sample replay (no future times)", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 40
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, beta = log(3), seed = 5)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = NULL, nsim = 10, seed = 5)

  expect_s3_class(plot(proj), "ggplot")
})

test_that("plot.dinarch_project() can be extended by the caller", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 40
  sim <- dinarch_simulate(n = n, b = 0.2, phi = 8, beta = log(3), seed = 5)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 10, seed = 5)

  p <- plot(proj) + ggplot2::ggtitle("Forecast")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Forecast")
})

test_that("plot.dinarch_project() step = TRUE returns a ggplot for both types", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 10, seed = 2)

  expect_s3_class(plot(proj, step = TRUE), "ggplot")
  expect_s3_class(plot(proj, type = "spaghetti", step = TRUE), "ggplot")
})

test_that("plot.dinarch_project() connects the projection to the last historical point", {
  testthat::skip_if_not_installed("ggplot2")

  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 10, seed = 2)
  hist_last <- utils::tail(dat$y, 1)

  ribbon_layer <- ggplot2::ggplot_build(plot(proj))$data[[2]]
  first_pt <- ribbon_layer[which.min(ribbon_layer$x), ]
  expect_equal(first_pt$ymin, hist_last)
  expect_equal(first_pt$ymax, hist_last)

  spaghetti_layer <- ggplot2::ggplot_build(plot(proj, type = "spaghetti"))$data[[2]]
  first_pts <- spaghetti_layer[spaghetti_layer$x == min(spaghetti_layer$x), ]
  expect_true(all(first_pts$y == hist_last))
})

test_that("summary.dinarch_project() returns one mean/interval row per (id, time)", {
  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 200, seed = 2)

  s <- summary(proj)
  expect_s3_class(s, "summary.dinarch_project")
  expect_identical(names(s$table), c("id", "time", "mean", "lower", "upper"))
  expect_equal(nrow(s$table), 5)
  expect_equal(s$table$time, (n + 1):(n + 5))
  expect_true(all(s$table$lower <= s$table$mean & s$table$mean <= s$table$upper))
  expect_output(print(s), "prediction interval")
})

test_that("summary.dinarch_project()'s level controls the interval width", {
  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  proj <- dinarch_project(fit, dat, horizon = 5, nsim = 500, seed = 2)

  s_narrow <- summary(proj, level = 0.5)
  s_wide <- summary(proj, level = 0.95)
  expect_true(all((s_wide$table$upper - s_wide$table$lower) >=
                     (s_narrow$table$upper - s_narrow$table$lower)))
})

test_that("summary.dinarch_project() works on predict.dinarch_fit()'s output too", {
  n <- 60
  sim <- dinarch_simulate(n = n, b = 0.3, phi = 8, beta = log(4), seed = 2)
  dat <- data.frame(y = sim$y, time = sim$time)
  fit <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)

  s <- summary(predict(fit, dat, horizon = 5, nsim = 50, seed = 1))
  expect_s3_class(s, "summary.dinarch_project")
  expect_equal(nrow(s$table), 5)
})
