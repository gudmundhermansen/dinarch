test_that("dinarch_transition_path returns the right length and endpoints", {
  path <- dinarch_transition_path(from = 0, to = 10, n = 25)
  expect_length(path, 25)
  expect_lt(path[1], 1)       # close to `from`, given the default front-loaded range
  expect_gt(path[25], 9.9)    # close to `to`
})

test_that("dinarch_transition_path with shape = 'linear' interpolates exactly", {
  path <- dinarch_transition_path(from = 0, to = 10, n = 5, shape = "linear")
  expect_equal(path, c(0, 2.5, 5, 7.5, 10))
})

test_that("a symmetric sigmoid range hits the midpoint exactly at x = 0", {
  path <- dinarch_transition_path(from = 0, to = 1, n = 3, range = c(-6, 6))
  expect_equal(path[2], 0.5)
  expect_lt(path[1], 0.01)
  expect_gt(path[3], 0.99)
})

test_that("dinarch_transition_path handles a decreasing path (to < from)", {
  path <- dinarch_transition_path(from = 5, to = 1, n = 10)
  expect_true(all(diff(path) <= 0))
  expect_equal(path, dinarch_transition_path(from = 5, to = 1, n = 10))
})

test_that("a transition path plugs directly into dinarch_simulate()'s time-varying phi", {
  phi_path <- dinarch_transition_path(from = 2, to = 20, n = 30)
  sim <- dinarch_simulate(n = 30, b = 0.2, phi = phi_path, beta = log(3), seed = 1)
  expect_equal(nrow(sim), 30)
})

# dinarch_transition ------------------------------------------------

test_that("dinarch_transition builds a new_para table ending exactly at fit_to's parameters", {
  sim_a <- dinarch_simulate(n = 100, b = 0.2, phi = 8, beta = log(3), seed = 1)
  sim_b <- dinarch_simulate(n = 100, b = 0.6, phi = 8, beta = log(10), seed = 2)
  dat_a <- data.frame(y = sim_a$y, time = seq_len(100))
  fit_a <- dinarch_fit_ml(dat_a, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  fit_b <- dinarch_fit_ml(
    data.frame(y = sim_b$y, time = seq_len(100)), y = "y", time = "time", n_lags = 1, vcov = FALSE
  )

  # shape = "linear" reaches its endpoint exactly (unlike the default
  # sigmoid, which only approaches it asymptotically) - this keeps the
  # test independent of how extreme a fitted phi happens to be, since phi
  # is often poorly identified and can land on very large estimates.
  np <- dinarch_transition(fit_a, fit_b, dat_a, horizon = 20, shape = "linear")

  expect_s3_class(np, "data.table")
  expect_identical(names(np), c("time", "b1", "phi", "(Intercept)"))
  expect_equal(nrow(np), 20)
  expect_equal(np$time, 101:120)
  expect_equal(np$b1[20], unname(fit_b$coefficients$b[1]))
  expect_equal(np$phi[20], fit_b$coefficients$phi)
  expect_equal(np[["(Intercept)"]][20], unname(fit_b$coefficients$beta["(Intercept)"]))
  expect_equal(np$b1[1], unname(fit_a$coefficients$b[1]))
})

test_that("dinarch_transition's output plugs directly into dinarch_project()'s new_para", {
  sim_a <- dinarch_simulate(n = 100, b = 0.2, phi = 8, beta = log(3), seed = 1)
  sim_b <- dinarch_simulate(n = 100, b = 0.6, phi = 8, beta = log(10), seed = 2)
  dat_a <- data.frame(y = sim_a$y, time = seq_len(100))
  fit_a <- dinarch_fit_ml(dat_a, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  fit_b <- dinarch_fit_ml(
    data.frame(y = sim_b$y, time = seq_len(100)), y = "y", time = "time", n_lags = 1, vcov = FALSE
  )

  np <- dinarch_transition(fit_a, fit_b, dat_a, horizon = 20)
  proj <- dinarch_project(fit_a, dat_a, horizon = 20, nsim = 10, new_para = np, seed = 1)
  future <- proj[proj$pred]
  expect_equal(nrow(future), 20 * 10)
})

test_that("dinarch_transition indexes each id from its own last observed time", {
  sim1 <- dinarch_simulate(n = 40, b = 0.3, phi = 6, beta = log(3), seed = 10)
  sim2 <- dinarch_simulate(n = 30, b = 0.3, phi = 6, beta = log(3), seed = 11)
  dat_g <- rbind(
    data.frame(y = sim1$y, time = seq_len(40), id = "A"),
    data.frame(y = sim2$y, time = seq_len(30), id = "B")
  )
  fit_g <- dinarch_fit_ml(dat_g, y = "y", time = "time", id = "id", n_lags = 1, vcov = FALSE)

  sim_b <- dinarch_simulate(n = 100, b = 0.6, phi = 8, beta = log(5), seed = 20)
  fit_b <- dinarch_fit_ml(
    data.frame(y = sim_b$y, time = seq_len(100)), y = "y", time = "time", n_lags = 1, vcov = FALSE
  )

  np <- dinarch_transition(fit_g, fit_b, dat_g, id = "id", horizon = 5)
  expect_identical(names(np)[1:2], c("id", "time"))
  expect_equal(np$time[np$id == "A"], 41:45)
  expect_equal(np$time[np$id == "B"], 31:35)
})

test_that("dinarch_transition errors on mismatched n_lags or beta names", {
  sim <- dinarch_simulate(n = 100, b = 0.2, phi = 8, beta = log(3), seed = 1)
  dat <- data.frame(y = sim$y, time = seq_len(100))
  fit_1lag <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 1, vcov = FALSE)
  fit_2lag <- dinarch_fit_ml(dat, y = "y", time = "time", n_lags = 2, vcov = FALSE)
  expect_error(dinarch_transition(fit_1lag, fit_2lag, dat, horizon = 5), "n_lags")

  gdp <- rnorm(100)
  sim_gdp <- dinarch_simulate(
    n = 100, b = 0.2, phi = 8, formula = ~gdp, covariates = data.frame(gdp = gdp),
    beta = c(log(3), 0.4), seed = 2
  )
  fit_gdp <- dinarch_fit_ml(
    data.frame(y = sim_gdp$y, time = seq_len(100), gdp = gdp),
    y = "y", time = "time", formula = ~gdp, n_lags = 1, vcov = FALSE
  )
  expect_error(dinarch_transition(fit_1lag, fit_gdp, dat, horizon = 5), "beta")
})

.fake_fit <- function(b, phi, beta, formula = ~1, method = "ml", posterior = NULL) {
  structure(
    list(method = method, coefficients = list(b = b, phi = phi, beta = beta),
         n_lags = length(b), formula = formula, posterior = posterior),
    class = "dinarch_fit"
  )
}

test_that("dinarch_transition holds the long-run mean constant between equal-mean regimes", {
  fit_a <- .fake_fit(b = 0.8, phi = 10, beta = c("(Intercept)" = 0))
  fit_b <- .fake_fit(b = 0, phi = 10, beta = c("(Intercept)" = log(5)))
  set.seed(1)
  dat <- data.frame(y = rpois(30, 5), time = 1:30)
  np <- dinarch_transition(fit_a, fit_b, dat, horizon = 11, shape = "linear")
  expect_equal(exp(np[["(Intercept)"]]) / (1 - np$b1), rep(5, 11))
  expect_equal(np$b1, seq(0.8, 0, length.out = 11))
})

test_that("dinarch_transition reproduces endpoints and interpolates sum(b) and 1/phi linearly", {
  fit_a <- .fake_fit(b = c(0.5, 0.2), phi = 2, beta = c("(Intercept)" = 0.3))
  fit_b <- .fake_fit(b = c(0.1, 0.3), phi = 200, beta = c("(Intercept)" = 1.2))
  dat <- data.frame(y = rpois(30, 5), time = 1:30)
  np <- dinarch_transition(fit_a, fit_b, dat, horizon = 9, shape = "linear")
  cols <- c("b1", "b2", "phi", "(Intercept)")
  expect_equal(unlist(np[1, cols, with = FALSE], use.names = FALSE), c(0.5, 0.2, 2, 0.3))
  expect_equal(unlist(np[9, cols, with = FALSE], use.names = FALSE), c(0.1, 0.3, 200, 1.2))
  expect_equal(diff(1 / np$phi), rep(diff(1 / np$phi)[1], 8))
  expect_equal(np$b1 + np$b2, seq(0.7, 0.4, length.out = 9))
})

test_that("dinarch_transition interpolates the long-run mean at each id's last covariate row", {
  fit_a <- .fake_fit(b = 0.4, phi = 5, beta = c("(Intercept)" = 0.5, x = 0.3), formula = ~x)
  fit_b <- .fake_fit(b = 0.1, phi = 5, beta = c("(Intercept)" = 1, x = -0.2), formula = ~x)
  dat <- rbind(
    data.frame(y = rpois(20, 3), time = 1:20, id = "A", x = c(rnorm(19), 0)),
    data.frame(y = rpois(15, 3), time = 1:15, id = "B", x = c(rnorm(14), 10))
  )
  h <- 7
  np <- dinarch_transition(fit_a, fit_b, dat, horizon = h, shape = "linear", id = "id")
  for (g in c("A", "B")) {
    x_g <- if (g == "A") 0 else 10
    r <- np[np$id == g]
    lm_t <- r[["(Intercept)"]] + x_g * r$x - log(1 - r$b1)
    lm_A <- 0.5 + x_g * 0.3 - log(1 - 0.4)
    lm_B <- 1 + x_g * -0.2 - log(1 - 0.1)
    expect_equal(lm_t, seq(lm_A, lm_B, length.out = h))
    expect_equal(r$x, seq(0.3, -0.2, length.out = h))
  }
})

test_that("dinarch_transition errors without an intercept or with nsim on ML fits", {
  dat <- data.frame(y = rpois(20, 3), time = 1:20, x = rnorm(20))
  f0 <- .fake_fit(b = 0.3, phi = 5, beta = c(x = 0.2), formula = ~0 + x)
  expect_error(dinarch_transition(f0, f0, dat, horizon = 3), "Intercept")

  f1 <- .fake_fit(b = 0.3, phi = 5, beta = c("(Intercept)" = 1))
  expect_error(dinarch_transition(f1, f1, dat, horizon = 3, nsim = 3), "Bayesian")
})

test_that("dinarch_transition with nsim builds one transition per posterior draw", {
  fit_from <- .fake_fit(
    b = 0.35, phi = 27.5, beta = c("(Intercept)" = 0.5), method = "bayes",
    posterior = list(
      b = matrix(c(0.1, 0.6), 2, 1), phi = c(5, 50),
      beta = matrix(c(0, 1), 2, 1, dimnames = list(NULL, "(Intercept)"))
    )
  )
  fit_to <- .fake_fit(b = 0.2, phi = 10, beta = c("(Intercept)" = 1))
  dat <- data.frame(y = rpois(30, 3), time = 1:30)
  np <- dinarch_transition(fit_from, fit_to, dat, horizon = 4, nsim = 20, seed = 1, shape = "linear")
  expect_identical(names(np), c("time", "sim", "b1", "phi", "(Intercept)"))
  expect_equal(sort(unique(np$sim)), 1:20)
  expect_setequal(np$b1[np$time == 31], c(0.1, 0.6))
  expect_equal(np$b1[np$time == 34], rep(0.2, 20))
  expect_identical(
    np,
    dinarch_transition(fit_from, fit_to, dat, horizon = 4, nsim = 20, seed = 1, shape = "linear")
  )
})
