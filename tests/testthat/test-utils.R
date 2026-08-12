test_that(".check_columns detects missing columns", {
  d <- data.frame(a = 1, b = 2)
  expect_silent(.check_columns(d, c("a", "b")))
  expect_error(.check_columns(d, c("a", "c")), "missing required column")
})

test_that(".recycle repeats scalars and validates vector length", {
  expect_equal(.recycle(2, 5), rep(2, 5))
  expect_equal(.recycle(c(1, 2, 3), 3), c(1, 2, 3))
  expect_error(.recycle(c(1, 2), 5), "length")
})

test_that(".as_lag_matrix builds an n x k matrix from scalar, vector, or matrix input", {
  m1 <- .as_lag_matrix(0.3, n = 4)
  expect_equal(dim(m1), c(4, 1))
  expect_true(all(m1 == 0.3))

  m2 <- .as_lag_matrix(c(0.3, 0.1), n = 4, k = 2)
  expect_equal(dim(m2), c(4, 2))
  expect_true(all(m2[, 1] == 0.3) && all(m2[, 2] == 0.1))

  b_mat <- matrix(runif(8), nrow = 4, ncol = 2)
  m3 <- .as_lag_matrix(b_mat, n = 4, k = 2)
  expect_identical(m3, b_mat)

  expect_error(.as_lag_matrix(c(0.3, 0.1, 0.2), n = 4, k = 2), "length")
})

test_that(".add_lag_columns adds correct within-group lags", {
  Y <- data.table::data.table(
    group = c(1, 1, 1, 2, 2),
    y     = c(10, 20, 30, 5, 6)
  )
  .add_lag_columns(Y, k = 2, by = "group")

  expect_equal(Y$y1, c(NA, 10, 20, NA, 5))
  expect_equal(Y$y2, c(NA, NA, 10, NA, NA))
})

test_that(".check_unit_interval and .check_positive validate ranges", {
  expect_silent(.check_unit_interval(0.5))
  expect_error(.check_unit_interval(1.2), "between 0 and 1")
  expect_silent(.check_positive(c(1, 2, 3)))
  expect_error(.check_positive(c(1, -1)), "positive")
})

test_that(".prepare_dinarch_data builds correct defaults with no group/population/covariates", {
  dat <- data.frame(y = c(5, 6, 7, 8), t = 1:4)
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                                 population = NULL, covariates = NULL, k = 1)

  expect_equal(prep$n, 3)  # one row lost to the lag
  expect_equal(prep$y_vec, c(6, 7, 8))
  expect_equal(prep$pop_vec, rep(1, 3))
  expect_equal(dim(prep$Xcov), c(3, 0))
  expect_false(prep$has_covariates)
  expect_equal(prep$n_beta, 0)
  expect_equal(as.numeric(prep$lag_mat), c(5, 6, 7))
})

test_that(".prepare_dinarch_data pools multiple groups without crossing boundaries", {
  dat <- data.frame(
    y   = c(10, 20, 30, 5, 6),
    t   = c(1, 2, 3, 1, 2),
    grp = c("A", "A", "A", "B", "B")
  )
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = "grp",
                                 population = NULL, covariates = NULL, k = 1)

  expect_equal(prep$n, 3)  # one row lost per group (2 groups, 5 rows -> 3)
  expect_equal(prep$y_vec, c(20, 30, 6))
  expect_equal(as.numeric(prep$lag_mat), c(10, 20, 5))
})

test_that(".prepare_dinarch_data picks up population and covariates by name", {
  dat <- data.frame(y = c(1, 2, 3), t = 1:3, pop = c(100, 100, 100), gdp = c(1, 2, 3))
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                                 population = "pop", covariates = "gdp", k = 1)

  expect_equal(prep$pop_vec, c(100, 100))
  expect_true(prep$has_covariates)
  expect_equal(prep$n_beta, 1)
  expect_equal(as.numeric(prep$Xcov), c(2, 3))
})

test_that(".prepare_dinarch_data errors on missing columns and on empty results", {
  dat <- data.frame(y = 1:3, t = 1:3)
  expect_error(
    .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                           population = "pop", covariates = NULL, k = 1),
    "missing required column"
  )

  dat2 <- data.frame(y = 1, t = 1)  # single row: no row survives a k = 1 lag
  expect_error(
    .prepare_dinarch_data(dat2, y = "y", index = "t", group = NULL,
                           population = NULL, covariates = NULL, k = 1),
    "No usable observations"
  )
})
