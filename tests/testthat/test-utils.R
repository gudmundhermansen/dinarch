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

test_that(".as_coef_matrix builds an n x k matrix from scalar, vector, or matrix input", {
  m1 <- .as_coef_matrix(0.3, n = 4)
  expect_equal(dim(m1), c(4, 1))
  expect_true(all(m1 == 0.3))

  m2 <- .as_coef_matrix(c(0.3, 0.1), n = 4, k = 2)
  expect_equal(dim(m2), c(4, 2))
  expect_true(all(m2[, 1] == 0.3) && all(m2[, 2] == 0.1))

  b_mat <- matrix(runif(8), nrow = 4, ncol = 2)
  m3 <- .as_coef_matrix(b_mat, n = 4, k = 2)
  expect_identical(m3, b_mat)

  expect_error(.as_coef_matrix(c(0.3, 0.1, 0.2), n = 4, k = 2), "length")
})

test_that(".add_lag_columns adds correct within-group lags", {
  Y <- data.table::data.table(
    group = c(1, 1, 1, 2, 2),
    y     = c(10, 20, 30, 5, 6)
  )
  .add_lag_columns(Y, n_lags = 2, by = "group")

  expect_equal(Y$y1, c(NA, 10, 20, NA, 5))
  expect_equal(Y$y2, c(NA, NA, 10, NA, NA))
})

test_that(".check_unit_interval and .check_positive validate ranges", {
  expect_silent(.check_unit_interval(0.5))
  expect_error(.check_unit_interval(1.2), "between 0 and 1")
  expect_silent(.check_positive(c(1, 2, 3)))
  expect_error(.check_positive(c(1, -1)), "positive")
})

test_that(".check_no_na and .check_nonneg_integer validate y-like inputs", {
  expect_silent(.check_no_na(c(1, 2, 3)))
  expect_error(.check_no_na(c(1, NA, 3)), "missing")

  expect_silent(.check_nonneg_integer(c(0, 1, 5)))
  expect_error(.check_nonneg_integer(c(1, -2, 3)), "non-negative integers")
  expect_error(.check_nonneg_integer(c(1, 2.5, 3)), "non-negative integers")
})

test_that(".check_has_variation rejects a constant series", {
  expect_silent(.check_has_variation(c(0, 1, 0, 2)))
  expect_error(.check_has_variation(c(0, 0, 0)), "constant")
  expect_error(.check_has_variation(c(3, 3, 3)), "constant")
})

test_that(".prepare_dinarch_data errors when y is constant (e.g. all zero)", {
  dat <- data.frame(y = rep(0, 30), index = seq_len(30))
  expect_error(
    dinarch_fit_ml(dat, y = "y", index = "index", n_lags = 1),
    "constant"
  )
})

test_that(".check_index_regularity accepts evenly-spaced index per group and rejects gaps/duplicates/decreases", {
  ok <- data.table::data.table(group = c(1, 1, 1, 2, 2), index = c(1, 2, 3, 10, 12))
  expect_silent(.check_index_regularity(ok))  # group 2's own step (2) just needs to be constant

  gap <- data.table::data.table(group = c(1, 1, 1), index = c(1, 2, 4))
  expect_error(.check_index_regularity(gap), "evenly spaced")

  dup <- data.table::data.table(group = c(1, 1, 1), index = c(1, 2, 2))
  expect_error(.check_index_regularity(dup), "strictly increasing")

  decreasing <- data.table::data.table(group = c(1, 1, 1), index = c(1, 3, 2))
  expect_error(.check_index_regularity(decreasing), "strictly increasing")
})

test_that(".prepare_dinarch_data builds correct defaults with no group and formula = ~1", {
  dat <- data.frame(y = c(5, 6, 7, 8), t = 1:4)
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                                 formula = ~1, n_lags = 1)

  expect_equal(prep$n, 3)  # one row lost to the lag
  expect_equal(prep$y_vec, c(6, 7, 8))
  expect_equal(dim(prep$Xcov), c(3, 1))
  expect_identical(colnames(prep$Xcov), "(Intercept)")
  expect_equal(prep$n_beta, 1)
  expect_equal(as.numeric(prep$lag_mat), c(5, 6, 7))
})

test_that(".prepare_dinarch_data pools multiple groups without crossing boundaries", {
  dat <- data.frame(
    y   = c(10, 20, 30, 5, 6),
    t   = c(1, 2, 3, 1, 2),
    grp = c("A", "A", "A", "B", "B")
  )
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = "grp",
                                 formula = ~1, n_lags = 1)

  expect_equal(prep$n, 3)  # one row lost per group (2 groups, 5 rows -> 3)
  expect_equal(prep$y_vec, c(20, 30, 6))
  expect_equal(as.numeric(prep$lag_mat), c(10, 20, 5))
})

test_that(".prepare_dinarch_data builds Xcov from formula, correctly aligned to surviving rows", {
  dat <- data.frame(y = c(1, 2, 3), t = 1:3, gdp = c(10, 2, 3))
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                                 formula = ~gdp, n_lags = 1)

  expect_equal(prep$n_beta, 2)  # intercept + gdp
  expect_identical(colnames(prep$Xcov), c("(Intercept)", "gdp"))
  # first row (t = 1) is dropped by the n_lags = 1 lag, so gdp column
  # should align to the surviving rows (t = 2, 3), not the first two raw rows
  expect_equal(unname(prep$Xcov[, "gdp"]), c(2, 3))
})

test_that(".prepare_dinarch_data supports dropping the intercept via formula", {
  dat <- data.frame(y = c(1, 2, 3), t = 1:3, gdp = c(10, 2, 3))
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                                 formula = ~ 0 + gdp, n_lags = 1)

  expect_identical(colnames(prep$Xcov), "gdp")
})

test_that(".prepare_dinarch_data errors on missing columns, empty results, and bad formulas", {
  dat <- data.frame(y = 1:3, t = 1:3)
  expect_error(
    .prepare_dinarch_data(dat, y = "yy", index = "t", group = NULL,
                           formula = ~1, n_lags = 1),
    "missing required column"
  )

  dat2 <- data.frame(y = 1, t = 1)  # single row: no row survives a n_lags = 1 lag
  expect_error(
    .prepare_dinarch_data(dat2, y = "y", index = "t", group = NULL,
                           formula = ~1, n_lags = 1),
    "No usable observations"
  )

  expect_error(
    .prepare_dinarch_data(dat, y = "y", index = "t", group = NULL,
                           formula = ~nonexistent, n_lags = 1)
  )
})

test_that(".prepare_dinarch_data errors on NA in y/index/group/covariates", {
  base <- data.frame(y = 1:5, t = 1:5, grp = "A", gdp = c(1, 2, 3, 4, 5))

  dat_na_y <- base; dat_na_y$y[3] <- NA
  expect_error(
    .prepare_dinarch_data(dat_na_y, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "missing"
  )

  dat_na_t <- base; dat_na_t$t[3] <- NA
  expect_error(
    .prepare_dinarch_data(dat_na_t, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "missing"
  )

  dat_na_grp <- base; dat_na_grp$grp[3] <- NA
  expect_error(
    .prepare_dinarch_data(dat_na_grp, y = "y", index = "t", group = "grp", formula = ~1, n_lags = 1),
    "missing"
  )

  dat_na_gdp <- base; dat_na_gdp$gdp[3] <- NA
  expect_error(
    .prepare_dinarch_data(dat_na_gdp, y = "y", index = "t", group = NULL, formula = ~gdp, n_lags = 1),
    "missing"
  )
})

test_that(".prepare_dinarch_data errors on non-count y and on irregular index", {
  dat_negative <- data.frame(y = c(1, -2, 3), t = 1:3)
  expect_error(
    .prepare_dinarch_data(dat_negative, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "non-negative integers"
  )

  dat_fractional <- data.frame(y = c(1, 2.5, 3), t = 1:3)
  expect_error(
    .prepare_dinarch_data(dat_fractional, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "non-negative integers"
  )

  dat_gap <- data.frame(y = 1:4, t = c(1, 2, 4, 5))  # skips 3
  expect_error(
    .prepare_dinarch_data(dat_gap, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "evenly spaced"
  )

  dat_dup <- data.frame(y = 1:4, t = c(1, 2, 2, 3))
  expect_error(
    .prepare_dinarch_data(dat_dup, y = "y", index = "t", group = NULL, formula = ~1, n_lags = 1),
    "strictly increasing"
  )
})

test_that("groups may each have their own (internally constant) index step size", {
  dat <- data.frame(
    y   = c(1, 2, 3, 4, 5, 6),
    t   = c(2000, 2001, 2002, 5, 10, 15),  # group A: annual; group B: step of 5
    grp = c("A", "A", "A", "B", "B", "B")
  )
  prep <- .prepare_dinarch_data(dat, y = "y", index = "t", group = "grp", formula = ~1, n_lags = 1)
  expect_equal(prep$n, 4)  # one row lost per group
})
