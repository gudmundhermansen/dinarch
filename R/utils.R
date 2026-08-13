# Internal helpers -------------------------------------------------------
# Not exported. Shared by dinarch_fit_ml(), dinarch_fit_bayes(), and
# dinarch_simulate().

#' Check that required columns are present in a data frame / data.table
#' @noRd
.check_columns <- function(data, cols) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "`data` is missing required column(s): %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }
  invisible(data)
}

#' Recycle a scalar to length n, or check that a vector already has length n
#'
#' Used throughout to let `phi`, `beta`, etc. be supplied either as a
#' single constant or as a full length-n path.
#' @noRd
.recycle <- function(x, n, arg_name = deparse(substitute(x))) {
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(sprintf(
      "`%s` must have length 1 or %d, not %d.", arg_name, n, length(x)
    ))
  }
  x
}

#' Coerce a coefficient-path argument (`b` or `beta`) to a canonical
#' n x k matrix
#'
#' Used for both `b` (k = n_lags) and `beta` (k = ncol(Xcov)) in
#' dinarch_simulate() - either may be supplied as: a scalar (k = 1,
#' constant over time), a length-k vector (constant coefficients over
#' time), or an n x k matrix (fully time-varying coefficients). Always
#' returns an n x k matrix.
#' @noRd
.as_coef_matrix <- function(x, n, k = NULL, arg_name = deparse(substitute(x))) {
  if (is.matrix(x)) {
    if (nrow(x) != n) {
      stop(sprintf("`%s` has %d rows but `n = %d`.", arg_name, nrow(x), n))
    }
    if (!is.null(k) && ncol(x) != k) {
      stop(sprintf("`%s` has %d columns but expected %d.", arg_name, ncol(x), k))
    }
    return(x)
  }

  k <- if (is.null(k)) length(x) else k
  if (length(x) != k) {
    stop(sprintf(
      "`%s` must have length %d, or be an n x %d matrix.", arg_name, k, k
    ))
  }
  matrix(x, nrow = n, ncol = k, byrow = TRUE)
}

#' Add n_lags lag columns of y (by group) to a data.table, in place
#'
#' Adds columns "y1", ..., "y<n_lags>" holding y shifted 1..n_lags steps
#' within each group. `Y` must already have a column literally named "y"
#' (the fit functions select/rename into this canonical layout before
#' calling this helper).
#' @noRd
.add_lag_columns <- function(Y, n_lags, by = "group") {
  lag_cols <- paste0("y", seq_len(n_lags))
  Y[, (lag_cols) := data.table::shift(y, n = seq_len(n_lags)), by = by]
  Y
}

#' Check values lie strictly in (0, 1)
#' @noRd
.check_unit_interval <- function(p, arg_name = deparse(substitute(p))) {
  if (any(p <= 0 | p >= 1, na.rm = TRUE)) {
    stop(sprintf("`%s` must be strictly between 0 and 1.", arg_name))
  }
  invisible(p)
}

#' Check values are strictly positive
#' @noRd
.check_positive <- function(x, arg_name = deparse(substitute(x))) {
  if (any(x <= 0, na.rm = TRUE)) {
    stop(sprintf("`%s` must be strictly positive.", arg_name))
  }
  invisible(x)
}

#' Pick (b, phi, beta) parameters from a dinarch_fit for one simulation
#' replicate
#'
#' `use_draw = NULL` (default) uses the point estimate (posterior mean for
#' Bayesian fits, MLE for ML fits). An integer pins one specific posterior
#' draw (Bayesian fits only). `"random"` samples a fresh posterior draw
#' (Bayesian fits only) - used by dinarch_project() to propagate full
#' posterior uncertainty across `nsim` replicates. Shared by
#' simulate.dinarch_fit() and dinarch_project() so the use_draw semantics
#' can't drift apart between the two.
#' @noRd
.dinarch_pick_params <- function(object, use_draw = NULL) {
  is_bayes <- identical(object$method, "bayes") && !is.null(object$posterior)

  if (is.null(use_draw)) {
    return(list(b = object$coefficients$b, phi = object$coefficients$phi, beta = object$coefficients$beta))
  }

  if (!is_bayes) {
    stop("`use_draw` is only available for Bayesian fits (this fit has no stored posterior draws).")
  }

  n_draws <- length(object$posterior$phi)
  draw <- if (identical(use_draw, "random")) sample.int(n_draws, 1) else use_draw
  if (!is.numeric(draw) || length(draw) != 1 || draw < 1 || draw > n_draws || draw != round(draw)) {
    stop(sprintf("`use_draw` must be \"random\" or an integer between 1 and %d.", n_draws))
  }

  list(
    b = object$posterior$b[draw, ],
    phi = object$posterior$phi[draw],
    beta = if (length(object$coefficients$beta) > 0) object$posterior$beta[draw, ] else numeric(0)
  )
}

#' Build the canonical panel data.table, add lags, drop rows without full
#' lag history, and build the covariate design matrix.
#'
#' Shared data-prep pipeline for dinarch_fit_ml() and dinarch_fit_bayes():
#' validates required columns, builds a data.table with canonical column
#' names (group/index/y), sorts by group then index, adds n_lags
#' within-group lag columns, drops rows with incomplete lag history, and
#' builds the covariate design matrix `Xcov` from `formula` via
#' stats::model.matrix() (so intercept/interactions/dropping-the-intercept
#' all follow ordinary R formula conventions, e.g. `~ pop_log + gdp`,
#' `~ 0 + pop_log + gdp`, `~ pop_log * gdp`). `formula`'s variables are
#' resolved against the *original* `data` (not the group/index/y-only `Y`
#' table), so an explicit `orig_row` column tracks row identity through
#' sorting and lag-based row-dropping. `formula`'s environment is reset to
#' `baseenv()` first, so a variable referenced by `formula` but missing
#' from `data` errors clearly instead of `model.matrix()` silently
#' falling back to a same-named variable wherever `formula` happened to
#' be written (base functions like `log()` remain usable). Also returns
#' `full_data` - group/index/y plus the raw (untransformed) covariate
#' columns `formula` references, for *every* row including the first
#' `n_lags` per group that get dropped from `Yc` - this is what fit
#' functions store as `fit$data` so `dinarch_project()` can seed forward
#' simulation or replay in-sample without the caller re-supplying
#' history. Returns the pieces both fit functions need in matrix/vector
#' form, so neither has to repeat this logic.
#' @noRd
.prepare_dinarch_data <- function(data, y, index, group, formula, n_lags) {
  environment(formula) <- baseenv()
  required_cols <- c(y, index)
  if (!is.null(group)) required_cols <- c(required_cols, group)
  .check_columns(data, required_cols)

  Y <- data.table::data.table(
    orig_row = seq_len(nrow(data)),
    group = if (!is.null(group)) data[[group]] else 0L,
    index = data[[index]],
    y     = data[[y]]
  )

  data.table::setorderv(Y, c("group", "index"))

  cov_vars <- all.vars(formula)
  missing_cov_vars <- setdiff(cov_vars, names(data))
  if (length(missing_cov_vars) > 0) {
    stop(sprintf(
      "`formula` references variable(s) not present in `data`: %s.",
      paste(missing_cov_vars, collapse = ", ")
    ))
  }
  for (v in cov_vars) {
    Y[, (v) := data[[v]][orig_row]]
  }
  full_data <- Y[, !"orig_row", with = FALSE]

  .add_lag_columns(Y, n_lags = n_lags, by = "group")

  lag_cols <- paste0("y", seq_len(n_lags))
  Yc <- stats::na.omit(Y, cols = lag_cols)

  if (nrow(Yc) == 0) {
    stop("No usable observations after constructing lags (check `n_lags` against series lengths).")
  }

  Xcov <- stats::model.matrix(formula, data = as.data.frame(data)[Yc$orig_row, , drop = FALSE])
  if (nrow(Xcov) != nrow(Yc)) {
    stop(sprintf(
      "`formula` produced %d rows but %d observations survive lagging - ",
      nrow(Xcov), nrow(Yc)
    ), "check `data` for missing values in the variables referenced by `formula`.")
  }
  if (!all(is.finite(Xcov))) {
    stop("`formula` produced non-finite values - check `data` for missing or infinite covariate values.")
  }

  list(
    Yc = Yc,
    full_data = full_data,
    lag_cols = lag_cols,
    n_beta = ncol(Xcov),
    y_vec = Yc$y,
    lag_mat = as.matrix(Yc[, lag_cols, with = FALSE]),
    Xcov = Xcov,
    n = nrow(Yc)
  )
}
