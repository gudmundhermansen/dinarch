# Internal helpers -------------------------------------------------------
# Not exported. Shared by dinarch_fit_ml(), dinarch_fit_bayes(), and
# dinarch_simulate().

#' Check that required columns are present in a data frame / data.table
#' @noRd
.check_columns <- function(data, cols, arg_name = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "`%s` is missing required column(s): %s.", arg_name,
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

#' Add n_lags lag columns of y (by id) to a data.table, in place
#'
#' Adds columns "y1", ..., "y<n_lags>" holding y shifted 1..n_lags steps
#' within each id. `Y` must already have a column literally named "y"
#' (the fit functions select/rename into this canonical layout before
#' calling this helper).
#' @noRd
.add_lag_columns <- function(Y, n_lags, by = "id") {
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

#' Check a vector has no missing (NA) values
#' @noRd
.check_no_na <- function(x, arg_name = deparse(substitute(x))) {
  if (anyNA(x)) {
    stop(sprintf("`%s` contains missing (NA) value(s).", arg_name))
  }
  invisible(x)
}

#' Check values are non-negative integers (a count series)
#' @noRd
.check_nonneg_integer <- function(x, arg_name = deparse(substitute(x)), tol = 1e-8) {
  if (any(x < 0) || any(abs(x - round(x)) > tol)) {
    stop(sprintf("`%s` must be non-negative integers (a count series).", arg_name))
  }
  invisible(x)
}

#' Check a series isn't constant (e.g. all zero) - a degenerate series
#' with no variation for the model to fit
#'
#' Only applies once there's more than one observation - a single-row
#' series is trivially "constant" but that's a distinct, more informative
#' error case (not enough data for any `n_lags`), raised separately.
#' @noRd
.check_has_variation <- function(x, arg_name = deparse(substitute(x))) {
  if (length(x) > 1 && length(unique(x)) <= 1) {
    stop(sprintf(
      "`%s` is constant (every value is %s) - there is nothing to fit.",
      arg_name, format(x[1])
    ))
  }
  invisible(x)
}

#' Check `time` is strictly increasing with a constant step size within
#' each id
#'
#' `Y` must already be sorted by id then time (as
#' .prepare_dinarch_data() does before calling this). Checked per id
#' via plain vector operations (not a data.table grouped `[, ..., by=]`
#' expression) to sidestep any risk of "no visible binding" NOTEs from
#' NSE column references. `has_id = FALSE` (single series without an id
#' column, internally `id = 0L`) drops the "(id ...)" mention from error
#' messages, since the user never supplied an `id` argument.
#' @noRd
.check_time_regularity <- function(Y, has_id = TRUE) {
  for (g in unique(Y$id)) {
    idx <- Y$time[Y$id == g]
    if (length(idx) < 2) next
    d <- diff(idx)
    id_note <- if (has_id) sprintf(" (id %s)", g) else ""
    if (any(d <= 0)) {
      stop(sprintf(
        "`time` must be strictly increasing%s - found a duplicate or decreasing value.",
        id_note
      ))
    }
    if (length(unique(d)) > 1) {
      stop(sprintf(
        "`time` must be evenly spaced with no gaps%s - found step size(s) %s.",
        id_note, paste(sort(unique(d)), collapse = ", ")
      ))
    }
  }
  invisible(Y)
}

#' Central-difference numerical Jacobian of f: R^n -> R^m at x
#'
#' Per-component step size scales with `abs(x)` (relative step, floored at
#' 1e-6) so it works reasonably across parameters of very different
#' magnitude without per-call tuning. Used to delta-method standard errors
#' from the internal unconstrained optimization scale to the natural
#' (b, phi, beta) scale in dinarch_fit_ml().
#' @noRd
.numerical_jacobian <- function(f, x, eps = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (is.null(eps)) eps <- pmax(abs(x) * 1e-4, 1e-6)
  f0 <- f(x)
  J <- matrix(0, nrow = length(f0), ncol = n)
  for (i in seq_len(n)) {
    dx <- numeric(n)
    dx[i] <- eps[i]
    J[, i] <- (f(x + dx) - f(x - dx)) / (2 * eps[i])
  }
  J
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

#' Validate raw `data` and build the canonical id/time/y/covariate panel
#'
#' Shared by .prepare_dinarch_data() (fit time) and dinarch_project()
#' (project time): checks that the `y`/`time`/`id`/`cov_vars` columns exist
#' in `data` (a single combined check, so a missing covariate and a missing
#' `y`/`time`/`id` column produce the same kind of error), contain no NA,
#' that `y` is a non-negative integer count series, and that `time` is
#' strictly increasing with a constant step size (no gaps or duplicates)
#' within each `id`. Returns a data.table with canonical column names
#' (`orig_row`, `id`, `time`, `y`, then `cov_vars`), sorted by id then time.
#' Does *not* check that `y` has variation - only fitting needs that, not
#' seeding/replaying a projection.
#' @noRd
.build_dinarch_panel <- function(data, y, time, id, cov_vars = character(0)) {
  required_cols <- c(y, time, cov_vars)
  if (!is.null(id)) required_cols <- c(required_cols, id)
  .check_columns(data, required_cols)

  .check_no_na(data[[y]], y)
  .check_no_na(data[[time]], time)
  if (!is.null(id)) .check_no_na(data[[id]], id)
  for (v in cov_vars) .check_no_na(data[[v]], v)
  .check_nonneg_integer(data[[y]], y)

  Y <- data.table::data.table(
    orig_row = seq_len(nrow(data)),
    id   = if (!is.null(id)) data[[id]] else 0L,
    time = data[[time]],
    y    = data[[y]]
  )
  for (v in cov_vars) Y[, (v) := data[[v]][orig_row]]

  data.table::setorderv(Y, c("id", "time"))
  .check_time_regularity(Y, has_id = !is.null(id))

  Y
}

#' Add lags, drop rows without full lag history, and build the covariate
#' design matrix from an already-validated panel.
#'
#' Fit-time-only continuation of .build_dinarch_panel(): additionally
#' requires `y` to have some variation (not every value identical, e.g. all
#' zero - there would be nothing to fit), adds n_lags within-id lag
#' columns, drops rows with incomplete lag history, and builds the
#' covariate design matrix `Xcov` from `formula` via stats::model.matrix()
#' (so intercept/interactions/dropping-the-intercept all follow ordinary R
#' formula conventions, e.g. `~ pop_log + gdp`, `~ 0 + pop_log + gdp`,
#' `~ pop_log * gdp`). `formula`'s variables are resolved against the
#' *original* `data` (not the id/time/y-only `Y` table), via the
#' `orig_row` column .build_dinarch_panel() tracked through sorting.
#' `formula`'s environment is reset to `baseenv()` first, so a variable
#' referenced by `formula` but missing from `data` errors clearly instead
#' of `model.matrix()` silently falling back to a same-named variable
#' wherever `formula` happened to be written (base functions like `log()`
#' remain usable).
#' @noRd
.prepare_dinarch_data <- function(data, y, time, id, formula, n_lags) {
  environment(formula) <- baseenv()
  cov_vars <- all.vars(formula)

  Y <- .build_dinarch_panel(data, y, time, id, cov_vars)
  .check_has_variation(Y$y, y)

  .add_lag_columns(Y, n_lags = n_lags, by = "id")

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
    lag_cols = lag_cols,
    n_beta = ncol(Xcov),
    y_vec = Yc$y,
    lag_mat = as.matrix(Yc[, lag_cols, with = FALSE]),
    Xcov = Xcov,
    n = nrow(Yc)
  )
}


#' Optimal alpha estimation
#'
#' Estimate the optimal alpha level per line based on counts of overlapping lines.
#'
#' @param dt Panel data used to plot lines
#' @param group The combination of variables that identify each unique line
#' @param x The x-axis variable
#' @param y The y-axis variable
#' @param c The max alpha target
.alpha_est <- function(dt, group = c("id", "sim"), x = "time", y = "y", c = 0.95) {

  x_bin <- cut(dt[[x]], breaks = 200, labels = FALSE)
  y_bin <- cut(dt[[y]], breaks = 200, labels = FALSE)

  overlap_counts <- unique(
    dt[, group, with = FALSE][, `:=`(x_bin = x_bin, y_bin = y_bin)]
  )[, .N, by = c("x_bin", "y_bin")]

  m <- stats::quantile(overlap_counts$N, c)
  1 - (1 - c)^(1/m)
}
