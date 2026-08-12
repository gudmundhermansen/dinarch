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
#' Used throughout to let `a`, `phi`, `population`, etc. be supplied either
#' as a single constant or as a full length-n path.
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

#' Coerce a lag-coefficient argument to a canonical n x k matrix
#'
#' `b` may be supplied as: a scalar (k = 1, constant over time), a length-k
#' vector (constant lag coefficients over time), or an n x k matrix (fully
#' time-varying lag coefficients). Always returns an n x k matrix.
#' @noRd
.as_lag_matrix <- function(b, n, k = NULL) {
  if (is.matrix(b)) {
    if (nrow(b) != n) {
      stop(sprintf("`b` has %d rows but `n = %d`.", nrow(b), n))
    }
    if (!is.null(k) && ncol(b) != k) {
      stop(sprintf("`b` has %d columns but `k = %d`.", ncol(b), k))
    }
    return(b)
  }

  k <- if (is.null(k)) length(b) else k
  if (length(b) != k) {
    stop(sprintf(
      "`b` must have length %d (= k), or be an n x k matrix.", k
    ))
  }
  matrix(b, nrow = n, ncol = k, byrow = TRUE)
}

#' Add k lag columns of y (by group) to a data.table, in place
#'
#' Adds columns "y1", ..., "yk" holding y shifted 1..k steps within each
#' group. `Y` must already have a column literally named "y" (the fit
#' functions select/rename into this canonical layout before calling this
#' helper). Uses data.table::shift(), a cleaner and faster replacement for
#' the manual shift-loop used in earlier drafts of this code.
#' @noRd
.add_lag_columns <- function(Y, k, by = "group") {
  lag_cols <- paste0("y", seq_len(k))
  Y[, (lag_cols) := data.table::shift(y, n = seq_len(k)), by = by]
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

#' Build the canonical panel data.table, add lags, and drop rows without
#' full lag history.
#'
#' Shared data-prep pipeline for dinarch_fit_ml() and dinarch_fit_bayes():
#' validates required columns, builds a data.table with canonical column
#' names (group/index/y/population/covariates), sorts by group then index,
#' adds k within-group lag columns, and drops rows with incomplete lag
#' history. Returns the pieces both fit functions need in matrix/vector
#' form, so neither has to repeat this logic.
#' @noRd
.prepare_dinarch_data <- function(data, y, index, group, population,
                                   covariates, k) {
  required_cols <- c(y, index)
  if (!is.null(group)) required_cols <- c(required_cols, group)
  if (!is.null(population)) required_cols <- c(required_cols, population)
  if (!is.null(covariates)) required_cols <- c(required_cols, covariates)
  .check_columns(data, required_cols)

  Y <- data.table::data.table(
    group = if (!is.null(group)) data[[group]] else 0L,
    index = data[[index]],
    y     = data[[y]]
  )
  Y[, population := if (!is.null(population)) data[[population]] else 1]

  has_covariates <- !is.null(covariates)
  if (has_covariates) {
    for (cov_name in covariates) {
      Y[, (cov_name) := data[[cov_name]]]
    }
  }

  data.table::setorderv(Y, c("group", "index"))
  .add_lag_columns(Y, k = k, by = "group")

  lag_cols <- paste0("y", seq_len(k))
  Yc <- stats::na.omit(Y, cols = lag_cols)

  if (nrow(Yc) == 0) {
    stop("No usable observations after constructing lags (check `k` against series lengths).")
  }

  n_beta <- if (has_covariates) length(covariates) else 0L

  list(
    Yc = Yc,
    lag_cols = lag_cols,
    has_covariates = has_covariates,
    n_beta = n_beta,
    y_vec = Yc$y,
    pop_vec = Yc$population,
    lag_mat = as.matrix(Yc[, lag_cols, with = FALSE]),
    Xcov = if (has_covariates) as.matrix(Yc[, covariates, with = FALSE]) else matrix(numeric(0), nrow = nrow(Yc), ncol = 0),
    n = nrow(Yc)
  )
}
