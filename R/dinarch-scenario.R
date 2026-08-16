# Scenario construction -----------------------------------------------------

#' Build a smooth parameter transition path
#'
#' Interpolates from `from` to `to` over `n` steps, for constructing
#' scenario parameter paths (e.g. `b` or `phi` trajectories) to feed into
#' [dinarch_simulate()] or [simulate.dinarch_fit()].
#'
#' @param from,to Starting and ending values (scalars). Works the same way
#'   whether `to > from` (increasing path) or `to < from` (decreasing
#'   path).
#' @param n Number of steps in the path.
#' @param shape `"sigmoid"` (default) for a smooth S-curve transition, or
#'   `"linear"` for straight-line interpolation.
#' @param range For `shape = "sigmoid"`, the domain over which the sigmoid
#'   is evaluated before rescaling to `[from, to]`. The default
#'   `c(-3, 9)` is asymmetric - it front-loads the transition, reaching
#'   close to `to` well before the final step. Use a symmetric range like
#'   `c(-6, 6)` for a transition that eases in and out evenly at both
#'   ends.
#'
#' @return A numeric vector of length `n`.
#'
#' @examples
#' # front-loaded transition (the default range)
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20)
#'
#' # symmetric sigmoid ramp
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20, range = c(-6, 6))
#'
#' # straight-line interpolation
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20, shape = "linear")
#'
#' @export
dinarch_transition_path <- function(from, to, n,
                                     shape = c("sigmoid", "linear"),
                                     range = c(-3, 9)) {
  shape <- match.arg(shape)
  stopifnot(is.numeric(n), length(n) == 1, n >= 1, n == round(n))
  stopifnot(length(from) == 1, length(to) == 1)

  if (shape == "linear") {
    return(seq(from, to, length.out = n))
  }

  x <- seq(range[1], range[2], length.out = n)
  w <- stats::plogis(x)
  from + (to - from) * w
}

#' Build a `new_para` table transitioning between two fits' parameters
#'
#' Interpolates every parameter (`b_1, ..., b_k`, `phi`, and each `beta`
#' coefficient) from `fit_from`'s fitted values to `fit_to`'s, over
#' `horizon` steps per group, via [dinarch_transition_path()] - producing
#' a table in exactly the shape [dinarch_project()]'s `new_para` expects,
#' ready to pass straight through, e.g.
#' `dinarch_project(fit_from, horizon = horizon, new_para =
#' dinarch_transition_new_para(fit_from, fit_to, horizon))`. At the last
#' forecast step, the parameters equal `fit_to`'s exactly.
#'
#' `fit_from` and `fit_to` must share the same `n_lags` and exactly the
#' same set of `beta` coefficient names (i.e. the same `formula`) - each
#' fit has exactly one parameter set (pooled across groups even for
#' multi-group fits, since `group` only ever affects lag boundaries at fit
#' time, not the parameters themselves), so "transition `b[2]` from A to
#' B" is only well-defined when A and B's parameter spaces line up
#' exactly. `fit_to`'s own `data`/`group` structure is never used - only
#' its `$coefficients`.
#'
#' Each group's forecast window is indexed exactly as
#' [dinarch_project()]'s own `horizon` forecasting does - starting at
#' that group's own `<last observed index> + 1` - so the same shape of
#' transition (steps `1:horizon`) lands at different absolute `index`
#' values across groups with different histories, but the returned table
#' still lines up with what `dinarch_project(fit_from, horizon = horizon,
#' ...)` actually simulates.
#'
#' @param fit_from,fit_to `"dinarch_fit"` objects to transition from/to.
#' @param horizon Number of forecast steps - must match the `horizon` you
#'   pass to [dinarch_project()].
#' @param shape,range Passed to [dinarch_transition_path()] - one shared
#'   interpolation shape for every parameter.
#'
#' @return A `data.table` with an `index` column (and a `group` column if
#'   `fit_from` has more than one group), plus `b1, ..., bk`, `phi`, and
#'   one column per `beta` coefficient name.
#'
#' @examples
#' sim_a <- dinarch_simulate(n = 100, b = 0.2, phi = 8, beta = log(3), seed = 1)
#' sim_b <- dinarch_simulate(n = 100, b = 0.6, phi = 8, beta = log(10), seed = 2)
#' fit_a <- dinarch_fit_ml(
#'   data.frame(y = sim_a$y, index = seq_len(100)), y = "y", index = "index",
#'   n_lags = 1, vcov = FALSE
#' )
#' fit_b <- dinarch_fit_ml(
#'   data.frame(y = sim_b$y, index = seq_len(100)), y = "y", index = "index",
#'   n_lags = 1, vcov = FALSE
#' )
#' new_para <- dinarch_transition_new_para(fit_a, fit_b, horizon = 20)
#' proj <- dinarch_project(fit_a, horizon = 20, nsim = 100, new_para = new_para)
#'
#' @export
dinarch_transition_new_para <- function(fit_from, fit_to, horizon,
                                         shape = c("sigmoid", "linear"),
                                         range = c(-3, 9)) {
  shape <- match.arg(shape)
  if (!inherits(fit_from, "dinarch_fit") || !inherits(fit_to, "dinarch_fit")) {
    stop("`fit_from` and `fit_to` must both be \"dinarch_fit\" objects.")
  }
  if (is.null(fit_from$data)) {
    stop(
      "`fit_from` has no stored `data` (fit with an older version of the ",
      "package?) - refit with the current dinarch_fit_ml()/",
      "dinarch_fit_bayes() to use dinarch_transition_new_para()."
    )
  }
  if (!identical(fit_from$n_lags, fit_to$n_lags)) {
    stop(sprintf(
      "`fit_from` and `fit_to` must have the same `n_lags` (%d vs %d).",
      fit_from$n_lags, fit_to$n_lags
    ))
  }
  beta_from <- names(fit_from$coefficients$beta)
  beta_to   <- names(fit_to$coefficients$beta)
  if (!setequal(beta_from, beta_to)) {
    stop(
      "`fit_from` and `fit_to` must have exactly the same `beta` coefficients ",
      "(i.e. the same `formula`) - found only in `fit_from`: ",
      paste(setdiff(beta_from, beta_to), collapse = ", "), "; only in `fit_to`: ",
      paste(setdiff(beta_to, beta_from), collapse = ", "), "."
    )
  }
  stopifnot(is.numeric(horizon), length(horizon) == 1, horizon >= 1, horizon == round(horizon))

  n_lags <- fit_from$n_lags
  b_cols <- paste0("b", seq_len(n_lags))
  path <- function(from, to) {
    dinarch_transition_path(from, to, n = horizon, shape = shape, range = range)
  }

  b_paths <- stats::setNames(
    lapply(seq_len(n_lags), function(j) path(fit_from$coefficients$b[j], fit_to$coefficients$b[j])),
    b_cols
  )
  phi_path <- path(fit_from$coefficients$phi, fit_to$coefficients$phi)
  beta_paths <- stats::setNames(
    lapply(beta_from, function(nm) path(fit_from$coefficients$beta[[nm]], fit_to$coefficients$beta[[nm]])),
    beta_from
  )

  groups <- unique(fit_from$data$group)
  multi_group <- length(groups) > 1

  out <- data.table::rbindlist(lapply(groups, function(g) {
    last_idx <- max(fit_from$data$index[fit_from$data$group == g])
    tbl <- data.table::data.table(index = last_idx + seq_len(horizon))
    if (multi_group) tbl[, "group" := g]
    for (nm in names(b_paths)) tbl[, (nm) := b_paths[[nm]]]
    tbl[, "phi" := phi_path]
    for (nm in names(beta_paths)) tbl[, (nm) := beta_paths[[nm]]]
    tbl
  }))

  data.table::setcolorder(out, if (multi_group) c("group", "index") else "index")
  out[]
}
