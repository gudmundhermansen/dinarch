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

#' Solve one transition path from interpolated target quantities
#'
#' Interpolates interpretable target quantities between two parameter sets
#' `pars_from`/`pars_to` (each `list(b, phi, beta)`, the shape returned by
#' .dinarch_pick_params(); `beta` named) with weights `w`, then solves back
#' to parameters: persistence `S = sum(b)` and lag shares `b / S` linear,
#' slopes linear, the log long-run mean `z_ref'beta - log(1 - S)` at
#' reference design row `z_ref` (ordered as `beta_names`) linear, `1/phi`
#' linear, and the `(Intercept)` solved to hold that long-run mean. At
#' `w = 0`/`w = 1` the result reproduces `pars_from`/`pars_to` exactly
#' (up to floating point), whatever `z_ref` is.
#' @return `list(b = <length(w) x k matrix>, phi = <length(w) vector>,
#'   beta = <length(w) x p matrix, colnames beta_names>)`.
#' @noRd
.dinarch_transition_params <- function(pars_from, pars_to, z_ref, w, beta_names) {
  h <- length(w)
  k <- length(pars_from$b)
  p <- length(beta_names)
  ip <- match("(Intercept)", beta_names)

  b_A <- unname(pars_from$b)
  b_B <- unname(pars_to$b)
  beta_A <- unname(pars_from$beta[beta_names])
  beta_B <- unname(pars_to$beta[beta_names])

  # persistence and lag shares
  S_A <- sum(b_A)
  S_B <- sum(b_B)
  S_t <- (1 - w) * S_A + w * S_B
  if (S_A > 0 && S_B > 0) {
    q_A <- b_A / S_A
    q_B <- b_B / S_B
  } else if (S_A > 0) {
    q_A <- q_B <- b_A / S_A
  } else if (S_B > 0) {
    q_A <- q_B <- b_B / S_B
  } else {
    q_A <- q_B <- rep(1 / k, k)
  }
  q_t <- outer(1 - w, q_A) + outer(w, q_B)
  b_t <- S_t * q_t

  # slopes linear; intercept solved to hold the interpolated long-run mean
  beta_t <- matrix(0, h, p, dimnames = list(NULL, beta_names))
  if (p > 1) {
    beta_t[, -ip] <- outer(1 - w, beta_A[-ip]) + outer(w, beta_B[-ip])
  }
  lm_A <- sum(z_ref * beta_A) - log(1 - S_A)
  lm_B <- sum(z_ref * beta_B) - log(1 - S_B)
  lm_t <- (1 - w) * lm_A + w * lm_B
  slope_part <- if (p > 1) drop(beta_t[, -ip, drop = FALSE] %*% z_ref[-ip]) else 0
  beta_t[, ip] <- lm_t + log(1 - S_t) - slope_part

  # dispersion on the 1/phi scale
  phi_t <- 1 / ((1 - w) / pars_from$phi + w / pars_to$phi)

  list(b = b_t, phi = phi_t, beta = beta_t)
}

#' Build a `new_para` table transitioning between two fits' parameters
#'
#' Moves the model from `fit_from`'s fitted parameters to `fit_to`'s over
#' `horizon` steps per id, producing a table in exactly the shape
#' [dinarch_project()]'s `new_para` expects, ready to pass straight
#' through, e.g. `dinarch_project(fit_from, data, horizon = horizon,
#' new_para = dinarch_transition(fit_from, fit_to, data, horizon))`.
#'
#' Rather than interpolating each raw coefficient on its own (which
#' produces artefacts - e.g. the long-run mean `exp(eta) / (1 - sum(b))`
#' dips mid-path between two regimes with the same mean, and `phi` on
#' its natural scale is badly behaved when it's huge), every intermediate
#' step interpolates a complete set of interpretable target quantities
#' with one shared weight path from [dinarch_transition_path()] (`0` to
#' `1`), then solves back to parameters:
#' * persistence `sum(b)`: linear;
#' * lag shares `b_j / sum(b)`: linear (an endpoint with `sum(b) = 0`
#'   borrows the other endpoint's shares);
#' * slopes (every `beta` except `(Intercept)`): linear;
#' * log long-run mean `eta - log(1 - sum(b))`, evaluated at each id's
#'   last observed covariate row in `data`: linear;
#' * `1 / phi` (the NB overdispersion): linear;
#' * `(Intercept)`: solved so the long-run mean at the reference row
#'   equals its interpolated value.
#'
#' Both endpoints reproduce each fit's parameters exactly regardless of
#' the reference row - it only affects the intermediate steps. Since the
#' intercept is what gets solved for, the fits' `formula` must include
#' one. With `shape = "linear"` the last forecast step equals `fit_to`'s
#' parameters exactly; the default sigmoid only gets as far as weight
#' `plogis(range[2])` (about 0.99988 for the default `range`).
#'
#' `fit_from` and `fit_to` must share the same `n_lags` and exactly the
#' same set of `beta` coefficient names (i.e. the same `formula`) - each
#' fit has exactly one parameter set (pooled across ids even for
#' multi-id fits, since `id` only ever affects lag boundaries at fit
#' time, not the parameters themselves), so "transition `b[2]` from A to
#' B" is only well-defined when A and B's parameter spaces line up
#' exactly. `fit_to`'s own coefficients are the only thing used from it -
#' `data` is `fit_from`'s historical data only, used to find each id's
#' last observed `time` to index the forecast window from, and its last
#' observed covariate row (the long-run-mean reference).
#'
#' Each id's forecast window is indexed exactly as
#' [dinarch_project()]'s own `horizon` forecasting does - starting at
#' that id's own `<last observed time> + 1` - so the same shape of
#' transition (steps `1:horizon`) lands at different absolute `time`
#' values across ids with different histories, but the returned table
#' still lines up with what `dinarch_project(fit_from, data, horizon =
#' horizon, ...)` actually simulates (using the same `data`).
#'
#' With `nsim` set, parameter uncertainty is propagated: for each
#' replicate `1:nsim`, an independent posterior draw is taken from each
#' Bayesian fit (an ML fit contributes its point estimate every time) and
#' a full transition is built from that pair, shared across all ids. The
#' output gains a `sim` column; pass the same `nsim` to
#' [dinarch_project()], which applies each replicate's rows to that
#' replicate only.
#'
#' @param fit_from,fit_to `"dinarch_fit"` objects to transition from/to.
#' @param data `fit_from`'s historical data - see [dinarch_project()]'s
#'   `data` argument (same validation, same `y`/`time`/`id` column-name
#'   convention). Must contain every variable the fits' `formula`
#'   references.
#' @param horizon Number of forecast steps - must match the `horizon` you
#'   pass to [dinarch_project()]; a mismatch makes [dinarch_project()]
#'   warn: a shorter transition reverts to `fit_from`'s parameters after
#'   its last step, and a longer one is cut off.
#' @param shape,range Passed to [dinarch_transition_path()] - one shared
#'   interpolation weight path for every target quantity.
#' @param nsim `NULL` (default) builds one transition from both fits'
#'   point estimates. A whole number `>= 1` builds one transition per
#'   replicate from independent posterior draws (requires at least one
#'   Bayesian fit) and adds a `sim` column - see Details.
#' @param seed Optional random seed, set once before any posterior draws.
#' @param y,time,id Column names (strings) in `data` - see
#'   [dinarch_project()].
#'
#' @return A `data.table` with a `time` column (and an `id` column if
#'   `data` has more than one id, and a `sim` column if `nsim` is set),
#'   plus `b1, ..., bk`, `phi`, and one column per `beta` coefficient name.
#'
#' @examples
#' sim_a <- dinarch_simulate(n = 100, b = 0.2, phi = 8, beta = log(3), seed = 1)
#' sim_b <- dinarch_simulate(n = 100, b = 0.6, phi = 8, beta = log(10), seed = 2)
#' dat_a <- data.frame(y = sim_a$y, time = seq_len(100))
#' fit_a <- dinarch_fit_ml(dat_a, y = "y", time = "time", n_lags = 1, vcov = FALSE)
#' fit_b <- dinarch_fit_ml(
#'   data.frame(y = sim_b$y, time = seq_len(100)), y = "y", time = "time",
#'   n_lags = 1, vcov = FALSE
#' )
#' new_para <- dinarch_transition(fit_a, fit_b, dat_a, horizon = 20)
#' proj <- dinarch_project(fit_a, dat_a, horizon = 20, nsim = 100, new_para = new_para)
#'
#' @export
dinarch_transition <- function(fit_from, fit_to, data, horizon,
                               shape = c("sigmoid", "linear"),
                               range = c(-3, 9),
                               nsim = NULL, seed = NULL,
                               y = "y", time = "time", id = NULL) {
  shape <- match.arg(shape)
  if (!inherits(fit_from, "dinarch_fit") || !inherits(fit_to, "dinarch_fit")) {
    stop("`fit_from` and `fit_to` must both be \"dinarch_fit\" objects.")
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
  if (!"(Intercept)" %in% beta_from) {
    stop("dinarch_transition() needs an `(Intercept)` term in the fits' `formula` - the intercept is solved for to hold the interpolated long-run mean.")
  }
  if (!is.null(nsim)) {
    stopifnot(is.numeric(nsim), length(nsim) == 1, nsim >= 1, nsim == round(nsim))
  }
  is_bayes <- function(f) identical(f$method, "bayes") && !is.null(f$posterior)
  if (!is.null(nsim) && !is_bayes(fit_from) && !is_bayes(fit_to)) {
    stop("`nsim` requires at least one Bayesian fit (ML fits have no posterior draws); leave `nsim = NULL` for ML fits.")
  }

  if (!is.null(seed)) set.seed(seed)

  n_lags <- fit_from$n_lags
  b_cols <- paste0("b", seq_len(n_lags))
  w <- dinarch_transition_path(0, 1, n = horizon, shape = shape, range = range)

  # reference design rows: each id's last observed row of `data`
  formula <- fit_from$formula
  environment(formula) <- baseenv()
  panel <- .build_dinarch_panel(data, y, time, id, all.vars(formula))
  X <- stats::model.matrix(formula, data = as.data.frame(data)[panel$orig_row, , drop = FALSE])
  if (!identical(colnames(X), beta_from)) {
    stop(sprintf(
      "`data` produces design columns (%s) that don't match the fits' `beta` names (%s).",
      paste(colnames(X), collapse = ", "), paste(beta_from, collapse = ", ")
    ))
  }
  ids <- unique(panel$id)
  multi_id <- length(ids) > 1
  id_ref <- lapply(ids, function(g) {
    last_row <- max(which(panel$id == g))
    list(id = g, last_time = panel$time[last_row], z_ref = X[last_row, ])
  })

  # posterior draws may come back unnamed (e.g. rstan::extract()), so name
  # `beta` positionally by the fit's own coefficient order
  pick <- function(f, use_draw) {
    pars <- .dinarch_pick_params(f, use_draw)
    names(pars$beta) <- names(f$coefficients$beta)
    pars
  }

  build <- function(pars_from, pars_to) {
    data.table::rbindlist(lapply(id_ref, function(r) {
      tp <- .dinarch_transition_params(pars_from, pars_to, r$z_ref, w, beta_from)
      tbl <- data.table::data.table(time = r$last_time + seq_len(horizon))
      if (multi_id) tbl[, "id" := r$id]
      for (j in seq_len(n_lags)) tbl[, (b_cols[j]) := tp$b[, j]]
      tbl[, "phi" := tp$phi]
      for (nm in beta_from) tbl[, (nm) := tp$beta[, nm]]
      tbl
    }), use.names = TRUE)
  }

  if (is.null(nsim)) {
    out <- build(pick(fit_from, NULL), pick(fit_to, NULL))
    data.table::setcolorder(out, if (multi_id) c("id", "time") else "time")
  } else {
    out <- data.table::rbindlist(lapply(seq_len(nsim), function(s) {
      pars_from <- pick(fit_from, if (is_bayes(fit_from)) "random" else NULL)
      pars_to   <- pick(fit_to, if (is_bayes(fit_to)) "random" else NULL)
      tbl <- build(pars_from, pars_to)
      tbl[, "sim" := s]
      tbl
    }), use.names = TRUE)
    data.table::setcolorder(out, c(if (multi_id) "id", "time", "sim"))
  }
  out[]
}
