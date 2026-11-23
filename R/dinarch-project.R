# Projection (forecast / in-sample replay) from a fitted model -----------

#' Simulate forward (or replay in-sample) from a fitted DINARCH model
#'
#' Unlike [simulate.dinarch_fit()] (a thin single-series wrapper around
#' [dinarch_simulate()]), this is id-aware and history-aware: it uses the
#' `data` you supply to seed lags and default covariates per id, so the
#' simplest call - `dinarch_project(fit, data, horizon = 12)` - "just
#' simulates forward with everything fixed" with no further arguments
#' needed. A `"dinarch_fit"` object is a fitted model only (no stored data),
#' so `data` is always explicit here - there's no ambiguity about which
#' data a projection is built from, and it's easy to project the same fit
#' forward from updated/extended data without refitting.
#'
#' @param object A `"dinarch_fit"` object, fit with [dinarch_fit_ml()]/
#'   [dinarch_fit_bayes()].
#' @param data A data.frame or data.table with the historical series to
#'   project from: at least `y` and `time` columns (and `id`, and any
#'   variables `object$formula` references, if used) - the same shape you'd
#'   pass to [dinarch_fit_ml()]/[dinarch_fit_bayes()]. Validated the same
#'   way: required columns present, no missing values, `y` a non-negative
#'   integer count series, `time` strictly increasing with a constant step
#'   size (no gaps or duplicates) within each `id`.
#' @param horizon `NULL` (default) replays **in-sample**: for each id,
#'   starting from its first `n_lags` real observations, `y` is simulated
#'   forward (unconditionally - the model's own simulated `y` feeds the
#'   lags, not the real observed values) over that id's actual
#'   historical time range and covariates, for a "does the fitted model
#'   look reasonable" check. A positive integer instead simulates
#'   **out-of-sample**: for each id, `horizon` periods starting at
#'   `time = <last observed> + 1`, seeded with that id's last
#'   `n_lags` real `y` values.
#' @param nsim Number of stochastic replicate simulations, matching
#'   [stats::simulate()]'s convention for "number of replicate response
#'   draws". Note [simulate.dinarch_fit()]'s `nsim` means something
#'   different there (horizon length) - don't assume the two match.
#' @param new_data Optional `data.frame`/`data.table` with a `time`
#'   column, an `id` column (required if the fit has more than one
#'   id, optional otherwise), and any of the raw covariate columns
#'   `object$formula` references, giving projected/alternative covariate
#'   values for specific `(id, time)` pairs. Any `(id, time)` in
#'   the simulated window not covered here falls back to the default: the
#'   last observed row of `data` held fixed (out-of-sample) or the real
#'   historical value (in-sample - lets you re-simulate history under an
#'   alternative covariate path, e.g. for counterfactual checks; `time` is
#'   *not* required to be after the last observed value for this reason).
#' @param new_para Optional `data.frame`/`data.table`, same `(id,
#'   time)` shape as `new_data`, with columns `phi`, `b1, ..., bk`
#'   (`k = object$n_lags`, always `b1` even at `n_lags = 1` - never bare
#'   `b`), and/or any of `names(object$coefficients$beta)` (e.g.
#'   `"(Intercept)"`, `"gdp"` - quote non-syntactic names when building
#'   the table, e.g. `` `(Intercept)` = ... `` or
#'   `` data.frame(check.names = FALSE, "(Intercept)" = ...) ``), giving
#'   `b`/`phi`/`beta` overrides for specific periods - e.g. built from
#'   [dinarch_transition_path()] to move `b` from one persistence regime
#'   to another over (most likely) the forecast horizon. Unspecified
#'   periods/parameters use the fitted point estimate (or that
#'   replicate's posterior draw, see `use_draw`). A column named the same
#'   as a raw covariate (e.g. `"gdp"`) means "override `gdp`'s
#'   *coefficient*" here, vs. "override `gdp`'s *value*" as a `new_data`
#'   column - the two tables are separate, so there's no ambiguity in
#'   code, just something to keep straight when naming columns.
#' @param use_draw For Bayesian fits only: `NULL` (default) uses the
#'   posterior mean for every replicate; `"random"` samples a fresh
#'   posterior draw per replicate, propagating full parameter uncertainty
#'   into the `nsim` replicates on top of the observation-level NB
#'   noise; an integer pins one specific draw for every replicate.
#' @param y_threshold Passed through to [dinarch_simulate()] for each
#'   id/replicate - see there. Applied identically across all ids
#'   and replicates (not id-specific). If you want a population-based
#'   limit, compute it yourself first (e.g. `fraction * population`) -
#'   there's no distinguished `population` concept here, matching the fit
#'   functions.
#' @param seed Optional random seed.
#' @param y,time,id Column names (strings) in `data`/`new_data`/`new_para` -
#'   same convention as [dinarch_fit_ml()]. `id` identifies separate
#'   series; if `NULL` (default), all rows of `data` are treated as one
#'   series.
#'
#' @return A `data.table` (also classed `"dinarch_project"`, so
#'   [plot.dinarch_project()] works, but otherwise behaves exactly like an
#'   ordinary `data.table` - nothing about the usual `data.table` syntax
#'   changes) with columns `id`, `time`, `sim` (replicate id, `1:nsim`),
#'   `pred` (`TRUE` for simulated/projected draws, `FALSE` for `data`'s
#'   real observations), `y`, and the raw covariate columns
#'   `object$formula` references (whatever values were actually used, so
#'   the output is self-documenting even when mixing defaults with
#'   `new_data` overrides). `data`'s historical rows are repeated once per
#'   `sim` replicate, so `out[out$sim == k, ]` is always one complete,
#'   continuously-time-indexed series (history through the simulated
#'   window) for every `id` - convenient for downstream plotting/analysis
#'   at the cost of `data` being duplicated `nsim` times in the result.
#'
#' @export
dinarch_project <- function(object, data,
                             horizon = NULL,
                             nsim = 1,
                             new_data = NULL,
                             new_para = NULL,
                             use_draw = NULL,
                             y_threshold = NULL,
                             seed = NULL,
                             y = "y", time = "time", id = NULL) {

  if (!inherits(object, "dinarch_fit")) {
    stop("`object` must be a \"dinarch_fit\" object.")
  }
  stopifnot(is.numeric(nsim), length(nsim) == 1, nsim >= 1, nsim == round(nsim))
  if (!is.null(horizon)) {
    stopifnot(is.numeric(horizon), length(horizon) == 1, horizon >= 1, horizon == round(horizon))
  }
  is_insample <- is.null(horizon)

  if (!is.null(seed)) set.seed(seed)

  formula  <- object$formula
  n_lags   <- object$n_lags
  cov_vars <- all.vars(formula)

  fit_data <- .build_dinarch_panel(data, y, time, id, cov_vars)
  fit_data <- fit_data[, !"orig_row", with = FALSE]

  ids <- unique(fit_data$id)
  multi_id <- length(ids) > 1

  .check_id_time_cols <- function(df, arg_name) {
    if (is.null(df)) return(invisible(NULL))
    if (!"time" %in% names(df)) {
      stop(sprintf("`%s` must have a `time` column.", arg_name))
    }
    if ("id" %in% names(df)) {
      bad <- setdiff(unique(df$id), ids)
      if (length(bad) > 0) {
        stop(sprintf(
          "`%s$id` contains value(s) not present in `data`: %s.",
          arg_name, paste(bad, collapse = ", ")
        ))
      }
    } else if (multi_id) {
      stop(sprintf(
        "`%s` must have an `id` column - this fit has more than one id.",
        arg_name
      ))
    }
    invisible(NULL)
  }
  .check_id_time_cols(new_data, "new_data")
  .check_id_time_cols(new_para, "new_para")

  new_data <- if (!is.null(new_data)) data.table::as.data.table(new_data) else NULL
  new_para <- if (!is.null(new_para)) data.table::as.data.table(new_para) else NULL
  b_cols <- paste0("b", seq_len(n_lags))
  beta_names <- names(object$coefficients$beta)
  n_beta <- length(beta_names)

  any_data_matched <- is.null(new_data)
  any_para_matched <- is.null(new_para)

  # --- per-id scaffolding: y_init, time sequence, default (+ new_data-
  # overridden) covariates, matched new_para overrides, and the full
  # historical block (to prepend, once per `sim`, to that id's output) -
  # none of this depends on which parameter draw a given replicate uses,
  # so it's built once and reused across all `nsim` replicates. -----------
  scaffold <- lapply(ids, function(g) {
    g_data <- fit_data[fit_data$id == g]

    if (is_insample) {
      if (nrow(g_data) <= n_lags) {
        warning(sprintf(
          "Id %s has %d observation(s), not enough for n_lags = %d; skipping.",
          g, nrow(g_data), n_lags
        ))
        return(NULL)
      }
      y_init <- as.numeric(utils::head(g_data$y, n_lags))
      sim_rows <- seq(n_lags + 1L, nrow(g_data))
      idx_seq <- g_data$time[sim_rows]
      cov_df <- data.table::copy(g_data[sim_rows, cov_vars, with = FALSE])
    } else {
      y_init <- as.numeric(utils::tail(g_data$y, n_lags))
      idx_seq <- max(g_data$time) + seq_len(horizon)
      last_row <- g_data[nrow(g_data)]
      cov_df <- data.table::copy(last_row[rep(1L, horizon), cov_vars, with = FALSE])
    }

    if (!is.null(new_data)) {
      nd_g <- if ("id" %in% names(new_data)) new_data[new_data$id == g] else new_data
      pos <- match(nd_g$time, idx_seq)
      matched <- !is.na(pos)
      if (any(matched)) {
        any_data_matched <<- TRUE
        for (col in intersect(cov_vars, names(nd_g))) {
          cov_df[pos[matched], (col) := nd_g[[col]][matched]]
        }
      }
    }

    np_match <- NULL
    if (!is.null(new_para)) {
      np_g <- if ("id" %in% names(new_para)) new_para[new_para$id == g] else new_para
      pos <- match(np_g$time, idx_seq)
      matched <- !is.na(pos)
      if (any(matched)) {
        any_para_matched <<- TRUE
        np_match <- list(pos = pos[matched], data = np_g[matched])
      }
    }

    list(id = g, y_init = y_init, idx_seq = idx_seq, cov_df = cov_df,
         n_sim_periods = length(idx_seq), np_match = np_match, hist = g_data)
  })
  scaffold <- Filter(Negate(is.null), scaffold)

  if (length(scaffold) == 0) {
    stop("No id had enough observations to simulate (see warnings above).")
  }
  if (!is.null(new_data) && !any_data_matched) {
    warning("`new_data` did not match any (id, time) pair in the simulated window; using defaults for all covariates.")
  }
  if (!is.null(new_para) && !any_para_matched) {
    warning("`new_para` did not match any (id, time) pair in the simulated window; using fitted parameters for all periods.")
  }

  # --- simulate: nsim replicates x ids, each prepended with the full
  # historical block for that id (marked pred = FALSE) so every
  # `sim`'s rows form one complete, continuously-indexed series. -----------
  replicate_results <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    pars <- .dinarch_pick_params(object, use_draw)

    id_sims <- lapply(scaffold, function(sc) {
      n_p <- sc$n_sim_periods
      b_path <- .as_coef_matrix(pars$b, n_p, n_lags)
      phi_path <- rep(pars$phi, n_p)
      beta_path <- .as_coef_matrix(pars$beta, n_p, n_beta)

      if (!is.null(sc$np_match)) {
        pos <- sc$np_match$pos
        np_data <- sc$np_match$data
        if ("phi" %in% names(np_data)) phi_path[pos] <- np_data$phi
        for (col in intersect(b_cols, names(np_data))) {
          j <- as.integer(sub("^b", "", col))
          b_path[pos, j] <- np_data[[col]]
        }
        for (col in intersect(beta_names, names(np_data))) {
          j <- match(col, beta_names)
          beta_path[pos, j] <- np_data[[col]]
        }
      }

      sim <- dinarch_simulate(
        n = n_p, b = b_path, phi = phi_path,
        formula = formula,
        covariates = if (length(cov_vars) > 0) sc$cov_df else NULL,
        beta = beta_path,
        y_init = sc$y_init,
        y_threshold = y_threshold
      )
      sim$time <- sc$idx_seq
      sim$id <- sc$id
      if (length(cov_vars) > 0) sim <- cbind(sim, sc$cov_df)
      sim$pred <- TRUE

      hist_block <- data.table::copy(sc$hist)
      hist_block[, "pred" := FALSE]

      data.table::rbindlist(list(hist_block, sim), use.names = TRUE)
    })

    combined <- data.table::rbindlist(id_sims)
    combined[, "sim" := s]
    replicate_results[[s]] <- combined
  }

  out <- data.table::rbindlist(replicate_results)
  data.table::setcolorder(out, c("id", "time", "sim", "pred", "y"))
  out <- .new_dinarch_project(out)
  out[]
}

#' Attach the "dinarch_project" class to a data.table (no validation;
#' idempotent - the class is never duplicated)
#' @noRd
.new_dinarch_project <- function(x) {
  data.table::setattr(x, "class", c("dinarch_project", setdiff(class(x), "dinarch_project")))
  x
}

#' Coerce to a DINARCH projection
#'
#' Validates that `x` has the shape [dinarch_project()]'s output has, and
#' returns it classed `"dinarch_project"` so `plot()`/`summary()` dispatch
#' on it. The main use is combining per-series projections (each from its
#' own fit and data) into one table: a plain `rbind()`/
#' `data.table::rbindlist()` drops the `"dinarch_project"` class, while
#' `as_dinarch_project()` takes the list directly, row-binds it
#' (`data.table::rbindlist(use.names = TRUE, fill = TRUE)` - columns are
#' matched by name, and a covariate present in only some projections is
#' filled with `NA` elsewhere), and reattaches the class.
#'
#' Checks that the result has columns `id`, `time`, `sim`, `pred`
#' and `y`, none containing missing values; that `pred` is logical,
#' `sim` whole numbers, and `y` non-negative integers; and that no
#' `(id, time, sim)` combination appears twice - which is what happens when
#' two combined projections share an `id` label (e.g. both left at the
#' single-id default `0`), and would otherwise silently merge the two
#' series in `plot()`/`summary()`. The input is never modified in place;
#' a copy is returned.
#'
#' @param x A `data.frame`/`data.table` shaped like [dinarch_project()]'s
#'   output, or a list of them.
#' @return A `data.table`, also classed `"dinarch_project"`.
#' @seealso [dinarch_project()]
#' @examples
#' sims <- lapply(1:2, function(i) {
#'   dinarch_simulate(n = 40, b = 0.3, phi = 8, beta = log(4), seed = i)
#' })
#' dats <- lapply(seq_along(sims), function(i) {
#'   data.frame(y = sims[[i]]$y, time = sims[[i]]$time, id = LETTERS[i])
#' })
#' projs <- lapply(dats, function(d) {
#'   fit <- dinarch_fit_ml(d, n_lags = 1, vcov = FALSE)
#'   dinarch_project(fit, d, id = "id", horizon = 5, nsim = 10, seed = 1)
#' })
#' combined <- as_dinarch_project(projs)
#' summary(combined)
#' @export
as_dinarch_project <- function(x) {
  if (is.data.frame(x)) {
    x <- data.table::as.data.table(x)
  } else if (is.list(x) && all(vapply(x, is.data.frame, logical(1)))) {
    x <- data.table::rbindlist(x, use.names = TRUE, fill = TRUE)
  } else {
    stop("`x` must be a data.frame/data.table, or a list of them.")
  }

  required_cols <- c("id", "time", "sim", "pred", "y")
  .check_columns(x, required_cols, arg_name = "x")
  for (col in required_cols) .check_no_na(x[[col]], col)

  if (!is.logical(x$pred)) {
    stop("`pred` must be logical (TRUE for simulated draws, FALSE for observed rows).")
  }
  if (!is.numeric(x$sim) || any(x$sim != round(x$sim))) {
    stop("`sim` must be whole numbers (replicate ids).")
  }
  if (!is.numeric(x$y)) {
    stop("`y` must be non-negative integers (a count series).")
  }
  .check_nonneg_integer(x$y, "y")
  if (anyDuplicated(x, by = c("id", "time", "sim")) > 0) {
    stop("`x` has duplicated (id, time, sim) rows - combined projections must each use distinct `id` values.")
  }

  .new_dinarch_project(x)[]
}

# Expand one id's sorted (time, lo, hi, med) rows into a step ("post",
# i.e. flat-then-jump - matching geom_step()'s default direction) polygon
# boundary for geom_ribbon(), which has no built-in step variant of its own.
.dinarch_stepify_ribbon <- function(dt) {
  data.table::setorderv(dt, "time")
  if (nrow(dt) <= 1) return(dt[, c("time", "lo", "hi", "med"), with = FALSE])
  idx <- c(dt$time[1], rep(dt$time[-1], each = 2))
  step_col <- function(v) c(rep(utils::head(v, -1), each = 2), utils::tail(v, 1))
  data.table::data.table(
    time = idx, lo = step_col(dt$lo), hi = step_col(dt$hi), med = step_col(dt$med)
  )
}

#' Plot a DINARCH projection
#'
#' A minimal `ggplot2` plot of a [dinarch_project()] result: the historical
#' `y` (`x$pred == FALSE`, deduplicated across replicates) as a single
#' line, plus either every simulated replicate as its own thin line
#' (`type = "spaghetti"`) or a shaded quantile ribbon with a median line
#' (`type = "ribbon"`, the default) summarizing the replicates at each
#' `time`. Each id's projected series is prepended with the last historical
#' point strictly before it, so the projection visually starts exactly
#' where the historical line ends instead of jumping in from a
#' disconnected point. One facet per id (skipped for a single-id/ungrouped
#' fit). Returns the `ggplot` object rather than just drawing it, so it
#' prints like a normal plot when called at the top level but can also be
#' captured and extended, e.g. `plot(proj) + ggplot2::ggtitle("Forecast")` -
#' titles, themes, and any other customization are deliberately left to the
#' caller rather than built in here.
#'
#' @param x A `"dinarch_project"` object, i.e. the output of
#'   [dinarch_project()].
#' @param type `"ribbon"` (default) or `"spaghetti"` - see above.
#' @param probs For `type = "ribbon"`, the lower/upper quantile
#'   probabilities of the shaded band (default the 10th/90th percentile).
#' @param step If `TRUE`, draw step functions (`ggplot2::geom_step()`,
#'   holding each value flat until the next `time`) instead of straight
#'   line segments between points - appropriate for count data, where a
#'   straight line between two integers implies values that were never
#'   observed. Default `FALSE`.
#' @param ... Unused (kept for `plot()` generic compatibility).
#'
#' @return A `ggplot` object.
#' @export
plot.dinarch_project <- function(object, type = c("ribbon", "spaghetti"), probs = c(0.1, 0.9),
                                  step = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.dinarch_project() but is not installed.")
  }

  stopifnot("dinarch_project" %in% class(object))

  type <- match.arg(type)
  stopifnot(length(probs) == 2, probs[1] < probs[2])
  stopifnot(is.logical(step), length(step) == 1, !is.na(step))

  line_geom <- if (step) ggplot2::geom_step else ggplot2::geom_line

  multi_id <- length(unique(object$id)) > 1

  # historical rows are identical across every `sim` replicate - dedupe to
  # the single historical line rather than drawing nsim overlapping copies.
  history <- unique(object[object$pred == FALSE, c("id", "time", "y"), with = FALSE])
  data.table::setorderv(history, c("id", "time"))
  proj <- object[object$pred == TRUE]

  # Last historical point strictly before each id's first projected
  # time - prepended below so the projection connects to, rather than
  # jumps in disconnected from, the end of the historical line.
  connector <- NULL
  if (nrow(history) > 0 && nrow(proj) > 0) {
    connector <- data.table::rbindlist(lapply(unique(proj$id), function(g) {
      first_idx <- min(proj$time[proj$id == g])
      h_g <- history[history$id == g & history$time < first_idx]
      if (nrow(h_g) == 0) return(NULL)
      h_g[which.max(h_g$time)]
    }))
  }
  has_connector <- !is.null(connector) && nrow(connector) > 0

  p <- ggplot2::ggplot()

  if (nrow(history) > 0) {
    p <- p + line_geom(
      data = history, ggplot2::aes(x = .data$time, y = .data$y)
    )
  }

  if (type == "spaghetti") {
    spaghetti_data <- proj[, c("id", "sim", "time", "y"), with = FALSE]
    if (has_connector) {
      sims <- unique(proj$sim)
      rep_pos <- rep(seq_len(nrow(connector)), each = length(sims))
      conn <- data.table::data.table(
        id = connector$id[rep_pos], sim = rep(sims, times = nrow(connector)),
        time = connector$time[rep_pos], y = connector$y[rep_pos]
      )
      spaghetti_data <- data.table::rbindlist(list(conn, spaghetti_data), use.names = TRUE)
    }
    data.table::setorderv(spaghetti_data, c("id", "sim", "time"))
    a <- .alpha_est(spaghetti_data, group = c("id", "sim"), x = "time")

    p <- p + line_geom(
      data = spaghetti_data,
      ggplot2::aes(x = .data$time, y = .data$y, group = interaction(.data$id, .data$sim)),
      colour = "steelblue", alpha = a
    )
  } else {
    ribbon_stats <- function(v) list(
      lo  = stats::quantile(v, probs[1]),
      hi  = stats::quantile(v, probs[2]),
      med = stats::median(v)
    )
    agg <- proj[, ribbon_stats(.SD[["y"]]), by = c("id", "time"), .SDcols = "y"]

    if (has_connector) {
      conn <- data.table::data.table(
        id = connector$id, time = connector$time,
        lo = connector$y, hi = connector$y, med = connector$y
      )
      agg <- data.table::rbindlist(list(conn, agg), use.names = TRUE)
    }
    data.table::setorderv(agg, c("id", "time"))

    ribbon_data <- if (step) agg[, .dinarch_stepify_ribbon(.SD), by = "id"] else agg

    p <- p +
      ggplot2::geom_ribbon(
        data = ribbon_data, ggplot2::aes(x = .data$time, ymin = .data$lo, ymax = .data$hi),
        fill = "steelblue", alpha = 0.3
      ) +
      line_geom(
        data = agg, ggplot2::aes(x = .data$time, y = .data$med),
        colour = "steelblue"
      )
  }

  if (multi_id) {
    p <- p + ggplot2::facet_wrap(~id, scales = "free_y")
  }

  p + ggplot2::labs(x = "time", y = "y")
}

#' Consecutive block index for `n_t` ordered periods split into `nsplit`
#' near-equal blocks (later blocks get the extra period); `NULL` = one block
#' per period. `nsplit` is capped at `n_t`.
#' @noRd
.dinarch_time_blocks <- function(n_t, nsplit) {
  if (is.null(nsplit)) return(seq_len(n_t))
  k <- min(nsplit, n_t)
  as.integer(ceiling(seq_len(n_t) * k / n_t))
}

#' Summarize a DINARCH projection/prediction
#'
#' Collapses a [dinarch_project()]/[predict.dinarch_fit()] result into the
#' mean and simulation quantiles of `y` over blocks of consecutive periods -
#' a compact numeric alternative to [plot()]. Simulated rows (`pred == TRUE`)
#' pool all `nsim` replicates within each block; historical rows
#' (`pred == FALSE`, real data repeated once per replicate) summarize the
#' observed series itself. Blocks are formed separately for each id, so ids
#' with different time ranges are split independently.
#'
#' @param object A `"dinarch_project"` object.
#' @param probs Quantiles of `y` to report (default `c(0.025, 0.5, 0.975)`),
#'   returned as columns `q2.5`, `q50`, `q97.5`, etc.
#' @param nsplit Number of consecutive blocks each id's simulated window is
#'   split into (default 1: the whole window; capped at the number of
#'   simulated periods). History is always one block per id. `NULL` gives one
#'   row per period, for history and simulated rows alike.
#' @param include_hist Whether to include rows summarizing the historical
#'   (`pred == FALSE`) data.
#' @param ... Unused.
#'
#' @return An object of class `"summary.dinarch_project"`: a list whose
#'   `table` element is a `data.table` with columns `id`, `pred`, `start`,
#'   `end` (inclusive first/last `time` of the block), `mean`, and one
#'   `q<percent>` column per `probs` entry, sorted by `id`, `pred`, `start`.
#' @export
summary.dinarch_project <- function(object, probs = c(0.025, 0.5, 0.975), nsplit = 1,
                                    include_hist = TRUE, ...) {
  stopifnot("dinarch_project" %in% class(object))
  stopifnot(
    is.numeric(probs), length(probs) >= 1,
    all(probs > 0 & probs < 1)
  )
  probs <- sort(unique(probs))
  if (!is.null(nsplit)) {
    stopifnot(is.numeric(nsplit), length(nsplit) == 1, nsplit >= 1, nsplit == round(nsplit))
  }
  stopifnot(is.logical(include_hist), length(include_hist) == 1, !is.na(include_hist))

  dt <- data.table::as.data.table(object)
  if (!any(dt$pred)) stop("`object` has no simulated rows (pred == TRUE) to summarize.")
  if (!include_hist) dt <- dt[dt$pred == TRUE, ]

  pct <- probs * 100
  pct_chr <- sub("\\.?0+$", "", formatC(pct, format = "f", digits = 6))
  qnames <- paste0("q", pct_chr)
  summarise_y <- function(y) {
    q <- stats::quantile(y, probs = probs, names = FALSE)
    c(list(mean = mean(y)), stats::setNames(as.list(q), qnames))
  }

  # One block index per distinct (id, pred, time), computed within each id so
  # ids with different time ranges are split independently.
  keys <- unique(dt[, c("id", "pred", "time")])
  data.table::setorderv(keys, c("id", "pred", "time"))
  pred_keys <- keys[keys$pred == TRUE, ]
  if (!is.null(nsplit) && nsplit > 1 && any(pred_keys[, .N, by = "id"][["N"]] <= nsplit)) {
    message("nsplit >= the number of simulated periods for at least one id; summarizing each of those periods separately.")
  }
  pred_keys[, "block" := .dinarch_time_blocks(.N, nsplit), by = "id"]
  hist_keys <- keys[keys$pred == FALSE, ]
  if (nrow(hist_keys) > 0L) {
    hist_keys[, "block" := .dinarch_time_blocks(.N, if (is.null(nsplit)) NULL else 1L), by = "id"]
    keys <- data.table::rbindlist(list(hist_keys, pred_keys), use.names = TRUE)
  } else {
    keys <- pred_keys
  }

  dt <- merge(dt, keys, by = c("id", "pred", "time"))
  agg <- dt[, c(
    list(start = min(.SD[["time"]]), end = max(.SD[["time"]])),
    summarise_y(.SD[["y"]])
  ), by = c("id", "pred", "block"), .SDcols = c("time", "y")]
  agg[, "block" := NULL]
  data.table::setorderv(agg, c("id", "pred", "start"))

  structure(
    list(table = agg, probs = probs, nsplit = nsplit, include_hist = include_hist),
    class = "summary.dinarch_project"
  )
}

#' @export
print.summary.dinarch_project <- function(x, ...) {
  cat("DINARCH projection\n")
  cat(strrep("-", 50), "\n", sep = "")
  print(as.data.frame(x$table), row.names = FALSE, digits = 4)
  invisible(x)
}
