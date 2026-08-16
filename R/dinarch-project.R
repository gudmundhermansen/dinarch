# Projection (forecast / in-sample replay) from a fitted model -----------

#' Simulate forward (or replay in-sample) from a fitted DINARCH model
#'
#' Unlike [simulate.dinarch_fit()] (a thin single-series wrapper around
#' [dinarch_simulate()]), this is group-aware and history-aware: it uses
#' `object$data` (stored by [dinarch_fit_ml()]/[dinarch_fit_bayes()]) to
#' seed lags and default covariates per group, so the simplest call -
#' `dinarch_project(fit, horizon = 12)` - "just simulates forward with
#' everything fixed" with no further arguments needed.
#'
#' @param object A `"dinarch_fit"` object, fit with a version of
#'   [dinarch_fit_ml()]/[dinarch_fit_bayes()] that stores `$data`.
#' @param horizon `NULL` (default) replays **in-sample**: for each group,
#'   starting from its first `n_lags` real observations, `y` is simulated
#'   forward (unconditionally - the model's own simulated `y` feeds the
#'   lags, not the real observed values) over that group's actual
#'   historical index range and covariates, for a "does the fitted model
#'   look reasonable" check. A positive integer instead simulates
#'   **out-of-sample**: for each group, `horizon` periods starting at
#'   `index = <last observed> + 1`, seeded with that group's last
#'   `n_lags` real `y` values.
#' @param nsim Number of stochastic replicate simulations, matching
#'   [stats::simulate()]'s convention for "number of replicate response
#'   draws". Note [simulate.dinarch_fit()]'s `nsim` means something
#'   different there (horizon length) - don't assume the two match.
#' @param new_data Optional `data.frame`/`data.table` with an `index`
#'   column, a `group` column (required if the fit has more than one
#'   group, optional otherwise), and any of the raw covariate columns
#'   `object$formula` references, giving projected/alternative covariate
#'   values for specific `(group, index)` pairs. Any `(group, index)` in
#'   the simulated window not covered here falls back to the default: the
#'   last observed row held fixed (out-of-sample) or the real historical
#'   value (in-sample - lets you re-simulate history under an alternative
#'   covariate path, e.g. for counterfactual checks; `index` is *not*
#'   required to be after the last observed value for this reason).
#' @param new_para Optional `data.frame`/`data.table`, same `(group,
#'   index)` shape as `new_data`, with columns `phi`, `b1, ..., bk`
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
#'   group/replicate - see there. Applied identically across all groups
#'   and replicates (not group-specific). If you want a population-based
#'   limit, compute it yourself first (e.g. `fraction * population`) -
#'   there's no distinguished `population` concept here, matching the fit
#'   functions.
#' @param seed Optional random seed.
#'
#' @return A `data.table` (also classed `"dinarch_project"`, so
#'   [plot.dinarch_project()] works, but otherwise behaves exactly like an
#'   ordinary `data.table` - nothing about the usual `data.table` syntax
#'   changes) with columns `group`, `index`, `sim` (replicate id,
#'   `1:nsim`), `y`, and the raw covariate columns `object$formula`
#'   references (whatever values were actually used, so the output is
#'   self-documenting even when mixing defaults with `new_data`
#'   overrides).
#'
#' @export
dinarch_project <- function(object,
                             horizon = NULL,
                             nsim = 1,
                             new_data = NULL,
                             new_para = NULL,
                             use_draw = NULL,
                             y_threshold = NULL,
                             seed = NULL) {

  if (!inherits(object, "dinarch_fit")) {
    stop("`object` must be a \"dinarch_fit\" object.")
  }
  if (is.null(object$data)) {
    stop(
      "This fit object has no stored `data` (fit with an older version of ",
      "the package?) - refit with the current dinarch_fit_ml()/",
      "dinarch_fit_bayes() to use dinarch_project()."
    )
  }
  stopifnot(is.numeric(nsim), length(nsim) == 1, nsim >= 1, nsim == round(nsim))
  if (!is.null(horizon)) {
    stopifnot(is.numeric(horizon), length(horizon) == 1, horizon >= 1, horizon == round(horizon))
  }
  is_insample <- is.null(horizon)

  if (!is.null(seed)) set.seed(seed)

  fit_data <- object$data
  formula  <- object$formula
  n_lags   <- object$n_lags
  cov_vars <- all.vars(formula)
  groups   <- unique(fit_data$group)
  multi_group <- length(groups) > 1

  .check_group_index_cols <- function(df, arg_name) {
    if (is.null(df)) return(invisible(NULL))
    if (!"index" %in% names(df)) {
      stop(sprintf("`%s` must have an `index` column.", arg_name))
    }
    if ("group" %in% names(df)) {
      bad <- setdiff(unique(df$group), groups)
      if (length(bad) > 0) {
        stop(sprintf(
          "`%s$group` contains value(s) not present in the fitted data: %s.",
          arg_name, paste(bad, collapse = ", ")
        ))
      }
    } else if (multi_group) {
      stop(sprintf(
        "`%s` must have a `group` column - this fit has more than one group.",
        arg_name
      ))
    }
    invisible(NULL)
  }
  .check_group_index_cols(new_data, "new_data")
  .check_group_index_cols(new_para, "new_para")

  new_data <- if (!is.null(new_data)) data.table::as.data.table(new_data) else NULL
  new_para <- if (!is.null(new_para)) data.table::as.data.table(new_para) else NULL
  b_cols <- paste0("b", seq_len(n_lags))
  beta_names <- names(object$coefficients$beta)
  n_beta <- length(beta_names)

  any_data_matched <- is.null(new_data)
  any_para_matched <- is.null(new_para)

  # --- per-group scaffolding: y_init, index sequence, default (+ new_data-
  # overridden) covariates, matched new_para overrides - none of this
  # depends on which parameter draw a given replicate uses, so it's built
  # once and reused across all `nsim` replicates. ---------------------------
  scaffold <- lapply(groups, function(g) {
    g_data <- fit_data[fit_data$group == g]
    data.table::setorderv(g_data, "index")

    if (is_insample) {
      if (nrow(g_data) <= n_lags) {
        warning(sprintf(
          "Group %s has %d observation(s), not enough for n_lags = %d; skipping.",
          g, nrow(g_data), n_lags
        ))
        return(NULL)
      }
      y_init <- as.numeric(utils::head(g_data$y, n_lags))
      sim_rows <- seq(n_lags + 1L, nrow(g_data))
      idx_seq <- g_data$index[sim_rows]
      cov_df <- data.table::copy(g_data[sim_rows, cov_vars, with = FALSE])
    } else {
      y_init <- as.numeric(utils::tail(g_data$y, n_lags))
      idx_seq <- max(g_data$index) + seq_len(horizon)
      last_row <- g_data[nrow(g_data)]
      cov_df <- data.table::copy(last_row[rep(1L, horizon), cov_vars, with = FALSE])
    }

    if (!is.null(new_data)) {
      nd_g <- if ("group" %in% names(new_data)) new_data[new_data$group == g] else new_data
      pos <- match(nd_g$index, idx_seq)
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
      np_g <- if ("group" %in% names(new_para)) new_para[new_para$group == g] else new_para
      pos <- match(np_g$index, idx_seq)
      matched <- !is.na(pos)
      if (any(matched)) {
        any_para_matched <<- TRUE
        np_match <- list(pos = pos[matched], data = np_g[matched])
      }
    }

    list(group = g, y_init = y_init, idx_seq = idx_seq, cov_df = cov_df,
         n_sim_periods = length(idx_seq), np_match = np_match)
  })
  scaffold <- Filter(Negate(is.null), scaffold)

  if (length(scaffold) == 0) {
    stop("No group had enough observations to simulate (see warnings above).")
  }
  if (!is.null(new_data) && !any_data_matched) {
    warning("`new_data` did not match any (group, index) pair in the simulated window; using defaults for all covariates.")
  }
  if (!is.null(new_para) && !any_para_matched) {
    warning("`new_para` did not match any (group, index) pair in the simulated window; using fitted parameters for all periods.")
  }

  # --- simulate: nsim replicates x groups -----------------------------------
  replicate_results <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    pars <- .dinarch_pick_params(object, use_draw)

    group_sims <- lapply(scaffold, function(sc) {
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
      sim$index <- sc$idx_seq
      sim$group <- sc$group
      if (length(cov_vars) > 0) sim <- cbind(sim, sc$cov_df)
      sim
    })

    combined <- data.table::rbindlist(group_sims)
    combined[, "sim" := s]
    replicate_results[[s]] <- combined
  }

  out <- data.table::rbindlist(replicate_results)
  data.table::setcolorder(out, c("group", "index", "sim", "y"))

  # Carried along for plot.dinarch_project() so `plot(proj)` works with no
  # further arguments - the fit's own historical (group, index, y), not
  # recomputed here since `fit_data` already has it.
  data.table::setattr(out, "history", fit_data[, c("group", "index", "y"), with = FALSE])
  data.table::setattr(out, "class", c("dinarch_project", class(out)))
  out[]
}

# Expand one group's sorted (index, lo, hi, med) rows into a step ("post",
# i.e. flat-then-jump - matching geom_step()'s default direction) polygon
# boundary for geom_ribbon(), which has no built-in step variant of its own.
.dinarch_stepify_ribbon <- function(dt) {
  data.table::setorderv(dt, "index")
  if (nrow(dt) <= 1) return(dt[, c("index", "lo", "hi", "med"), with = FALSE])
  idx <- c(dt$index[1], rep(dt$index[-1], each = 2))
  step_col <- function(v) c(rep(utils::head(v, -1), each = 2), utils::tail(v, 1))
  data.table::data.table(
    index = idx, lo = step_col(dt$lo), hi = step_col(dt$hi), med = step_col(dt$med)
  )
}

#' Plot a DINARCH projection
#'
#' A minimal `ggplot2` plot of a [dinarch_project()] result: the fitted
#' model's historical `y` (if available) as a single line, plus either
#' every simulated replicate as its own thin line (`type = "spaghetti"`)
#' or a shaded quantile ribbon with a median line (`type = "ribbon"`,
#' the default) summarizing the replicates at each `index`. Each group's
#' projected series is prepended with the last historical point strictly
#' before it, so the projection visually starts exactly where the
#' historical line ends instead of jumping in from a disconnected point.
#' One facet per group (skipped for a single-group/ungrouped fit). Returns
#' the `ggplot` object rather than just drawing it, so it prints like a
#' normal plot when called at the top level but can also be captured and
#' extended, e.g. `plot(proj) + ggplot2::ggtitle("Forecast")` - titles,
#' themes, and any other customization are deliberately left to the
#' caller rather than built in here.
#'
#' @param x A `"dinarch_project"` object, i.e. the output of
#'   [dinarch_project()].
#' @param type `"ribbon"` (default) or `"spaghetti"` - see above.
#' @param probs For `type = "ribbon"`, the lower/upper quantile
#'   probabilities of the shaded band (default the 10th/90th percentile).
#' @param step If `TRUE`, draw step functions (`ggplot2::geom_step()`,
#'   holding each value flat until the next `index`) instead of straight
#'   line segments between points - appropriate for count data, where a
#'   straight line between two integers implies values that were never
#'   observed. Default `FALSE`.
#' @param ... Unused (kept for `plot()` generic compatibility).
#'
#' @return A `ggplot` object.
#' @export
plot.dinarch_project <- function(x, type = c("ribbon", "spaghetti"), probs = c(0.1, 0.9),
                                  step = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.dinarch_project() but is not installed.")
  }
  type <- match.arg(type)
  stopifnot(length(probs) == 2, probs[1] < probs[2])
  stopifnot(is.logical(step), length(step) == 1, !is.na(step))

  line_geom <- if (step) ggplot2::geom_step else ggplot2::geom_line

  history <- attr(x, "history")
  multi_group <- length(unique(c(x$group, history$group))) > 1

  # Last historical point strictly before each group's first projected
  # index - prepended below so the projection connects to, rather than
  # jumps in disconnected from, the end of the historical line.
  connector <- NULL
  if (!is.null(history) && nrow(history) > 0) {
    connector <- data.table::rbindlist(lapply(unique(x$group), function(g) {
      first_idx <- min(x$index[x$group == g])
      h_g <- history[history$group == g & history$index < first_idx]
      if (nrow(h_g) == 0) return(NULL)
      h_g[which.max(h_g$index)]
    }))
  }
  has_connector <- !is.null(connector) && nrow(connector) > 0

  p <- ggplot2::ggplot()

  if (!is.null(history) && nrow(history) > 0) {
    p <- p + line_geom(
      data = history, ggplot2::aes(x = .data$index, y = .data$y)
    )
  }

  if (type == "spaghetti") {
    spaghetti_data <- x[, c("group", "sim", "index", "y"), with = FALSE]
    if (has_connector) {
      sims <- unique(x$sim)
      rep_pos <- rep(seq_len(nrow(connector)), each = length(sims))
      conn <- data.table::data.table(
        group = connector$group[rep_pos], sim = rep(sims, times = nrow(connector)),
        index = connector$index[rep_pos], y = connector$y[rep_pos]
      )
      spaghetti_data <- data.table::rbindlist(list(conn, spaghetti_data), use.names = TRUE)
    }
    data.table::setorderv(spaghetti_data, c("group", "sim", "index"))

    p <- p + line_geom(
      data = spaghetti_data,
      ggplot2::aes(x = .data$index, y = .data$y, group = interaction(.data$group, .data$sim)),
      colour = "steelblue", alpha = 0.15
    )
  } else {
    agg <- x[, list(
      lo  = stats::quantile(y, probs[1]),
      hi  = stats::quantile(y, probs[2]),
      med = stats::median(y)
    ), by = c("group", "index")]

    if (has_connector) {
      conn <- data.table::data.table(
        group = connector$group, index = connector$index,
        lo = connector$y, hi = connector$y, med = connector$y
      )
      agg <- data.table::rbindlist(list(conn, agg), use.names = TRUE)
    }
    data.table::setorderv(agg, c("group", "index"))

    ribbon_data <- if (step) agg[, .dinarch_stepify_ribbon(.SD), by = "group"] else agg

    p <- p +
      ggplot2::geom_ribbon(
        data = ribbon_data, ggplot2::aes(x = .data$index, ymin = .data$lo, ymax = .data$hi),
        fill = "steelblue", alpha = 0.3
      ) +
      line_geom(
        data = agg, ggplot2::aes(x = .data$index, y = .data$med),
        colour = "steelblue"
      )
  }

  if (multi_group) {
    p <- p + ggplot2::facet_wrap(~group, scales = "free_y")
  }

  p + ggplot2::labs(x = "index", y = "y")
}

#' Summarize a DINARCH projection/prediction
#'
#' Collapses a [dinarch_project()]/[predict.dinarch_fit()] result's `nsim`
#' replicates at each `(group, index)` into a mean and a simulation-based
#' prediction interval - a compact numeric alternative to [plot()]. `index`
#' carries its usual meaning: for an out-of-sample call (`horizon` set),
#' it's periods ahead of the last observed period; for an in-sample call,
#' it's the historical period being replayed.
#'
#' @param object A `"dinarch_project"` object.
#' @param level Prediction interval level (default 0.95, i.e. the
#'   simulated draws' 2.5th/97.5th percentiles at each period).
#' @param ... Unused.
#'
#' @return An object of class `"summary.dinarch_project"`, printing as a
#'   table with columns `group`, `index`, `mean`, `lower`, `upper`.
#' @export
summary.dinarch_project <- function(object, level = 0.95, ...) {
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1)
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  agg <- object[, list(
    mean  = mean(y),
    lower = stats::quantile(y, probs[1]),
    upper = stats::quantile(y, probs[2])
  ), by = c("group", "index")]
  data.table::setorderv(agg, c("group", "index"))

  structure(list(table = agg, level = level), class = "summary.dinarch_project")
}

#' @export
print.summary.dinarch_project <- function(x, ...) {
  cat(sprintf("DINARCH projection - %.0f%% simulation-based prediction interval\n", 100 * x$level))
  cat(strrep("-", 50), "\n", sep = "")
  print(as.data.frame(x$table), row.names = FALSE, digits = 4)
  invisible(x)
}
