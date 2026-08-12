# ML estimation -----------------------------------------------------------

#' Fit a DINARCH model by maximum likelihood
#'
#' Fits \deqn{\mu_t = \exp(a) \cdot population_t + \sum_{j=1}^k b_j y_{t-j}
#'   + \sum_l \exp(\beta_l) covariate_{t,l}}, \deqn{y_t \sim
#'   \mathrm{NegBinom}(r_t = \phi \mu_t,\; p_t = \phi/(1+\phi))}
#' by maximizing the (optionally penalized) log-likelihood via
#' [stats::optim()].
#'
#' @param data A data.frame or data.table containing at least the `y` and
#'   `index` columns (and `group`/`population`/`covariates` columns if
#'   used).
#' @param y,index,group,population Column names (strings) in `data`. `group`
#'   identifies separate series (lags are never taken across a group
#'   boundary); if `NULL`, all rows are treated as one series. `population`
#'   is optional; if `NULL` the baseline term reduces to `exp(a)`.
#' @param covariates Optional character vector of additional covariate
#'   column names in `data`.
#' @param k Number of autoregressive lags.
#' @param lambda Optional regularization strength discouraging `b` near its
#'   upper boundary and `phi` near 0. `NULL` (default) disables it.
#' @param start Optional custom starting values for the optimizer, on the
#'   internal unconstrained scale. If `NULL` (default), data-driven
#'   starting values are used: `a` is initialized from the observed
#'   mean(y)/mean(population) ratio (rather than a fixed value), which
#'   matters a lot when `population` is on a very different scale than
#'   the count series - a mismatched fixed start is a common cause of
#'   Nelder-Mead's simplex degenerating (`optim()` convergence code 10).
#' @param optimizer `"nelder-mead+bfgs"` (default) runs Nelder-Mead then
#'   refines with a BFGS pass; `"nelder-mead"` skips the refinement.
#' @param control Optional list passed to [stats::optim()]'s `control`
#'   argument (e.g. `list(maxit = 5000, reltol = 1e-10)`), merged over the
#'   package's defaults - use this to tune convergence behavior if the
#'   defaults struggle on a particular dataset.
#' @param max_restarts If Nelder-Mead reports a non-zero convergence code
#'   (most commonly code 10, simplex degeneracy), automatically restart
#'   from the last result up to this many times - the standard remedy for
#'   this failure mode, since a degenerate simplex doesn't mean the
#'   optimum itself is bad, just that particular simplex shape got stuck.
#'
#' @return An object of class `"dinarch_fit"`: a list with elements
#'   `method`, `coefficients` (`a`, `b`, `phi`, `beta`), `k`,
#'   `population_used`, `covariates_used`, `loglik`, `n_obs`,
#'   `convergence`, `optim_restarts`, and `call`.
#'
#' @export
dinarch_fit_ml <- function(data,
                            y = "y", index = "index", group = NULL,
                            population = NULL, covariates = NULL,
                            k = 1, lambda = NULL,
                            start = NULL,
                            optimizer = c("nelder-mead+bfgs", "nelder-mead"),
                            control = list(),
                            max_restarts = 3) {

  optimizer <- match.arg(optimizer)
  stopifnot(is.numeric(k), length(k) == 1, k >= 1, k == round(k))

  prep <- .prepare_dinarch_data(data, y, index, group, population, covariates, k)

  y_vec   <- prep$y_vec
  pop_vec <- prep$pop_vec
  lag_mat <- prep$lag_mat
  Xcov    <- prep$Xcov
  n_beta  <- prep$n_beta
  n       <- prep$n

  # --- parameter (un)packing ----------------------------------------------
  unpack <- function(par) {
    pos <- 1L
    a <- par[pos]

    b <- stats::plogis(par[pos + seq_len(k)])
    pos <- pos + k

    pos <- pos + 1L
    phi <- exp(par[pos])

    beta <- numeric(0)
    if (n_beta > 0) {
      beta <- par[pos + seq_len(n_beta)]
    }

    list(a = a, b = b, phi = phi, beta = beta)
  }

  n_par <- 1L + k + 1L + n_beta

  if (is.null(start)) {
    # Data-driven starting values, rather than a fixed a = 0 (baseline
    # multiplier = 1) that can be wildly off when `population` is on a
    # different scale than the count series.
    mean_y   <- mean(y_vec)
    mean_pop <- mean(pop_vec)
    a_start  <- log(max(mean_y / mean_pop, 1e-6))

    # Rough negative-binomial method-of-moments estimate of phi, as a
    # starting point (Var = mean + mean^2/phi); falls back to a moderate
    # default if the data doesn't show overdispersion.
    var_y <- stats::var(y_vec)
    phi_start <- if (is.finite(var_y) && var_y > mean_y) {
      mean_y^2 / (var_y - mean_y)
    } else {
      5
    }
    phi_start <- min(max(phi_start, 0.5), 50)

    start <- c(a_start, rep(-1, k), log(phi_start), rep(0, n_beta))
  } else if (length(start) != n_par) {
    stop(sprintf("`start` must have length %d, not %d.", n_par, length(start)))
  }

  # --- (negative, optionally penalized) log-likelihood ---------------------
  neg_loglik <- function(par) {
    pars <- unpack(par)

    baseline <- exp(pars$a) * pop_vec
    lag_term <- as.numeric(lag_mat %*% pars$b)
    cov_term <- if (n_beta > 0) as.numeric(Xcov %*% exp(pars$beta)) else 0
    mu <- baseline + lag_term + cov_term

    if (any(!is.finite(mu)) || any(mu <= 0) || sum(pars$b) >= 1 || pars$phi <= 0) {
      return(1e10)
    }

    r <- pars$phi * mu
    p <- pars$phi / (1 + pars$phi)
    ll <- sum(stats::dnbinom(y_vec, size = r, prob = p, log = TRUE))

    penalty <- 0
    if (!is.null(lambda)) {
      penalty <- (-lambda * log(1 - max(pars$b) + 1e-5) -
        10 * lambda * log(min(pars$phi, 1) + 1e-5)) / n
    }

    -ll + penalty
  }

  # --- optimize, with automatic restarts on non-convergence -----------------
  nm_control <- utils::modifyList(list(maxit = 2000), control)

  run_nm <- function(par) {
    stats::optim(par = par, fn = neg_loglik, method = "Nelder-Mead", control = nm_control)
  }

  fit <- run_nm(start)

  restarts_used <- 0L
  while (fit$convergence != 0 && restarts_used < max_restarts) {
    # Nelder-Mead's simplex can degenerate (convergence code 10); the
    # standard remedy is to restart from the last result with a fresh
    # simplex, rather than treating this as a hard failure.
    fit <- run_nm(fit$par)
    restarts_used <- restarts_used + 1L
  }

  if (optimizer == "nelder-mead+bfgs") {
    bfgs_control <- utils::modifyList(list(maxit = 2000), control)
    fit2 <- try(
      stats::optim(par = fit$par, fn = neg_loglik, method = "BFGS", control = bfgs_control),
      silent = TRUE
    )
    if (!inherits(fit2, "try-error") && fit2$value <= fit$value) {
      fit <- fit2
    }
  }

  if (fit$convergence != 0) {
    warning(
      "optim() did not report successful convergence after ", restarts_used,
      " restart(s) (code ", fit$convergence, "). Consider passing custom ",
      "`start` values, `control = list(maxit = ...)`, or check whether ",
      "the data genuinely supports this model (e.g. too few non-zero ",
      "observations)."
    )
  }

  final <- unpack(fit$par)

  structure(
    list(
      method = "ml",
      coefficients = list(a = final$a, b = final$b, phi = final$phi, beta = final$beta),
      k = k,
      population_used = !is.null(population),
      covariates_used = covariates,
      loglik = -fit$value,
      n_obs = n,
      convergence = fit$convergence,
      optim_restarts = restarts_used,
      call = match.call()
    ),
    class = "dinarch_fit"
  )
}
