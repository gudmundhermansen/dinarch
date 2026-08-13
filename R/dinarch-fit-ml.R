# ML estimation -----------------------------------------------------------

#' Fit a DINARCH model by maximum likelihood
#'
#' Fits \deqn{\eta_t = \sum_j \beta_j z_{t,j}}
#' \deqn{\mu_t = \exp(\eta_t) + \sum_{j=1}^k b_j y_{t-j}}
#' \deqn{y_t \sim \mathrm{NegBinom}(r_t = \phi \mu_t,\; p_t = \phi/(1+\phi))}
#' by maximizing the log-likelihood via [stats::optim()], where `k` above
#' is `n_lags` below, and \eqn{z_{t,\cdot}} is the row of
#' `stats::model.matrix(formula, data)` for observation `t` (see
#' `formula` below).
#'
#' @param data A data.frame or data.table containing at least the `y` and
#'   `index` columns (and `group` and any variables referenced by
#'   `formula` if used).
#' @param y,index,group Column names (strings) in `data`. `group`
#'   identifies separate series (lags are never taken across a group
#'   boundary); if `NULL`, all rows are treated as one series.
#' @param formula A one-sided formula specifying the linear predictor
#'   `eta` inside `exp(eta)` - e.g. `~ pop_log + gdp`, `~ pop_log * gdp`
#'   for main effects plus their interaction, `~ dem:pop_log` for an
#'   interaction only, or `~ 0 + pop_log` to drop the intercept. Built via
#'   [stats::model.matrix()], so ordinary R formula conventions apply
#'   (interactions via `:`, `*`; intercept dropped via `- 1` or `0 +`).
#'   Variables are resolved against `data`, so inline transforms like
#'   `~ log(population) + gdp` work directly. Default `~1` (intercept
#'   only, i.e. `mu_t = exp(beta_0) + lag term`, the model's simplest
#'   case). There is no separate `population` argument - if `mu` should
#'   scale proportionally with population the way `stats::glm()`'s
#'   `offset()` scales a Poisson/NB mean, include `log(population)` in
#'   `formula` (a freely estimated coefficient here nests that as the
#'   special case where it equals 1, i.e. `population^1`, while also
#'   allowing sub/super-linear scaling if the data supports it).
#' @param n_lags Number of autoregressive lags.
#' @param start Optional custom starting values for the optimizer, on the
#'   internal unconstrained scale: `n_lags` stick-breaking break-fraction
#'   logits for `b` (see Details), then `log(phi)`, then `ncol(Xcov)`
#'   values for `beta` (`Xcov` from `formula`), in that order. If `NULL`
#'   (default), data-driven starting values are used: the intercept
#'   coefficient (if `formula` includes one) is initialized from
#'   `log(mean(y))` (other `beta` entries start at 0), rather than a fixed
#'   value - matters when the response is far from 1 on its natural
#'   scale, a common cause of `optim()` starting far from any reasonable
#'   optimum.
#' @param optim_method Optimization method passed to [stats::optim()] - any
#'   of `"Nelder-Mead"` (default), `"BFGS"`, `"CG"`, `"L-BFGS-B"`, `"SANN"`,
#'   or `"Brent"`. Named `optim_method` rather than `method` so it can be
#'   forwarded through [dinarch_fit()]'s `...` without colliding with that
#'   function's own `method` argument (which picks `"ml"` vs `"bayes"`).
#'   `"Nelder-Mead"` is derivative-free and a robust default; gradient-based
#'   methods can be faster once near the optimum.
#' @param refine If `TRUE` (default) and `optim_method` is not already
#'   `"BFGS"`,
#'   follow the primary optimization with a `"BFGS"` polish step, keeping
#'   it only if it improves on the primary result. This is a standard way
#'   to sharpen a derivative-free search's result once it's in the basin of
#'   the optimum.
#' @param n_starts Number of starting points to optimize from (default 1).
#'   For `n_starts > 1`, the first start is `start` (or the data-driven
#'   default) and the remaining `n_starts - 1` are `start` perturbed by
#'   Gaussian noise (`start_jitter_sd`) on the unconstrained scale; the
#'   result with the best (lowest) negative log-likelihood across all
#'   starts is returned. This is a standard guard against local optima -
#'   more valuable the more `n_lags` and the number of covariates grow.
#' @param start_jitter_sd Standard deviation of the Gaussian noise used to
#'   perturb `start` for `n_starts > 1`. Ignored if `n_starts == 1`.
#' @param control Optional list passed to [stats::optim()]'s `control`
#'   argument (e.g. `list(maxit = 5000, reltol = 1e-10)`), merged over the
#'   package's defaults - use this to tune convergence behavior if the
#'   defaults struggle on a particular dataset.
#' @param max_restarts If `optim()` reports a non-zero convergence code
#'   (most commonly code 10 for Nelder-Mead, simplex degeneracy),
#'   automatically restart from the last result up to this many times -
#'   the standard remedy for this failure mode, since a degenerate simplex
#'   doesn't mean the optimum itself is bad, just that particular simplex
#'   shape got stuck.
#' @param seed Optional random seed, used only to make the `start_jitter_sd`
#'   perturbation in multi-start (`n_starts > 1`) reproducible.
#'
#' @details `b` is recovered from a stick-breaking transform: `n_lags`
#'   unconstrained values are mapped via `plogis()` to break fractions
#'   `z_1, ..., z_{n_lags}` in `(0, 1)`, and
#'   `b_j = z_j * prod_{i<j}(1 - z_i)` is the fraction of a unit-length
#'   "stick" broken off at step `j`. This guarantees
#'   `sum(b) = 1 - prod_i(1 - z_i) < 1` - the model's stationarity
#'   condition (see [dinarch_simulate()]) - for *every* real parameter
#'   vector, with no infeasible region for `optim()` to land in.
#'
#'   `mu = exp(eta) + lag term` is strictly positive for *any* real `eta`
#'   and any real covariate values in `formula` - `exp()` can't be
#'   negative or zero for finite input - so there is no positivity
#'   requirement on covariate values and no infeasible region there
#'   either. The log-likelihood only falls back to a `-1e10` hard wall for
#'   the residual, essentially unreachable case of floating-point
#'   underflow (`mu` or `phi` rounding to exactly 0).
#'
#' @return An object of class `"dinarch_fit"`: a list with elements
#'   `method`, `coefficients` (`b`, `phi`, `beta` - `beta` named from
#'   `colnames(model.matrix(formula, data))`), `n_lags`, `formula`, `data`
#'   (a `data.table` with `group`, `index`, `y`, and the raw covariate
#'   columns `formula` references, for *every* row of the input `data`
#'   (including the first `n_lags` per group, which aren't part of the
#'   likelihood but are needed to seed lags) - used by
#'   [dinarch_project()] to seed forward simulation or in-sample replay
#'   without needing the original `data` again), `loglik`, `n_obs`,
#'   `convergence`, `optim_method`, `optim_restarts`, `n_starts`, and
#'   `call`.
#'
#' @export
dinarch_fit_ml <- function(data,
                            y = "y", index = "index", group = NULL,
                            formula = ~1,
                            n_lags = 1,
                            start = NULL,
                            optim_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                            refine = TRUE,
                            n_starts = 1,
                            start_jitter_sd = 0.5,
                            control = list(),
                            max_restarts = 3,
                            seed = NULL) {

  optim_method <- match.arg(optim_method)
  formula <- stats::as.formula(formula)
  stopifnot(is.numeric(n_lags), length(n_lags) == 1, n_lags >= 1, n_lags == round(n_lags))
  stopifnot(is.numeric(n_starts), length(n_starts) == 1, n_starts >= 1, n_starts == round(n_starts))

  prep <- .prepare_dinarch_data(data, y, index, group, formula, n_lags)

  y_vec      <- prep$y_vec
  lag_mat    <- prep$lag_mat
  Xcov       <- prep$Xcov
  n_beta     <- prep$n_beta
  n          <- prep$n
  beta_names <- colnames(Xcov)

  # --- parameter (un)packing ----------------------------------------------
  # b is recovered from a stick-breaking transform: z = plogis(par) gives
  # n_lags break fractions in (0, 1), and b_j = z_j * prod_{i<j}(1 - z_i)
  # is the fraction of the still-unbroken "stick" (starting at length 1)
  # broken off at step j. This makes sum(b) = 1 - prod_i(1 - z_i) < 1 for
  # *any* real `par` - the stationarity condition (see this function's
  # Details / dinarch_simulate()'s stationarity check) holds by
  # construction, with no infeasible region for optim() to land in or get
  # stuck against.
  unpack <- function(par) {
    pos <- 0L
    z <- stats::plogis(par[pos + seq_len(n_lags)])
    remaining <- cumprod(c(1, 1 - z))[seq_len(n_lags)]
    b <- z * remaining
    pos <- pos + n_lags

    pos <- pos + 1L
    phi <- exp(par[pos])

    beta <- par[pos + seq_len(n_beta)]

    list(b = b, phi = phi, beta = beta)
  }

  n_par <- n_lags + 1L + n_beta

  if (is.null(start)) {
    # Data-driven starting value for the intercept (if formula includes
    # one), rather than a fixed 0, which can be wildly off when the
    # response isn't O(1) on its natural scale.
    beta_start <- rep(0, n_beta)
    intercept_pos <- which(beta_names == "(Intercept)")
    if (length(intercept_pos) == 1) {
      beta_start[intercept_pos] <- log(max(mean(y_vec), 1e-6))
    }

    # Rough negative-binomial method-of-moments estimate of phi, as a
    # starting point (Var = mean + mean^2/phi); falls back to a moderate
    # default if the data doesn't show overdispersion.
    mean_y <- mean(y_vec)
    var_y  <- stats::var(y_vec)
    phi_start <- if (is.finite(var_y) && var_y > mean_y) {
      mean_y^2 / (var_y - mean_y)
    } else {
      5
    }
    phi_start <- min(max(phi_start, 0.5), 50)

    # Constant break fraction z0 chosen so the starting sum(b) is always
    # ~0.27 (= plogis(-1)) regardless of n_lags, rather than growing with
    # n_lags: 1 - (1 - z0)^n_lags = 0.27.
    b_sum_start <- stats::plogis(-1)
    z0 <- 1 - (1 - b_sum_start)^(1 / n_lags)
    b_start <- stats::qlogis(z0)

    start <- c(rep(b_start, n_lags), log(phi_start), beta_start)
  } else if (length(start) != n_par) {
    stop(sprintf("`start` must have length %d, not %d.", n_par, length(start)))
  }

  # --- negative log-likelihood ---------------------------------------------
  neg_loglik <- function(par) {
    pars <- unpack(par)

    eta <- as.numeric(Xcov %*% pars$beta)
    lag_term <- as.numeric(lag_mat %*% pars$b)
    mu <- exp(eta) + lag_term

    # mu > 0 holds by construction for any finite `par` (see Details);
    # this only guards the residual case of floating-point underflow.
    if (any(!is.finite(mu)) || any(mu <= 0) || pars$phi <= 0) {
      return(1e10)
    }

    r <- pars$phi * mu
    p <- pars$phi / (1 + pars$phi)
    -sum(stats::dnbinom(y_vec, size = r, prob = p, log = TRUE))
  }

  # --- one optimization run: primary method, restarts on non-convergence,
  # optional BFGS polish -----------------------------------------------------
  opt_control <- utils::modifyList(list(maxit = 2000), control)

  run_once <- function(par0) {
    fit <- stats::optim(par = par0, fn = neg_loglik, method = optim_method, control = opt_control)

    restarts_used <- 0L
    while (fit$convergence != 0 && restarts_used < max_restarts) {
      # optim()'s reported non-convergence (e.g. Nelder-Mead's simplex
      # degeneracy, code 10) doesn't mean the optimum itself is bad, just
      # that the search got stuck - the standard remedy is to restart from
      # the last result with a fresh search.
      fit <- stats::optim(par = fit$par, fn = neg_loglik, method = optim_method, control = opt_control)
      restarts_used <- restarts_used + 1L
    }

    if (refine && optim_method != "BFGS") {
      fit2 <- try(
        stats::optim(par = fit$par, fn = neg_loglik, method = "BFGS", control = opt_control),
        silent = TRUE
      )
      if (!inherits(fit2, "try-error") && is.finite(fit2$value) && fit2$value <= fit$value) {
        fit <- fit2
      }
    }

    list(fit = fit, restarts_used = restarts_used)
  }

  # --- multi-start: guard against local optima by trying several starting
  # points and keeping the best result --------------------------------------
  if (!is.null(seed)) set.seed(seed)

  starts <- vector("list", n_starts)
  starts[[1]] <- start
  if (n_starts > 1) {
    for (s in 2:n_starts) {
      starts[[s]] <- start + stats::rnorm(n_par, sd = start_jitter_sd)
    }
  }

  runs <- lapply(starts, run_once)
  values <- vapply(runs, function(r) r$fit$value, numeric(1))
  best <- runs[[which.min(values)]]
  fit <- best$fit
  restarts_used <- best$restarts_used

  if (fit$convergence != 0) {
    warning(
      "optim() did not report successful convergence after ", restarts_used,
      " restart(s) across ", n_starts, " starting point(s) (code ", fit$convergence,
      "). Consider passing custom `start` values, a larger `n_starts`, ",
      "`control = list(maxit = ...)`, or check whether the data genuinely ",
      "supports this model (e.g. too few non-zero observations)."
    )
  }

  final <- unpack(fit$par)
  names(final$beta) <- beta_names

  structure(
    list(
      method = "ml",
      coefficients = list(b = final$b, phi = final$phi, beta = final$beta),
      n_lags = n_lags,
      formula = formula,
      data = prep$full_data,
      loglik = -fit$value,
      n_obs = n,
      convergence = fit$convergence,
      optim_method = optim_method,
      optim_restarts = restarts_used,
      n_starts = n_starts,
      call = match.call()
    ),
    class = "dinarch_fit"
  )
}
