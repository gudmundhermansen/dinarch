# Simulation --------------------------------------------------------------

#' Simulate from a DINARCH process
#'
#' Simulates a Negative-Binomial Dynamic INnovation ARCH (DINARCH) count
#' time series:
#' \deqn{\eta_t = \sum_j \beta_{t,j} z_{t,j}}
#' \deqn{\mu_t = \exp(\eta_t) + \sum_{j=1}^k b_{t,j} y_{t-j}}
#' \deqn{y_t \sim \mathrm{NegBinom}(r_t = \phi_t \mu_t,\; p_t = \phi_t/(1+\phi_t))}
#' `k` above is the number of lags, inferred from `b` (see below) rather
#' than passed separately, and \eqn{z_{t,\cdot}} is the row of
#' `stats::model.matrix(formula, covariates)` for period `t` (see
#' `formula` below). Every parameter (`b`, `phi`, `beta`) follows the same
#' convention: a scalar or matching-length vector (constant over time), or
#' a full `n`-row path (`n x k` matrix for `b`/`beta`, length-`n` vector
#' for `phi`) for fully time-varying parameters - see [dinarch_project()]'s
#' `new_para` for how a fitted model's forward `b`/`phi`/`beta` path can
#' be overridden the same way. This is exactly the mean structure
#' [dinarch_fit_ml()]/[dinarch_fit_bayes()] fit - there is no distinguished
#' `population` term here either; if `mu` should depend on population,
#' include it as a `formula`/`covariates` term (e.g. `~ log(population)`)
#' like any other covariate. See `y_threshold` below for the *only* place
#' population-like quantities still get special treatment in this
#' function, and note it does not feed back into `mu` at all.
#'
#' @param n Number of periods to simulate.
#' @param b Autoregressive (lag) coefficients, each in `[0, 1)`. A scalar
#'   (`n_lags = 1`, constant over time), a length-`n_lags` vector
#'   (`n_lags` lags, constant over time), or an `n x n_lags` matrix
#'   (`n_lags` lags, fully time-varying).
#' @param phi Dispersion/concentration parameter (> 0). Scalar or
#'   length-n vector. Larger `phi` means less overdispersion.
#' @param formula A one-sided formula specifying the linear predictor
#'   `eta` inside `exp(eta)` - e.g. `~ gdp`, `~ gdp * dem` for main
#'   effects plus their interaction, `~ log(population)` to give the mean
#'   a proportional-to-population baseline, or `~ 0 + gdp` to drop the
#'   intercept. Built via [stats::model.matrix()] against `covariates`,
#'   so ordinary R formula conventions apply. Default `~1` (intercept
#'   only, i.e. `mu_t = exp(beta) + lag term` when `covariates` is
#'   `NULL`) - matches [dinarch_fit_ml()].
#' @param covariates Optional data.frame (or matrix) with `n` rows holding
#'   the *raw* covariate paths referenced by `formula` (e.g. a `gdp`, or
#'   `population`, column) - not a pre-built design matrix. `NULL`
#'   (default) simulates with `formula`'s intercept only.
#' @param beta Coefficients for `model.matrix(formula, covariates)`, each
#'   in the linear predictor `eta` (see Details) - length
#'   `ncol(model.matrix(formula, covariates))` (includes the intercept
#'   coefficient if `formula` has one - e.g. for the default
#'   `formula = ~1`, `covariates = NULL`, `beta` is a single number: the
#'   baseline log-rate), or an `n x ncol(model.matrix(...))` matrix for
#'   fully time-varying
#'   coefficients (e.g. a covariate *effect* that ramps up over a
#'   scenario horizon, as opposed to a time-varying covariate *value*,
#'   which instead varies `covariates` itself).
#' @param y_init Optional vector of the last `n_lags` observed y values,
#'   oldest to newest, seeding the autoregressive lags. If `NULL`, the
#'   process is warmed up stochastically for `burn_in` extra steps under
#'   the period-1 parameter values (covariates also held at their
#'   period-1 values), and the last `n_lags` draws are used as the seed.
#' @param burn_in Number of stochastic warm-up steps used when `y_init` is
#'   `NULL`. Ignored if `y_init` is supplied.
#' @param y_threshold Optional list with elements `window` (scalar) and
#'   `limit` (scalar, or length-n vector for a time-varying cap, matching
#'   the `b`/`phi`/`beta` convention). Once at least `window` periods of
#'   history exist: if the sum of the last `window` simulated values
#'   exceeds `limit_t`, the current period's simulated value is set to 0.
#'   `NULL` (default) disables this.
#'
#'   `limit` is a plain number, not a population-derived quantity - this
#'   is deliberately independent of `mu`/`formula`/`covariates`, so
#'   changing one can never silently change the other. If the limit you
#'   actually want is "X% of population," compute it yourself first (e.g.
#'   `limit = fraction * population`, either a single number or a
#'   length-n vector) and pass the result directly.
#'
#'   This is a *reactive*, rolling-window rule, not a hard per-period cap:
#'   it compares history *before* the current draw against the limit, so
#'   it stops sustained escalation once a bad stretch has already
#'   happened, but doesn't constrain any single period's draw directly -
#'   a single very large outlier can still occur before the rule reacts
#'   (though a conservatively chosen `limit` means a single serious
#'   outlier is often enough to trigger it on the *next* period). With
#'   `n_lags > 1` a zeroed period still enters the lag history like any
#'   other value, which can interact with the rule in less obvious ways.
#' @param seed Optional random seed.
#'
#' @return A `data.table` with columns `index` (`1:n`) and `y`.
#'
#' @examples
#' # simple constant-parameter AR(1) count series
#' dinarch_simulate(n = 10, b = 0.3, phi = 10, beta = log(5), seed = 1)
#'
#' # population as an ordinary covariate (proportional effect on the mean,
#' # via log(population) with beta = 1), plus a cap forcing y to 0 once
#' # the last 5 periods sum past 3% of population
#' pop <- rep(1000, 10)
#' dinarch_simulate(
#'   n = 10, b = 0.3, phi = 5,
#'   formula = ~ log(population), covariates = data.frame(population = pop),
#'   beta = c(log(0.01), 1),
#'   y_threshold = list(window = 5, limit = 0.03 * pop),
#'   seed = 1
#' )
#'
#' @export
dinarch_simulate <- function(n,
                              b, phi,
                              formula = ~1,
                              covariates = NULL,
                              beta,
                              y_init = NULL,
                              burn_in = 0,
                              y_threshold = NULL,
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.numeric(n), length(n) == 1, n >= 1, n == round(n))
  stopifnot(
    is.numeric(burn_in), length(burn_in) == 1, burn_in >= 0,
    burn_in == round(burn_in)
  )

  # baseenv() so a variable referenced by `formula` but missing from
  # `covariates` errors clearly, instead of model.matrix() silently
  # falling back to a same-named variable wherever `formula` was written
  # (base functions like log() remain usable).
  formula <- stats::as.formula(formula)
  environment(formula) <- baseenv()
  phi <- .recycle(phi, n)
  .check_positive(phi)

  B <- .as_coef_matrix(b, n)
  n_lags <- ncol(B)
  if (any(B < 0 | B >= 1, na.rm = TRUE)) {
    stop("`b` must lie in [0, 1).")
  }
  if (any(rowSums(B) >= 1)) {
    warning(
      "Some rows of `b` sum to >= 1; the simulated process may be ",
      "explosive/non-stationary."
    )
  }

  # --- covariate design matrix (formula) ----------------------------------
  cov_df <- if (is.null(covariates)) data.frame(row.names = seq_len(n)) else as.data.frame(covariates)
  if (nrow(cov_df) != n) {
    stop(sprintf("`covariates` must have %d rows (= n), not %d.", n, nrow(cov_df)))
  }
  Xcov <- stats::model.matrix(formula, data = cov_df)
  if (nrow(Xcov) != n) {
    stop("`formula` produced fewer rows than `n` - check `covariates` for missing values.")
  }
  n_beta <- ncol(Xcov)
  if (missing(beta)) {
    stop(sprintf(
      "`beta` must be supplied: length %d (= ncol(model.matrix(formula, covariates))), ",
      n_beta
    ), sprintf("or an n x %d matrix for time-varying coefficients.", n_beta))
  }
  Beta <- .as_coef_matrix(beta, n, n_beta)
  eta <- rowSums(Xcov * Beta)

  # --- y_threshold rule ------------------------------------------------
  # `limit` may be a scalar (constant) or a length-n vector (time-varying
  # cap), matching the b/phi/beta convention. Deliberately independent of
  # mu/formula/covariates - if `limit` should be "X% of population",
  # compute that yourself and pass the resulting number(s) directly.
  if (!is.null(y_threshold)) {
    if (!all(c("window", "limit") %in% names(y_threshold))) {
      stop("`y_threshold` must be a list with elements `window` and `limit`.")
    }
    y_threshold$limit <- .recycle(y_threshold$limit, n, "y_threshold$limit")
  }

  # --- seed the n_lags autoregressive lags --------------------------------
  if (!is.null(y_init)) {
    if (length(y_init) != n_lags) {
      stop(sprintf("`y_init` must have length %d (= n_lags).", n_lags))
    }
    y_history <- as.numeric(y_init)
  } else {
    y_history <- .warm_up_dinarch(
      n_lags = n_lags, burn_in = burn_in,
      eta0 = eta[1], b0 = B[1, ], phi0 = phi[1]
    )
  }

  # --- main simulation loop -----------------------------------------------
  y <- numeric(n)

  for (t in seq_len(n)) {

    lag_term <- sum(B[t, ] * rev(utils::tail(y_history, n_lags)))
    mu_t <- exp(eta[t]) + lag_term

    if (!is.finite(mu_t) || mu_t <= 0) {
      stop(sprintf(
        "Non-positive or non-finite mean at period %d (mu = %s). Check `beta`, `covariates`, or `b`.",
        t, format(mu_t)
      ))
    }

    r_t <- phi[t] * mu_t
    p_t <- phi[t] / (1 + phi[t])
    y[t] <- stats::rnbinom(1, size = r_t, prob = p_t)

    if (!is.null(y_threshold) && length(y_history) >= y_threshold$window) {
      recent <- sum(utils::tail(y_history, y_threshold$window))
      if (recent > y_threshold$limit[t]) {
        y[t] <- 0
      }
    }

    y_history <- c(y_history, y[t])
  }

  data.table::data.table(index = seq_len(n), y = y)
}

#' Stochastic warm-up to seed the initial n_lags lags under constant
#' (period-1) parameter values. Starts from an all-zero pre-history and
#' returns the last n_lags generated values.
#' @noRd
.warm_up_dinarch <- function(n_lags, burn_in, eta0, b0, phi0) {
  total <- n_lags + burn_in
  history <- rep(0, n_lags)

  for (t in seq_len(total)) {
    lag_vals <- rev(utils::tail(history, n_lags))
    mu <- exp(eta0) + sum(b0 * lag_vals)

    if (!is.finite(mu) || mu <= 0) {
      stop("Non-positive or non-finite mean during warm-up; check `beta`, `covariates`, or `b`.")
    }

    r <- phi0 * mu
    p <- phi0 / (1 + phi0)
    history <- c(history, stats::rnbinom(1, size = r, prob = p))
  }

  utils::tail(history, n_lags)
}
