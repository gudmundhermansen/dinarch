# Simulation --------------------------------------------------------------

#' Simulate from a DINARCH process
#'
#' Simulates a Negative-Binomial Dynamic INnovation ARCH (DINARCH) count
#' time series:
#' \deqn{\mu_t = \exp(a_t) \cdot population_t + \sum_{j=1}^k b_{t,j} y_{t-j}
#'   + \sum_l \exp(\beta_l) covariate_{t,l}}
#' \deqn{y_t \sim \mathrm{NegBinom}(r_t = \phi_t \mu_t,\; p_t = \phi_t/(1+\phi_t))}
#'
#' @param n Number of periods to simulate.
#' @param a Baseline log-rate. Scalar or length-n vector.
#' @param b Autoregressive (lag) coefficients, each in `[0, 1)`. A scalar
#'   (k = 1, constant over time), a length-k vector (k lags, constant over
#'   time), or an n x k matrix (k lags, fully time-varying).
#' @param phi Dispersion/concentration parameter (> 0). Scalar or
#'   length-n vector. Larger `phi` means less overdispersion.
#' @param population Optional population path: `NULL` (default), a scalar,
#'   or a length-n vector. When `NULL`, the baseline term reduces to
#'   `exp(a)`.
#' @param dynamic_population If `TRUE`, cumulative simulated deaths (within
#'   the `n` simulated periods) are subtracted from `population` as the
#'   simulation progresses, floored at 0. Requires `population`.
#' @param covariates Optional n x p matrix (or data.frame coercible to one)
#'   of additional covariates, each entering as
#'   `exp(beta_l) * covariates[, l]`. Off by default.
#' @param beta Coefficients for `covariates`, length p. Required if
#'   `covariates` is supplied.
#' @param y_init Optional vector of the last k observed y values, oldest to
#'   newest, seeding the autoregressive lags. If `NULL`, the process is
#'   warmed up stochastically for `burn_in` extra steps under the period-1
#'   parameter values (population/covariates also held at their period-1
#'   values), and the last k draws are used as the seed.
#' @param burn_in Number of stochastic warm-up steps used when `y_init` is
#'   `NULL`. Ignored if `y_init` is supplied.
#' @param threshold Optional list with elements `window`, `fraction`, and
#'   `pop_unit`. Once at least `window` periods of history exist: if the
#'   sum of the last `window` simulated deaths exceeds
#'   `fraction * pop_unit * population_t` (the *current*, possibly
#'   depleted, population), the current period's simulated value is set to
#'   0. `NULL` (default) disables this. Requires `population`.
#' @param seed Optional random seed.
#'
#' @return A `data.table` with columns `index` (`1:n`), `y`, and
#'   `population` (the population path actually used, after any dynamic
#'   depletion; equal to the input `population` if `dynamic_population` is
#'   `FALSE`). If population depletes to 0 with no recent autoregressive
#'   momentum, the mean drops to exactly 0 and subsequent periods
#'   deterministically simulate `y = 0` (an absorbing "process has ended"
#'   state) rather than erroring.
#'
#' @examples
#' # simple constant-parameter AR(1) count series, no population
#' dinarch_simulate(n = 10, a = log(5), b = 0.3, phi = 10, seed = 1)
#'
#' # with a population path and dynamic depletion
#' dinarch_simulate(
#'   n = 10, a = log(0.01), b = 0.3, phi = 5,
#'   population = rep(1000, 10), dynamic_population = TRUE,
#'   threshold = list(window = 5, fraction = 0.03, pop_unit = 1e6),
#'   seed = 1
#' )
#'
#' @export
dinarch_simulate <- function(n,
                              a, b, phi,
                              population = NULL,
                              dynamic_population = FALSE,
                              covariates = NULL,
                              beta = NULL,
                              y_init = NULL,
                              burn_in = 0,
                              threshold = NULL,
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.numeric(n), length(n) == 1, n >= 1, n == round(n))
  stopifnot(
    is.numeric(burn_in), length(burn_in) == 1, burn_in >= 0,
    burn_in == round(burn_in)
  )

  # --- baseline parameters ---------------------------------------------
  a   <- .recycle(a, n)
  phi <- .recycle(phi, n)
  .check_positive(phi)

  B <- .as_lag_matrix(b, n)
  k <- ncol(B)
  if (any(B < 0 | B >= 1, na.rm = TRUE)) {
    stop("`b` must lie in [0, 1).")
  }
  if (any(rowSums(B) >= 1)) {
    warning(
      "Some rows of `b` sum to >= 1; the simulated process may be ",
      "explosive/non-stationary."
    )
  }

  # --- population --------------------------------------------------------
  has_population <- !is.null(population)
  if (has_population) {
    pop <- .recycle(population, n)
    .check_positive(pop)
  } else {
    if (dynamic_population) {
      stop("`dynamic_population = TRUE` requires `population` to be supplied.")
    }
    if (!is.null(threshold)) {
      stop("`threshold` requires `population` to be supplied.")
    }
    pop <- rep(1, n)
  }

  # --- additional covariates ---------------------------------------------
  if (!is.null(covariates)) {
    Xcov <- as.matrix(covariates)
    if (nrow(Xcov) != n) {
      stop(sprintf("`covariates` must have %d rows (= n), not %d.", n, nrow(Xcov)))
    }
    if (is.null(beta) || length(beta) != ncol(Xcov)) {
      stop("`beta` must be supplied with length equal to ncol(covariates).")
    }
    cov_term <- as.numeric(Xcov %*% exp(beta))
  } else {
    cov_term <- rep(0, n)
  }

  # --- threshold rule ------------------------------------------------
  if (!is.null(threshold)) {
    if (!all(c("window", "fraction", "pop_unit") %in% names(threshold))) {
      stop("`threshold` must be a list with elements `window`, `fraction`, `pop_unit`.")
    }
  }

  # --- seed the k autoregressive lags -------------------------------------
  if (!is.null(y_init)) {
    if (length(y_init) != k) {
      stop(sprintf("`y_init` must have length %d (= k).", k))
    }
    y_history <- as.numeric(y_init)
  } else {
    y_history <- .warm_up_dinarch(
      k = k, burn_in = burn_in,
      a0 = a[1], b0 = B[1, ], phi0 = phi[1],
      pop0 = pop[1], cov0 = cov_term[1]
    )
  }

  # --- main simulation loop -----------------------------------------------
  y <- numeric(n)
  pop_used <- numeric(n)
  cum_deaths <- 0

  for (t in seq_len(n)) {

    pop_t <- if (dynamic_population) max(pop[t] - cum_deaths, 0) else pop[t]
    pop_used[t] <- pop_t

    lag_term <- sum(B[t, ] * rev(utils::tail(y_history, k)))
    mu_t <- exp(a[t]) * pop_t + lag_term + cov_term[t]

    if (!is.finite(mu_t) || mu_t < 0) {
      stop(sprintf(
        "Negative or non-finite mean at period %d (mu = %s). Check `a`, `population`, `b`, or `covariates`.",
        t, format(mu_t)
      ))
    }

    if (mu_t == 0) {
      # e.g. population fully depleted with no recent autoregressive
      # momentum: the only sensible outcome is no event this period.
      y[t] <- 0
    } else {
      r_t <- phi[t] * mu_t
      p_t <- phi[t] / (1 + phi[t])
      y[t] <- stats::rnbinom(1, size = r_t, prob = p_t)
    }

    if (!is.null(threshold) && length(y_history) >= threshold$window) {
      recent <- sum(utils::tail(y_history, threshold$window))
      if (recent > threshold$fraction * threshold$pop_unit * pop_t) {
        y[t] <- 0
      }
    }

    if (dynamic_population) cum_deaths <- cum_deaths + y[t]
    y_history <- c(y_history, y[t])
  }

  data.table::data.table(index = seq_len(n), y = y, population = pop_used)
}

#' Stochastic warm-up to seed the initial k lags under constant (period-1)
#' parameter values. Starts from an all-zero pre-history and returns the
#' last k generated values.
#' @noRd
.warm_up_dinarch <- function(k, burn_in, a0, b0, phi0, pop0, cov0 = 0) {
  total <- k + burn_in
  history <- rep(0, k)

  for (t in seq_len(total)) {
    lag_vals <- rev(utils::tail(history, k))
    mu <- exp(a0) * pop0 + sum(b0 * lag_vals) + cov0

    if (!is.finite(mu) || mu < 0) {
      stop("Negative or non-finite mean during warm-up; check `a`, `population`, or `b`.")
    }

    if (mu == 0) {
      y_t <- 0
    } else {
      r <- phi0 * mu
      p <- phi0 / (1 + phi0)
      y_t <- stats::rnbinom(1, size = r, prob = p)
    }
    history <- c(history, y_t)
  }

  utils::tail(history, k)
}
