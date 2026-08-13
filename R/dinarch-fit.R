# Unified fit dispatcher ----------------------------------------------------

#' Fit a DINARCH model
#'
#' A thin dispatcher over [dinarch_fit_ml()] and [dinarch_fit_bayes()],
#' giving a single entry point that switches on `method`. Method-specific
#' arguments (e.g. `optim_method`/`start`/`n_starts`/`control`/
#' `max_restarts` for `"ml"`; `prior`/`iter`/`chains`/`seed` for
#' `"bayes"`) are not duplicated here - they are documented on, and
#' forwarded via `...` straight through to, the underlying function, so
#' the two methods' argument sets can't drift out of sync with their real
#' implementations.
#'
#' @param data A data.frame or data.table; see [dinarch_fit_ml()].
#' @param ... Further arguments passed to [dinarch_fit_ml()] (if
#'   `method = "ml"`) or [dinarch_fit_bayes()] (if `method = "bayes"`),
#'   e.g. `y`, `index`, `group`, `formula`, `n_lags`, plus any
#'   method-specific arguments.
#' @param method `"ml"` (default) for maximum-likelihood estimation via
#'   [dinarch_fit_ml()], or `"bayes"` for Bayesian estimation via
#'   [dinarch_fit_bayes()].
#'
#' @return An object of class `"dinarch_fit"`; see [dinarch_fit_ml()] /
#'   [dinarch_fit_bayes()] for the full return value.
#'
#' @examples
#' sim <- dinarch_simulate(n = 200, b = 0.3, phi = 8, beta = log(4), seed = 1)
#' dat <- data.frame(y = sim$y, index = sim$index)
#'
#' fit <- dinarch_fit(dat, n_lags = 1, method = "ml")
#' fit
#'
#' @export
dinarch_fit <- function(data, ..., method = c("ml", "bayes")) {
  method <- match.arg(method)
  switch(
    method,
    ml    = dinarch_fit_ml(data, ...),
    bayes = dinarch_fit_bayes(data, ...)
  )
}
