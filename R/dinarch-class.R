# dinarch_fit class methods ------------------------------------------------

#' @export
print.dinarch_fit <- function(x, ...) {
  cat(sprintf("<dinarch_fit>  method = %s   k = %d\n", x$method, x$k))

  cat(sprintf("  a   = %s\n", format(x$coefficients$a, digits = 4)))
  cat(sprintf("  b   = %s\n", paste(format(x$coefficients$b, digits = 4), collapse = ", ")))
  cat(sprintf("  phi = %s\n", format(x$coefficients$phi, digits = 4)))

  if (length(x$coefficients$beta) > 0) {
    cat(sprintf("  beta = %s\n", paste(format(x$coefficients$beta, digits = 4), collapse = ", ")))
  }

  if (!is.null(x$loglik)) {
    cat(sprintf("Log-likelihood: %.2f  (n = %d)\n", x$loglik, x$n_obs))
  }

  invisible(x)
}

#' @export
coef.dinarch_fit <- function(object, ...) {
  object$coefficients
}

#' Summarize a DINARCH fit
#'
#' @param object A `"dinarch_fit"` object.
#' @param level Credible interval level for Bayesian fits (ignored for ML
#'   fits, which do not currently report uncertainty - see Details).
#' @param ... Unused.
#'
#' @details ML fits do not yet report standard errors (no Hessian is
#'   currently computed by [dinarch_fit_ml()]), only point estimates plus
#'   log-likelihood/AIC/BIC. Bayesian fits report posterior means with
#'   credible intervals computed from the stored posterior draws.
#'
#' @return An object of class `"summary.dinarch_fit"`.
#' @export
summary.dinarch_fit <- function(object, level = 0.95, ...) {

  coef_names <- c(
    "a",
    paste0("b[", seq_len(object$k), "]"),
    "phi",
    if (length(object$coefficients$beta) > 0) {
      paste0("beta[", seq_along(object$coefficients$beta), "]")
    }
  )

  estimate <- c(
    object$coefficients$a,
    object$coefficients$b,
    object$coefficients$phi,
    object$coefficients$beta
  )

  if (identical(object$method, "bayes") && !is.null(object$posterior)) {
    alpha <- 1 - level
    probs <- c(alpha / 2, 1 - alpha / 2)

    post_mat <- cbind(
      matrix(object$posterior$a, ncol = 1),
      object$posterior$b,
      matrix(object$posterior$phi, ncol = 1),
      if (length(object$coefficients$beta) > 0) object$posterior$beta
    )
    ci <- t(apply(post_mat, 2, stats::quantile, probs = probs))

    coef_table <- data.frame(
      term = coef_names, estimate = estimate,
      lower = ci[, 1], upper = ci[, 2],
      row.names = NULL
    )
  } else {
    coef_table <- data.frame(term = coef_names, estimate = estimate, row.names = NULL)
  }

  n_par <- length(estimate)
  aic <- if (!is.na(object$loglik)) -2 * object$loglik + 2 * n_par else NA_real_
  bic <- if (!is.na(object$loglik)) -2 * object$loglik + log(object$n_obs) * n_par else NA_real_

  structure(
    list(
      method = object$method, k = object$k,
      n_obs = object$n_obs, coefficients = coef_table,
      loglik = object$loglik, aic = aic, bic = bic,
      convergence = object$convergence, level = level
    ),
    class = "summary.dinarch_fit"
  )
}

#' @export
print.summary.dinarch_fit <- function(x, ...) {
  cat(sprintf("DINARCH fit (%s)   k = %d   n = %d\n", x$method, x$k, x$n_obs))
  cat(strrep("-", 50), "\n", sep = "")
  print(x$coefficients, row.names = FALSE, digits = 4)
  cat(strrep("-", 50), "\n", sep = "")

  if (identical(x$method, "ml")) {
    cat(sprintf("Log-likelihood: %.2f   AIC: %.2f   BIC: %.2f\n", x$loglik, x$aic, x$bic))
    if (!is.na(x$convergence) && x$convergence != 0) {
      cat("Note: optimizer did not report successful convergence.\n")
    }
  } else {
    cat(sprintf("Posterior %.0f%% credible intervals shown.\n", 100 * x$level))
    if (!is.na(x$convergence)) {
      note <- if (x$convergence > 1.1) "  (consider more iterations/chains)" else ""
      cat(sprintf("Max Rhat: %.3f%s\n", x$convergence, note))
    }
  }

  invisible(x)
}

#' Simulate forward from a fitted DINARCH model
#'
#' A convenience wrapper around [dinarch_simulate()] that pulls `a`, `b`,
#' `phi`, and `beta` directly from a `"dinarch_fit"` object, so you don't
#' have to copy coefficients out by hand.
#'
#' The parameter is named `nsim` (not `n`) specifically to match
#' [stats::simulate()]'s generic signature - calling `simulate(fit, n = 20)`
#' instead would have R's partial argument matching silently (and
#' incorrectly) bind `20` to `nsim` on the *generic* before dispatch even
#' happens.
#'
#' @param object A `"dinarch_fit"` object.
#' @param nsim Number of periods to simulate forward.
#' @param seed Optional random seed.
#' @param population,covariates,dynamic_population,threshold,y_init,burn_in
#'   Forward paths/settings for the simulated horizon - see
#'   [dinarch_simulate()]. These describe the *future*, so they are not
#'   taken from the fit automatically.
#' @param use_draw For Bayesian fits only: an integer index into the stored
#'   posterior draws. If supplied, that specific draw's parameters are used
#'   instead of the posterior mean - call this repeatedly with different
#'   indices to propagate full posterior uncertainty into a simulation
#'   ensemble. `NULL` (default) uses the point estimate (posterior mean for
#'   Bayesian fits, MLE for ML fits).
#' @param ... Unused.
#'
#' @return A `data.table`, as for [dinarch_simulate()].
#' @export
simulate.dinarch_fit <- function(object, nsim = 10, seed = NULL,
                                  population = NULL, covariates = NULL,
                                  dynamic_population = FALSE,
                                  threshold = NULL,
                                  y_init = NULL, burn_in = 0,
                                  use_draw = NULL, ...) {

  if (!is.null(use_draw)) {
    if (is.null(object$posterior)) {
      stop("`use_draw` is only available for Bayesian fits (this fit has no stored posterior draws).")
    }
    n_draws <- length(object$posterior$a)
    if (use_draw < 1 || use_draw > n_draws || use_draw != round(use_draw)) {
      stop(sprintf("`use_draw` must be an integer between 1 and %d.", n_draws))
    }
    a   <- object$posterior$a[use_draw]
    b   <- object$posterior$b[use_draw, ]
    phi <- object$posterior$phi[use_draw]
    beta <- if (length(object$coefficients$beta) > 0) object$posterior$beta[use_draw, ] else numeric(0)
  } else {
    a    <- object$coefficients$a
    b    <- object$coefficients$b
    phi  <- object$coefficients$phi
    beta <- object$coefficients$beta
  }

  if (isTRUE(object$population_used) && is.null(population)) {
    warning(
      "This fit used a `population` covariate; no forward `population` ",
      "path was supplied, so the baseline term will use population = 1, ",
      "which is likely not what you want."
    )
  }
  if (length(object$covariates_used) > 0 && is.null(covariates)) {
    warning(
      "This fit used additional covariates (",
      paste(object$covariates_used, collapse = ", "),
      "); no forward `covariates` were supplied, so that term will be ",
      "dropped from the simulation."
    )
  }

  dinarch_simulate(
    n = nsim,
    a = a, b = b, phi = phi,
    population = population,
    dynamic_population = dynamic_population,
    covariates = covariates,
    beta = if (length(beta) > 0) beta else NULL,
    y_init = y_init,
    burn_in = burn_in,
    threshold = threshold,
    seed = seed
  )
}
