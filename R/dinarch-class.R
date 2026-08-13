# dinarch_fit class methods ------------------------------------------------

#' @export
print.dinarch_fit <- function(x, ...) {
  cat(sprintf("<dinarch_fit>  method = %s   n_lags = %d\n", x$method, x$n_lags))

  cat(sprintf("  b   = %s\n", paste(format(x$coefficients$b, digits = 4), collapse = ", ")))
  cat(sprintf("  phi = %s\n", format(x$coefficients$phi, digits = 4)))
  cat(sprintf("  beta = %s\n", paste(
    sprintf("%s = %s", names(x$coefficients$beta), format(x$coefficients$beta, digits = 4)),
    collapse = ", "
  )))

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
    paste0("b[", seq_len(object$n_lags), "]"),
    "phi",
    names(object$coefficients$beta)
  )

  estimate <- c(
    object$coefficients$b,
    object$coefficients$phi,
    object$coefficients$beta
  )

  if (identical(object$method, "bayes") && !is.null(object$posterior)) {
    alpha <- 1 - level
    probs <- c(alpha / 2, 1 - alpha / 2)

    post_mat <- cbind(
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
      method = object$method, n_lags = object$n_lags,
      n_obs = object$n_obs, coefficients = coef_table,
      loglik = object$loglik, aic = aic, bic = bic,
      convergence = object$convergence, level = level
    ),
    class = "summary.dinarch_fit"
  )
}

#' @export
print.summary.dinarch_fit <- function(x, ...) {
  cat(sprintf("DINARCH fit (%s)   n_lags = %d   n = %d\n", x$method, x$n_lags, x$n_obs))
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
#' A convenience wrapper around [dinarch_simulate()] that pulls `b`,
#' `phi`, and `beta` (plus the fit's `formula`) directly from a
#' `"dinarch_fit"` object, so you don't have to copy coefficients out by
#' hand.
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
#' @param covariates,y_threshold,y_init,burn_in Forward paths/settings for
#'   the simulated horizon - see [dinarch_simulate()]. These describe the
#'   *future*, so they are not taken from the fit automatically;
#'   `covariates` must supply whatever variables the fit's `formula`
#'   references (a warning is issued if it references any and
#'   `covariates` is `NULL`).
#' @param use_draw For Bayesian fits only: an integer index into the stored
#'   posterior draws, or `"random"` to sample a fresh draw. If supplied,
#'   that draw's parameters are used instead of the posterior mean - call
#'   this repeatedly (e.g. with `"random"`) to propagate full posterior
#'   uncertainty into a simulation ensemble. `NULL` (default) uses the
#'   point estimate (posterior mean for Bayesian fits, MLE for ML fits).
#' @param ... Unused.
#'
#' @return A `data.table`, as for [dinarch_simulate()].
#' @export
simulate.dinarch_fit <- function(object, nsim = 10, seed = NULL,
                                  covariates = NULL,
                                  y_threshold = NULL,
                                  y_init = NULL, burn_in = 0,
                                  use_draw = NULL, ...) {

  pars <- .dinarch_pick_params(object, use_draw)
  b <- pars$b; phi <- pars$phi; beta <- pars$beta

  formula_terms <- attr(stats::terms(object$formula), "term.labels")
  if (length(formula_terms) > 0 && is.null(covariates)) {
    warning(
      "This fit's `formula` references covariate(s) (",
      paste(formula_terms, collapse = ", "),
      "); no forward `covariates` were supplied, so `model.matrix()` will ",
      "fail unless those variables happen to not be needed (e.g. an ",
      "interaction-only formula) - supply a `covariates` data.frame with ",
      "the forward paths for these variables."
    )
  }

  dinarch_simulate(
    n = nsim,
    b = b, phi = phi,
    formula = object$formula,
    covariates = covariates,
    beta = beta,
    y_init = y_init,
    burn_in = burn_in,
    y_threshold = y_threshold,
    seed = seed
  )
}
