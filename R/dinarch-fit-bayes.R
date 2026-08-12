# Bayesian estimation (Stan) -----------------------------------------------

.dinarch_stan_code <- "
data {
  int<lower=1> N;
  int<lower=0> y[N];
  vector<lower=0>[N] population;
  int<lower=1> k;
  matrix[N, k] ylags;
  int<lower=0> p;
  matrix[N, p] Xcov;

  real prior_a_mean;
  real<lower=0> prior_a_sd;
  real<lower=0> prior_b_alpha;
  real<lower=0> prior_b_beta;
  real<lower=0> prior_phi_shape;
  real<lower=0> prior_phi_rate;
  real prior_beta_mean;
  real<lower=0> prior_beta_sd;
}
parameters {
  real a;
  vector<lower=0, upper=1>[k] b;
  real<lower=0> phi;
  vector[p] beta;
}
transformed parameters {
  vector[N] mu;
  vector[N] r;
  for (n in 1:N) {
    mu[n] = exp(a) * population[n];
    for (j in 1:k) {
      mu[n] += b[j] * ylags[n, j];
    }
    if (p > 0) {
      for (l in 1:p) {
        mu[n] += exp(beta[l]) * Xcov[n, l];
      }
    }
    r[n] = phi * mu[n];
  }
}
model {
  a ~ normal(prior_a_mean, prior_a_sd);
  b ~ beta(prior_b_alpha, prior_b_beta);
  phi ~ gamma(prior_phi_shape, prior_phi_rate);
  if (p > 0) {
    beta ~ normal(prior_beta_mean, prior_beta_sd);
  }
  for (n in 1:N) {
    // Stan's neg_binomial(alpha, beta) has beta = odds = prob/(1-prob).
    // With p = phi/(phi+1), that odds is exactly phi - NOT phi/(phi+1).
    y[n] ~ neg_binomial(r[n], phi);
  }
}
"

.dinarch_env <- new.env(parent = emptyenv())

#' Get (and lazily compile + cache) the DINARCH Stan model
#' @noRd
.get_dinarch_stan_model <- function() {
  if (is.null(.dinarch_env$stan_model)) {
    .dinarch_env$stan_model <- rstan::stan_model(model_code = .dinarch_stan_code)
  }
  .dinarch_env$stan_model
}

#' Fit a DINARCH model by Bayesian estimation (Stan)
#'
#' Same model as [dinarch_fit_ml()], estimated via Hamiltonian Monte Carlo
#' with [rstan::sampling()]. Requires the (Suggested) `rstan` package.
#'
#' @inheritParams dinarch_fit_ml
#' @param prior A list with elements `a` (mean, sd for a Normal prior),
#'   `b` (shape1, shape2 for a Beta prior, applied to every lag
#'   coefficient), `phi` (shape, rate for a Gamma prior), and `beta`
#'   (mean, sd for a Normal prior, only used if `covariates` is supplied).
#'   Any omitted elements fall back to the defaults.
#' @param iter,chains Passed to [rstan::sampling()].
#' @param seed Optional random seed, passed to [rstan::sampling()].
#' @param ... Further arguments passed to [rstan::sampling()].
#'
#' @return An object of class `"dinarch_fit"`, as for [dinarch_fit_ml()],
#'   with `method = "bayes"`, plus `stanfit` (the raw `stanfit` object) and
#'   `posterior` (a list of posterior draws for `a`, `b`, `phi`, `beta`).
#'
#' @export
dinarch_fit_bayes <- function(data,
                               y = "y", index = "index", group = NULL,
                               population = NULL, covariates = NULL,
                               k = 1,
                               prior = list(),
                               iter = 5000, chains = 3, seed = NULL, ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for dinarch_fit_bayes() but is not installed.")
  }
  stopifnot(is.numeric(k), length(k) == 1, k >= 1, k == round(k))

  default_prior <- list(a = c(0, 3), b = c(1, 3), phi = c(1, 10), beta = c(0, 3))
  prior <- utils::modifyList(default_prior, prior)

  prep <- .prepare_dinarch_data(data, y, index, group, population, covariates, k)
  Yc <- prep$Yc
  n_beta <- prep$n_beta

  stan_data <- list(
    N = prep$n,
    y = prep$y_vec,
    population = prep$pop_vec,
    k = k,
    ylags = prep$lag_mat,
    p = n_beta,
    Xcov = prep$Xcov,
    prior_a_mean = prior$a[1], prior_a_sd = prior$a[2],
    prior_b_alpha = prior$b[1], prior_b_beta = prior$b[2],
    prior_phi_shape = prior$phi[1], prior_phi_rate = prior$phi[2],
    prior_beta_mean = prior$beta[1], prior_beta_sd = prior$beta[2]
  )

  stan_model_obj <- .get_dinarch_stan_model()

  sampling_args <- c(
    list(object = stan_model_obj, data = stan_data, iter = iter, chains = chains),
    if (!is.null(seed)) list(seed = seed),
    list(...)
  )
  stanfit <- do.call(rstan::sampling, sampling_args)

  # --- convergence diagnostic --------------------------------------------
  summary_pars <- c("a", "b", "phi", if (n_beta > 0) "beta")
  post_summary <- rstan::summary(stanfit, pars = summary_pars)$summary
  max_rhat <- suppressWarnings(max(post_summary[, "Rhat"], na.rm = TRUE))
  if (is.finite(max_rhat) && max_rhat > 1.1) {
    warning(sprintf(
      "Maximum Rhat = %.3f across monitored parameters exceeds 1.1; consider more iterations/chains.",
      max_rhat
    ))
  }

  # --- posterior means, via extract() (robust to k = 1 / p = 0 edge cases) -
  draws <- rstan::extract(stanfit, pars = c("a", "b", "phi", if (n_beta > 0) "beta"))

  a_hat   <- mean(draws$a)
  b_hat   <- unname(colMeans(draws$b))
  phi_hat <- mean(draws$phi)
  beta_hat <- if (n_beta > 0) unname(colMeans(draws$beta)) else numeric(0)

  structure(
    list(
      method = "bayes",
      coefficients = list(a = a_hat, b = b_hat, phi = phi_hat, beta = beta_hat),
      k = k,
      population_used = !is.null(population),
      covariates_used = covariates,
      loglik = NA_real_,
      n_obs = nrow(Yc),
      convergence = max_rhat,
      stanfit = stanfit,
      posterior = list(a = draws$a, b = draws$b, phi = draws$phi,
                        beta = if (n_beta > 0) draws$beta else numeric(0)),
      call = match.call()
    ),
    class = "dinarch_fit"
  )
}
