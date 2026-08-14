# Bayesian estimation (Stan) -----------------------------------------------

.dinarch_stan_code <- "
data {
  int<lower=1> N;
  int<lower=0> y[N];
  int<lower=1> n_lags;
  matrix[N, n_lags] ylags;
  int<lower=0> p;
  matrix[N, p] Xcov;

  vector<lower=0>[n_lags + 1] prior_b_alpha;
  real<lower=0> prior_phi_shape;
  real<lower=0> prior_phi_rate;
  vector[p] prior_beta_mean;
  vector<lower=0>[p] prior_beta_sd;
}
parameters {
  // b_full = (b_1, ..., b_n_lags, slack), slack = 1 - sum(b). Constraining
  // this to a simplex guarantees sum(b) < 1 for every posterior draw - the
  // stationarity condition for this model (see dinarch_fit_ml()'s
  // Details and dinarch_simulate()'s stationarity warning) - rather than
  // relying on n_lags independent b[j] ~ beta() priors, which leave
  // sum(b) >= 1 unconstrained and reachable for n_lags > 1.
  simplex[n_lags + 1] b_full;
  real<lower=0> phi;
  vector[p] beta;
}
transformed parameters {
  vector[n_lags] b = head(b_full, n_lags);
  // eta = Xcov * beta, the linear predictor for exp(eta): folds the
  // intercept and all covariates into one term (no separate baseline
  // parameter), matching dinarch_fit_ml(). Guarded by p > 0 since Stan's
  // matrix-vector product isn't meaningful for a zero-column Xcov.
  vector[N] eta = rep_vector(0, N);
  vector[N] mu;
  vector[N] r;
  if (p > 0) {
    eta = Xcov * beta;
  }
  for (n in 1:N) {
    real lag_term = 0;
    for (j in 1:n_lags) {
      lag_term += b[j] * ylags[n, j];
    }
    // mu = exp(eta) + lag_term is > 0 for any real eta/b - no positivity
    // requirement on covariate values.
    mu[n] = exp(eta[n]) + lag_term;
    r[n] = phi * mu[n];
  }
}
model {
  b_full ~ dirichlet(prior_b_alpha);
  phi ~ gamma(prior_phi_shape, prior_phi_rate);
  if (p > 0) {
    // prior_beta_mean/sd are precomputed per-column on the R side (see
    // dinarch_fit_bayes()): a fixed base prior on beta's *standardized*
    // scale, divided by each Xcov column's own empirical sd, so a
    // large-raw-scale covariate automatically gets a correspondingly
    // tight prior instead of the same width regardless of scale.
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
#' @param prior A list with elements `b` (shape1, shape2 - two numbers,
#'   regardless of `n_lags`), `phi` (shape, rate for a Gamma prior), and
#'   `beta` (mean, sd for a Normal prior - see below for how this is
#'   applied). Any omitted elements fall back to the defaults, which are
#'   intended to be non-informative/weakly-informative as a standard
#'   choice, not tuned to any particular dataset.
#'
#'   `b`'s two numbers parameterize a joint Dirichlet prior over
#'   `(b_1, ..., b_{n_lags}, 1 - sum(b))` - `shape1` repeated for each of
#'   the `n_lags` lag components, `shape2` for the `(n_lags + 1)`th
#'   ("slack") component - which guarantees `sum(b) < 1` (the model's
#'   stationarity condition, see [dinarch_fit_ml()]/[dinarch_simulate()])
#'   for every posterior draw, for any `n_lags`. At `n_lags = 1` this is
#'   exactly a `Beta(shape1, shape2)` prior on `b` (a 2-outcome Dirichlet
#'   is a Beta). The default `c(1, 1)` is the flat/uniform prior on the
#'   simplex (every point equally likely) - it still shrinks the prior
#'   mean of each individual lag coefficient towards 0 as `n_lags` grows,
#'   but that's a consequence of more lags sharing the same
#'   `sum(b) < 1` budget under a uniform allocation, not of the specific
#'   shape numbers chosen.
#'
#'   `beta`'s `sd` (default `2.5`, the standard weakly-informative choice
#'   for a log-link GLM coefficient) is treated as the prior width on
#'   *beta's standardized scale*, not applied literally to every raw
#'   coefficient: internally it's divided by each `model.matrix(formula,
#'   data)` column's own empirical standard deviation, so a covariate on
#'   a large raw scale (e.g. population in the millions) automatically
#'   gets a correspondingly tight prior on its (necessarily tiny)
#'   coefficient, instead of the same width regardless of scale (`mean`
#'   is divided by the same factor, so `mean = 0`, the default, stays 0).
#'   The intercept column is constant, so it has no scale to adapt to and
#'   keeps `c(mean, sd)` unscaled.
#' @param iter,chains Passed to [rstan::sampling()].
#' @param seed Optional random seed, passed to [rstan::sampling()].
#' @param ... Further arguments passed to [rstan::sampling()] - e.g. a
#'   custom `init` to override the default described in Details.
#'
#' @details **Covariate scale matters more here than for [dinarch_fit_ml()].**
#'   `mu = exp(eta) + lag term`, so a `formula` covariate on a large raw
#'   scale (e.g. a country's population in the tens of millions) needs a
#'   correspondingly tiny coefficient - fine once the sampler has found
#'   that region, but Stan's default initialization draws every
#'   unconstrained parameter independently from `Uniform(-2, 2)`, so a
#'   coefficient anywhere near that range multiplied by a covariate that
#'   large overflows `exp(eta)` to `Inf`/`0` at essentially every
#'   attempted starting point - `"Initialization failed"`. To avoid this
#'   regardless of covariate scale, `dinarch_fit_bayes()` uses its own
#'   data-driven default `init` (every non-intercept `beta` starts at
#'   exactly 0, so `eta_init` never depends on a covariate's raw scale;
#'   only the intercept and `phi` get a modest data-driven start,
#'   overridable via `init` in `...`). The default `beta` prior (see
#'   `prior` above) is separately scaled per covariate for the same
#'   reason - together these avoid the hard failure, but for the sampler
#'   to actually mix well (and not report a large Rhat), still prefer
#'   transformed/rescaled covariates over raw large-magnitude ones - e.g.
#'   `~ log(population) + gdppc_thousands` rather than
#'   `~ population + gdppc`.
#'
#' @return An object of class `"dinarch_fit"`, as for [dinarch_fit_ml()]
#'   (including `data`, see there), with `method = "bayes"`, plus
#'   `stanfit` (the raw `stanfit` object) and `posterior` (a list of
#'   posterior draws for `b`, `phi`, `beta`).
#'
#' @export
dinarch_fit_bayes <- function(data,
                               y = "y", index = "index", group = NULL,
                               formula = ~1,
                               n_lags = 1,
                               prior = list(),
                               iter = 5000, chains = 3, seed = NULL, ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for dinarch_fit_bayes() but is not installed.")
  }
  formula <- stats::as.formula(formula)
  stopifnot(is.numeric(n_lags), length(n_lags) == 1, n_lags >= 1, n_lags == round(n_lags))

  default_prior <- list(b = c(1, 1), phi = c(1, 0.05), beta = c(0, 2.5))
  prior <- utils::modifyList(default_prior, prior)

  # (b_1, ..., b_n_lags, slack) ~ dirichlet(c(rep(shape1, n_lags), shape2))
  # - see the Stan model and @param prior above for why this enforces
  # sum(b) < 1.
  prior_b_alpha_vec <- c(rep(prior$b[1], n_lags), prior$b[2])

  prep <- .prepare_dinarch_data(data, y, index, group, formula, n_lags)
  Yc <- prep$Yc
  n_beta <- prep$n_beta
  beta_names <- colnames(prep$Xcov)

  # `prior$beta = c(mean, sd)` is the reference prior on beta's
  # *standardized* scale - divide by each Xcov column's own empirical sd
  # to get a prior of comparable width regardless of a covariate's raw
  # scale (mean and sd divide by the same factor, so mean = 0 stays 0
  # unaffected). The intercept column is constant (sd = 0), so it keeps
  # the base prior unscaled - it has no covariate scale to adapt to.
  cov_sd <- if (n_beta > 0) apply(prep$Xcov, 2, stats::sd) else numeric(0)
  scale_factor <- ifelse(cov_sd > 0, cov_sd, 1)
  prior_beta_mean_vec <- prior$beta[1] / scale_factor
  prior_beta_sd_vec <- prior$beta[2] / scale_factor

  stan_data <- list(
    N = prep$n,
    y = prep$y_vec,
    n_lags = n_lags,
    ylags = prep$lag_mat,
    p = n_beta,
    Xcov = prep$Xcov,
    prior_b_alpha = prior_b_alpha_vec,
    prior_phi_shape = prior$phi[1], prior_phi_rate = prior$phi[2],
    # array(): avoids rstan collapsing a length-1 vector to a bare scalar,
    # which Stan then rejects for a declared vector[p] (see the same fix
    # for `beta_init` below).
    prior_beta_mean = array(prior_beta_mean_vec, dim = n_beta),
    prior_beta_sd = array(prior_beta_sd_vec, dim = n_beta)
  )

  stan_model_obj <- .get_dinarch_stan_model()

  # Data-driven default init (one call per chain - rstan calls a function
  # `init` once per chain, so the rnorm() jitters below naturally differ
  # across chains). Every non-intercept `beta` entry starts at exactly 0,
  # regardless of prior/data, so `eta_init` never depends on a covariate's
  # raw scale: Stan's own default ("random" ~ Uniform(-2, 2) per
  # parameter) multiplies that randomness straight into whatever scale the
  # covariate happens to be on - a covariate like raw population (10^6-10^9)
  # makes exp(eta) overflow to Inf/0 at essentially every attempted
  # initialization, which is exactly "Initialization failed" (see Details
  # in the source for a fuller explanation). Only the intercept (if
  # present) and `phi` get a modest jitter, since both are already on a
  # safe, roughly-O(1)-after-transform scale.
  mean_y <- mean(prep$y_vec)
  var_y  <- stats::var(prep$y_vec)
  phi_start <- if (is.finite(var_y) && var_y > mean_y) mean_y^2 / (var_y - mean_y) else 5
  phi_start <- min(max(phi_start, 0.5), 50)
  b_sum_start <- stats::plogis(-1)
  intercept_pos <- which(beta_names == "(Intercept)")

  default_init <- function() {
    beta_init <- rep(0, n_beta)
    if (length(intercept_pos) == 1) {
      beta_init[intercept_pos] <- log(max(mean_y, 1e-6)) + stats::rnorm(1, sd = 0.1)
    }
    list(
      b_full = c(rep(b_sum_start / n_lags, n_lags), 1 - b_sum_start),
      phi = phi_start * exp(stats::rnorm(1, sd = 0.1)),
      # array(): rstan/Stan's data marshalling collapses a plain length-1 R
      # vector to a bare scalar, which Stan then rejects for a declared
      # vector[n_beta] parameter ("dims declared=(1); dims found=()") -
      # array(..., dim = n_beta) keeps the explicit vector-ness for n_beta == 1.
      beta = array(beta_init, dim = n_beta)
    )
  }

  sampling_args <- utils::modifyList(
    list(object = stan_model_obj, data = stan_data, iter = iter, chains = chains, init = default_init),
    c(if (!is.null(seed)) list(seed = seed), list(...))
  )
  stanfit <- do.call(rstan::sampling, sampling_args)

  if (identical(stanfit@mode, 2L)) {
    stop(
      "Stan sampling produced no samples (see the messages above from ",
      "rstan::sampling() for the underlying cause). A common cause is a ",
      "`formula` covariate on a very different numeric scale from the ",
      "others (e.g. raw population in the millions) - try log-transforming ",
      "or otherwise rescaling large-magnitude covariates, or pass a custom ",
      "`init` via `...`."
    )
  }

  # --- convergence diagnostic --------------------------------------------
  summary_pars <- c("b", "phi", if (n_beta > 0) "beta")
  post_summary <- rstan::summary(stanfit, pars = summary_pars)$summary
  max_rhat <- suppressWarnings(max(post_summary[, "Rhat"], na.rm = TRUE))
  if (is.finite(max_rhat) && max_rhat > 1.1) {
    warning(sprintf(
      "Maximum Rhat = %.3f across monitored parameters exceeds 1.1; consider more iterations/chains.",
      max_rhat
    ))
  }

  # --- posterior means, via extract() (robust to n_lags = 1 / p = 0 edge cases)
  draws <- rstan::extract(stanfit, pars = c("b", "phi", if (n_beta > 0) "beta"))

  b_hat   <- unname(colMeans(draws$b))
  phi_hat <- mean(draws$phi)
  beta_hat <- if (n_beta > 0) unname(colMeans(draws$beta)) else numeric(0)
  names(beta_hat) <- beta_names
  if (n_beta > 0) colnames(draws$beta) <- beta_names

  structure(
    list(
      method = "bayes",
      coefficients = list(b = b_hat, phi = phi_hat, beta = beta_hat),
      n_lags = n_lags,
      formula = formula,
      data = prep$full_data,
      loglik = NA_real_,
      n_obs = nrow(Yc),
      convergence = max_rhat,
      stanfit = stanfit,
      posterior = list(b = draws$b, phi = draws$phi,
                        beta = if (n_beta > 0) draws$beta else numeric(0)),
      call = match.call()
    ),
    class = "dinarch_fit"
  )
}
