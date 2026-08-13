#' dinarch: Dynamic Innovation Negative Binomial ARCH Models
#'
#' Tools for fitting, simulating, and analysing Dynamic Innovation Negative
#' Binomial ARCH (NB-DINARCH) models for count time series. Supports maximum
#' likelihood estimation, optional Bayesian estimation via Stan, a formula/
#' `model.matrix()`-based covariate interface, dynamic (time-varying)
#' parameters, and scenario-based simulation.
#'
#' @section Model:
#' For a series `y_t` with `n_lags` autoregressive lags (`k` below) and a
#' covariate linear predictor `eta_t` built from `formula` via
#' [stats::model.matrix()]:
#'
#' \deqn{\eta_t = \sum_j \beta_j z_{t,j}}
#' \deqn{\mu_t = \exp(\eta_t) + \sum_{j=1}^k b_{t,j} \, y_{t-j}}
#' \deqn{r_t = \phi_t \, \mu_t, \quad p_t = \phi_t / (1 + \phi_t)}
#' \deqn{y_t \sim \mathrm{NegBinom}(r_t, p_t)}
#'
#' `formula`'s default `~1` gives an intercept-only `eta_t` (a single
#' `beta_0` playing the role a fixed baseline term might otherwise play).
#' There is no distinguished `population` argument anywhere in the
#' package, including [dinarch_simulate()] - if `mu` should scale
#' (sub/super-)proportionally with a population/exposure variable,
#' include `log(population)` as one of `formula`'s terms, same as fitting.
#' [dinarch_simulate()] separately takes an optional `y_threshold`
#' argument that caps the rolling sum of simulated `y` at a plain numeric
#' limit - deliberately unrelated to `mu`/`formula`/`covariates`, so it
#' can never be confused with (or silently coupled to) a population
#' covariate.
#'
#' @keywords internal
#' @importFrom data.table data.table setDT setorder is.data.table shift :=
#' @importFrom stats optim dnbinom rnbinom plogis qlogis simulate
"_PACKAGE"

## Columns referenced via data.table's non-standard evaluation get added to
## this list as they're introduced in fit/simulate code, to keep
## R CMD check quiet about "no visible binding for global variable".
## "y" is used bare in utils.R's .add_lag_columns().
utils::globalVariables(c("y"))
