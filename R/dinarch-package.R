#' dinarch: Dynamic Innovation Negative Binomial ARCH Models
#'
#' Tools for fitting, simulating, and analysing Dynamic Innovation Negative
#' Binomial ARCH (NB-DINARCH) models for count time series. Supports maximum
#' likelihood estimation, optional Bayesian estimation via Stan, optional
#' population scaling and additional covariates, dynamic (time-varying)
#' parameters, and scenario-based simulation.
#'
#' @section Model:
#' For a series `y_t` with `k` autoregressive lags, optional population
#' scaling, and optional additional covariates:
#'
#' \deqn{\mu_t = \exp(a_t) \cdot population_t +
#'   \sum_{j=1}^k b_{t,j} \, y_{t-j} +
#'   \sum_l \exp(\beta_l) \cdot covariate_{t,l}}
#' \deqn{r_t = \phi_t \, \mu_t, \quad p_t = \phi_t / (1 + \phi_t)}
#' \deqn{y_t \sim \mathrm{NegBinom}(r_t, p_t)}
#'
#' If `population` is not supplied, `population_t \equiv 1`, recovering a
#' pure-intercept baseline. If `covariates` is not supplied, that term is
#' dropped entirely.
#'
#' @keywords internal
#' @importFrom data.table data.table setDT setorder is.data.table shift :=
#' @importFrom stats optim dnbinom rnbinom plogis qlogis simulate
"_PACKAGE"

## Columns referenced via data.table's non-standard evaluation get added to
## this list as they're introduced in fit/simulate code, to keep
## R CMD check quiet about "no visible binding for global variable".
## "y" is used bare in utils.R's .add_lag_columns(); "population" is used
## bare as a := assignment target in dinarch_fit_ml()/dinarch_fit_bayes().
utils::globalVariables(c("y", "population"))
