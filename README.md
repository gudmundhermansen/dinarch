# dinarch

Dynamic Innovation Negative Binomial ARCH (NB-DINARCH) models for count
time series: fitting (maximum likelihood or Bayesian via Stan),
simulation, and forecasting. Intended for applications such as conflict
research and event-count modelling.

The model: for a series `y_t` with `n_lags` autoregressive lags and a
covariate linear predictor `eta_t` (specified via a model formula),

```
eta_t = beta_0 + beta_1 * z_1t + ...
mu_t  = exp(eta_t) + sum_j b_j * y_{t-j}
y_t  ~ NegBinom(mean = mu_t, dispersion = phi)
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("gudmundhermansen/dinarch")
```

Bayesian estimation (`dinarch_fit_bayes()`) additionally requires
[`rstan`](https://mc-stan.org/rstan/).

## Quick example

```r
library(dinarch)

# simulate a count series with one covariate
set.seed(2)
n <- 200
gdp <- rnorm(n)
sim <- dinarch_simulate(
  n = n, b = 0.3, phi = 8,
  formula = ~gdp, covariates = data.frame(gdp = gdp),
  beta = c(log(4), 0.5),
  seed = 2
)
dat <- data.frame(y = sim$y, index = sim$index, gdp = gdp)

# fit
fit <- dinarch_fit_ml(dat, y = "y", index = "index", formula = ~gdp, n_lags = 1)
summary(fit)

# forecast 10 periods ahead, 200 stochastic replicates
forecast <- dinarch_project(fit, horizon = 10, nsim = 200, seed = 2)
```

See `vignette("dinarch")` for a short walkthrough.

## Status

Early-stage and under active development; the API may still change.

## License

MIT
