# Scenario construction -----------------------------------------------------

#' Build a smooth parameter transition path
#'
#' Interpolates from `from` to `to` over `n` steps, for constructing
#' scenario parameter paths (e.g. `b` or `phi` trajectories) to feed into
#' [dinarch_simulate()] or [simulate.dinarch_fit()].
#'
#' @param from,to Starting and ending values (scalars). Works the same way
#'   whether `to > from` (increasing path) or `to < from` (decreasing
#'   path).
#' @param n Number of steps in the path.
#' @param shape `"sigmoid"` (default) for a smooth S-curve transition, or
#'   `"linear"` for straight-line interpolation.
#' @param range For `shape = "sigmoid"`, the domain over which the sigmoid
#'   is evaluated before rescaling to `[from, to]`. The default
#'   `c(-3, 9)` is asymmetric - it front-loads the transition, reaching
#'   close to `to` well before the final step. Use a symmetric range like
#'   `c(-6, 6)` for a transition that eases in and out evenly at both
#'   ends.
#'
#' @return A numeric vector of length `n`.
#'
#' @examples
#' # front-loaded transition (the default range)
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20)
#'
#' # symmetric sigmoid ramp
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20, range = c(-6, 6))
#'
#' # straight-line interpolation
#' dinarch_transition_path(from = 0.1, to = 0.5, n = 20, shape = "linear")
#'
#' @export
dinarch_transition_path <- function(from, to, n,
                                     shape = c("sigmoid", "linear"),
                                     range = c(-3, 9)) {
  shape <- match.arg(shape)
  stopifnot(is.numeric(n), length(n) == 1, n >= 1, n == round(n))
  stopifnot(length(from) == 1, length(to) == 1)

  if (shape == "linear") {
    return(seq(from, to, length.out = n))
  }

  x <- seq(range[1], range[2], length.out = n)
  w <- stats::plogis(x)
  from + (to - from) * w
}
