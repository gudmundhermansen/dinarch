test_that("dinarch_transition_path returns the right length and endpoints", {
  path <- dinarch_transition_path(from = 0, to = 10, n = 25)
  expect_length(path, 25)
  expect_lt(path[1], 1)       # close to `from`, given the default front-loaded range
  expect_gt(path[25], 9.9)    # close to `to`
})

test_that("dinarch_transition_path with shape = 'linear' interpolates exactly", {
  path <- dinarch_transition_path(from = 0, to = 10, n = 5, shape = "linear")
  expect_equal(path, c(0, 2.5, 5, 7.5, 10))
})

test_that("a symmetric sigmoid range hits the midpoint exactly at x = 0", {
  path <- dinarch_transition_path(from = 0, to = 1, n = 3, range = c(-6, 6))
  expect_equal(path[2], 0.5)
  expect_lt(path[1], 0.01)
  expect_gt(path[3], 0.99)
})

test_that("dinarch_transition_path handles a decreasing path (to < from)", {
  path <- dinarch_transition_path(from = 5, to = 1, n = 10)
  expect_true(all(diff(path) <= 0))
  expect_equal(path, dinarch_transition_path(from = 5, to = 1, n = 10))
})

test_that("a transition path plugs directly into dinarch_simulate()'s time-varying a", {
  a_path <- dinarch_transition_path(from = log(1), to = log(10), n = 30)
  sim <- dinarch_simulate(n = 30, a = a_path, b = 0.2, phi = 5, seed = 1)
  expect_equal(nrow(sim), 30)
})
