# exact_design(): integer allocation built on the approximate optimum.

test_that("linear regression D-optimal: even n splits evenly, efficiency 1", {
  lin_info <- function(x) tcrossprod(c(1, x[1]))
  res <- exact_design(n = 10, info_matrix = lin_info,
                      design_box = list(c(-1, 1)),
                      step_sequence = c(0.1, 0.05), p = 0, seed = 1)
  expect_equal(sum(res$counts), 10L)
  sp <- sort(as.numeric(res$support))
  expect_equal(length(sp), 2L)
  expect_lt(abs(sp[1] + 1), 1e-6)
  expect_lt(abs(sp[2] - 1), 1e-6)
  expect_equal(sort(res$counts), c(5L, 5L))
  expect_gt(res$efficiency, 1 - 1e-6)         # exact optimum reached
  expect_lte(res$efficiency, 1 + 1e-8)
})

test_that("bi-exponential D-optimal: counts sum to n, efficiency in (0,1]", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  for (N in c(8L, 20L, 100L)) {
    res <- exact_design(n = N, candidate_set = X, theta = theta_biexp,
                        info_matrix = biexp_info, p = 0, seed = 42)
    expect_equal(sum(res$counts), N)
    expect_true(all(res$counts >= 1L))
    expect_true(is.finite(res$criterion))
    expect_gt(res$efficiency, 0)
    expect_lte(res$efficiency, 1 + 1e-8)
  }
})

test_that("efficiency improves (weakly) as n grows", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  e_small <- exact_design(n = 8L,  candidate_set = X, theta = theta_biexp,
                          info_matrix = biexp_info, p = 0, seed = 7)$efficiency
  e_large <- exact_design(n = 200L, candidate_set = X, theta = theta_biexp,
                          info_matrix = biexp_info, p = 0, seed = 7)$efficiency
  expect_gt(e_large, e_small - 1e-6)
  expect_gt(e_large, 0.98)                     # large n -> near the optimum
})

test_that("A-optimal subset path works", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  res <- exact_design(n = 30L, candidate_set = X, theta = theta_biexp,
                      info_matrix = biexp_info, p = 1, subset = c(2, 4),
                      seed = 3)
  expect_equal(sum(res$counts), 30L)
  expect_gt(res$efficiency, 0)
  expect_lte(res$efficiency, 1 + 1e-8)
})

test_that("multistage (existing design) exact allocation of the new stage", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  res <- exact_design(n = 40L, candidate_set = X, theta = theta_biexp,
                      info_matrix = biexp_info, p = 0,
                      xi0_points = matrix(c(0, 1, 2, 3), ncol = 1),
                      xi0_weights = rep(0.25, 4), n0 = 40, n1 = 80, seed = 1)
  expect_equal(sum(res$counts), 40L)
  expect_true(is.finite(res$criterion))
  expect_lte(res$efficiency, 1 + 1e-8)
})

test_that("same seed is reproducible; vector input matches matrix input", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  a <- exact_design(n = 25L, candidate_set = X, theta = theta_biexp,
                    info_matrix = biexp_info, p = 0, seed = 99)
  b <- exact_design(n = 25L, candidate_set = X, theta = theta_biexp,
                    info_matrix = biexp_info, p = 0, seed = 99)
  expect_equal(a$counts, b$counts)
  expect_equal(a$support, b$support)

  v <- exact_design(n = 25L, candidate_set = X, theta = theta_biexp,
                    info_vector = biexp_vec, p = 0, seed = 99)
  expect_lt(abs(v$criterion - a$criterion), 1e-6)
})

test_that("print_result() dispatches to the exact-design printout", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  res <- exact_design(n = 25L, candidate_set = X, theta = theta_biexp,
                      info_matrix = biexp_info, p = 0, seed = 1)
  out <- capture.output(print_result(res))
  expect_true(any(grepl("Exact design", out)))
  expect_true(any(grepl("efficiency", out)))
  expect_true(any(grepl("count", out)))
})

test_that("n below the minimum support errors", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  expect_error(
    exact_design(n = 3L, candidate_set = X, theta = theta_biexp,
                 info_matrix = biexp_info, p = 0),
    "minimum support")
})
