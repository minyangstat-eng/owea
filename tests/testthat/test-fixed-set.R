# owea() on a fixed candidate set (ported from the original _selftest.R).

X <- matrix(seq(0, 3, length.out = 601), ncol = 1)

test_that("D-optimal nonlinear bi-exponential converges", {
  prob <- DesignProblem(X = X, theta = theta_biexp, p = 0, info_matrix = biexp_info)
  res  <- owea(prob)
  expect_true(res$converged)
  expect_lte(res$max_d, 1e-6)
  expect_lt(abs(sum(res$weights) - 1), 1e-8)
  expect_true(all(res$weights > 0))
  # On a discrete grid the optimum can be split across adjacent cells, so the
  # raw support has at least k = 4 points; merging collapses it to exactly k.
  expect_gte(nrow(res$support), 4L)
  res_m <- owea(prob, merge = TRUE, merge_atol = 0.05)
  expect_equal(nrow(res_m$support), 4L)
  expect_lt(abs(res_m$criterion - res$criterion), 1e-4)
})

test_that("A-optimal subset == explicit grad_g", {
  prob_sub <- DesignProblem(X = X, theta = theta_biexp, p = 1,
                            info_matrix = biexp_info, subset = c(2, 4))
  res_sub  <- owea(prob_sub)
  expect_true(res_sub$converged)
  expect_lte(res_sub$max_d, 1e-6)

  gg <- function(th) matrix(c(0,1,0,0, 0,0,0,1), nrow = 2, byrow = TRUE)
  res_gg <- owea(DesignProblem(X = X, theta = theta_biexp, p = 1,
                               info_matrix = biexp_info, grad_g = gg))
  expect_lt(abs(res_sub$criterion - res_gg$criterion), 1e-7)
})

test_that("two-stage design converges (existing design)", {
  prob <- DesignProblem(X = X, theta = theta_biexp, p = 0, info_matrix = biexp_info,
                        xi0_points = matrix(c(0,1,2,3), ncol = 1),
                        xi0_weights = rep(0.25, 4), n0 = 40, n1 = 80)
  res <- owea(prob)
  expect_true(res$converged)
  expect_lte(res$max_d, 1e-6)
})

test_that("vector input on a fixed set matches matrix input", {
  prob_v <- DesignProblem(X = X, theta = theta_biexp, p = 0, info_vector = biexp_vec)
  prob_m <- DesignProblem(X = X, theta = theta_biexp, p = 0, info_matrix = biexp_info)
  rv <- owea(prob_v); rm_ <- owea(prob_m)
  expect_lt(abs(rv$criterion - rm_$criterion), 1e-6)
})

test_that("known D-optimal design for simple linear regression", {
  lin_info <- function(x) tcrossprod(c(1, x[1]))
  Xlin <- matrix(seq(-1, 1, length.out = 401), ncol = 1)
  res  <- owea(DesignProblem(X = Xlin, theta = c(0, 0), p = 0, info_matrix = lin_info))
  sp   <- sort(as.numeric(res$support))
  expect_equal(length(sp), 2L)
  expect_lt(abs(sp[1] + 1), 1e-6)
  expect_lt(abs(sp[2] - 1), 1e-6)
  expect_true(all(abs(res$weights - 0.5) < 1e-6))
})

test_that("rejects bad subset / dual specification", {
  expect_error(DesignProblem(X = matrix(c(-1, 1), ncol = 1), theta = c(0, 0),
                             info_matrix = function(x) tcrossprod(c(1, x[1])),
                             subset = 1L, grad_g = function(th) diag(2)))
  expect_error(DesignProblem(X = matrix(c(-1, 1), ncol = 1), theta = c(0, 0),
                             info_matrix = function(x) tcrossprod(c(1, x[1])),
                             subset = c(1L, 5L)))
})
