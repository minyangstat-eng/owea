# Convergence warnings, the multistage LOCAL-convergence warning, and check_global.

lin_info <- function(x) tcrossprod(c(1, x[1]))

test_that("multistage emits a LOCAL-convergence warning", {
  expect_warning(
    optimal_design(info_matrix = lin_info, design_box = list(c(-1, 1)),
                   step_sequence = c(0.5, 0.25), p = 0),
    "LOCAL")
})

test_that("check_global verifies the design over the whole box", {
  res <- optimal_design(info_matrix = lin_info, design_box = list(c(-1, 1)),
                        step_sequence = c(0.5, 0.25, 0.1), p = 0,
                        check_global = TRUE)
  expect_false(is.na(res$global_max_d))
  expect_true(isTRUE(res$global_check))      # linear D-opt is globally optimal
  expect_lte(res$global_max_d, 1e-6)
})

test_that("owea warns when the design does not converge", {
  prob <- DesignProblem(X = matrix(seq(0, 3, length.out = 601), ncol = 1),
                        theta = theta_biexp, p = 0, info_matrix = biexp_info)
  expect_warning(r <- owea(prob, max_outer = 1L), "did NOT converge")
  expect_false(r$converged)
})

test_that("candidate_set carries global_check (whole set is the domain)", {
  X <- make_grid(0, 3, 0.05)
  res <- optimal_design(info_matrix = biexp_info, theta = theta_biexp,
                        candidate_set = X, p = 0)
  expect_true(res$converged)
  expect_true(isTRUE(res$global_check))
})
