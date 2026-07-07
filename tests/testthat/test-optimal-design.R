# optimal_design() multistage / partial-parameter / robustness.

test_that("multistage logistic GLM A-optimal converges", {
  res <- suppressWarnings(optimal_design(info_matrix = logistic_info,
                        design_box = list(c(-2,2), c(-1,1), c(-3,3)),
                        step_sequence = c(0.2, 0.1, 0.05, 0.02, 0.01),
                        p = 1))
  expect_true(res$converged)
  expect_lte(res$max_d, 1e-6)
  expect_lt(abs(sum(res$weights) - 1), 1e-8)
  expect_equal(ncol(as.matrix(res$support)), 3L)
  expect_true(is.finite(res$criterion))
})

test_that("subset == grad_g in optimal_design()", {
  seq5 <- c(0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
  res_sub <- suppressWarnings(optimal_design(info_matrix = biexp_info, design_box = list(c(0, 3)),
                            step_sequence = seq5, p = 1, theta = theta_biexp,
                            subset = c(2, 4)))
  gg <- function(th) matrix(c(0,1,0,0, 0,0,0,1), nrow = 2, byrow = TRUE)
  res_gg <- suppressWarnings(optimal_design(info_matrix = biexp_info, design_box = list(c(0, 3)),
                           step_sequence = seq5, p = 1, theta = theta_biexp,
                           grad_g = gg))
  expect_true(res_sub$converged)
  expect_lte(res_sub$max_d, 1e-6)
  expect_lt(abs(res_sub$criterion - res_gg$criterion), 1e-6)
})

test_that("singular subset problem does not crash (pseudo-inverse path)", {
  hard <- suppressWarnings(try(owea(DesignProblem(X = make_grid(c(-2,-1,-3), c(2,1,3), 0.5),
                                 theta = c(1, -0.5, 0.5, 1), p = 1,
                                 info_matrix = logistic_info, subset = c(2, 4)),
                   max_outer = 40L), silent = TRUE))
  expect_false(inherits(hard, "try-error"))
  expect_true(is.finite(hard$criterion))
})

test_that("init_method strategies reach the same optimum", {
  Xg <- make_grid(c(-2, -1, -3), c(2, 1, 3), 0.4)
  crit <- vapply(c("minmax", "minmaxmedian", "random"), function(m) {
    r <- owea(DesignProblem(X = Xg, theta = c(1,-0.5,0.5,1), p = 0,
                            info_matrix = logistic_info), init_method = m)
    c(conv = as.numeric(isTRUE(r$converged)), crit = r$criterion)
  }, numeric(2))
  expect_true(all(crit["conv", ] == 1))
  expect_lt(diff(range(crit["crit", ])), 1e-6)
  expect_error(owea(DesignProblem(X = Xg, theta = c(1,-0.5,0.5,1), p = 0,
                                  info_matrix = logistic_info),
                    init_method = "bogus"))
})
