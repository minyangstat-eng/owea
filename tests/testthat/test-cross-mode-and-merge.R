# Vector/matrix cross-equivalence, criterion conventions, and merging.

test_that("vector and matrix input give the same optimal design", {
  box  <- list(c(-2, 2), c(-1, 1), c(-3, 3))
  step <- c(0.2, 0.1)
  rv <- suppressWarnings(optimal_design(info_vector = logistic_vec, theta = beta_logit,
                       design_box = box, step_sequence = step, p = 0))
  rm_ <- suppressWarnings(optimal_design(info_matrix = logistic_info,
                        design_box = box, step_sequence = step, p = 0))
  # The optimum (criterion) is mode-invariant; the discrete support may split
  # grid cells differently between the IBOSS (vector) and minmax (matrix)
  # starts, so we compare the certified criterion, not raw support rows.
  expect_true(rv$converged && rm_$converged)
  expect_lt(abs(rv$criterion - rm_$criterion), 1e-6)
})

test_that("normalised criterion equals -log det(M)/k for full-parameter D-opt", {
  X <- matrix(seq(0, 3, length.out = 301), ncol = 1)
  res <- owea(DesignProblem(X = X, theta = theta_biexp, p = 0,
                            info_matrix = biexp_info))
  M <- design_information(res$support, res$weights, biexp_info)
  expect_lt(abs(res$criterion - (-log(det(M)) / 4)), 1e-7)
})

test_that("phi_value matches owea() value (vector path)", {
  x_grid <- seq(-1, 1, by = 0.01)
  IVA <- rbind(1, x_grid)
  res <- appro_opt(pp = 0, wb = diag(2), infor_vec_all = IVA, tol = 1e-8)
  v   <- phi_value(0, res$index, res$weight, IVA, diag(2))
  expect_lt(abs(v - res$value), 1e-9)
})

test_that("merge = FALSE leaves the engine result unchanged; merge default off", {
  expect_identical(formals(owea)$merge, FALSE)
  expect_identical(formals(optimal_design)$merge, FALSE)
  X <- matrix(seq(0, 3, length.out = 601), ncol = 1)
  prob <- DesignProblem(X = X, theta = theta_biexp, p = 0, info_matrix = biexp_info)
  r0 <- owea(prob, merge = FALSE)
  r1 <- owea(prob, merge = TRUE, merge_atol = 1e-2)
  # criterion should be essentially unchanged by merging (designs are equivalent)
  expect_lt(abs(r0$criterion - r1$criterion), 1e-4)
})

test_that("merge_close_points performs weighted-centroid merging", {
  sup <- matrix(c(0, 0.001, 1.0), ncol = 1); wts <- c(0.3, 0.2, 0.5)
  mc <- merge_close_points(sup, wts, atol = 0.01)
  expect_equal(nrow(mc$support), 2L)
  expect_lt(abs(sum(mc$weights) - 1), 1e-12)
})
