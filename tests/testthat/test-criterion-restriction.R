# Only D (p = 0) and A (p = 1) optimality are allowed; every public entry point
# that takes a criterion must reject anything else.

lin_info <- function(x) tcrossprod(c(1, x[1]))
ivall    <- rbind(1, seq(-1, 1, by = 0.1))

test_that("p = 0 and p = 1 are accepted by optimal_design", {
  expect_error(suppressWarnings(optimal_design(info_matrix = lin_info,
    design_box = list(c(-1, 1)), step_sequence = c(0.2, 0.1), p = 0)), NA)
  expect_error(suppressWarnings(optimal_design(info_matrix = lin_info,
    design_box = list(c(-1, 1)), step_sequence = c(0.2, 0.1), p = 1)), NA)
})

test_that("optimal_design rejects p >= 2", {
  expect_error(optimal_design(info_matrix = lin_info, design_box = list(c(-1, 1)),
    step_sequence = c(0.2, 0.1), p = 2), "Only D-optimality")
})

test_that("exact_design rejects p >= 2", {
  expect_error(exact_design(n = 10, info_matrix = lin_info,
    design_box = list(c(-1, 1)), step_sequence = c(0.2, 0.1), p = 3),
    "Only D-optimality")
})

test_that("DesignProblem (and thus owea) rejects p >= 2", {
  X <- candidate_grid(list(c(-1, 1)), 0.1)
  expect_error(DesignProblem(X, info_matrix = lin_info, p = 2), "Only D-optimality")
})

test_that("appro_opt and phi_value reject pp >= 2", {
  expect_error(appro_opt(pp = 2, wb = diag(2), infor_vec_all = ivall),
    "Only D-optimality")
  expect_error(phi_value(pp = 4, index = 1:2, weight = c(0.5, 0.5),
    infor_vec_all = cbind(c(1, -1), c(1, 1)), wb = diag(2)),
    "Only D-optimality")
})

test_that("verify_optimality rejects p >= 2", {
  expect_error(verify_optimality(matrix(c(-1, 1), 2, 1), c(0.5, 0.5),
    info_matrix = lin_info, design_box = list(c(-1, 1)), step = 0.1, p = 2),
    "Only D-optimality")
})

test_that("non-integer / NA criteria are rejected", {
  expect_error(optimal_design(info_matrix = lin_info, design_box = list(c(-1, 1)),
    step_sequence = c(0.2, 0.1), p = NA), "Only D-optimality")
})
