# verify_optimality(): equivalence-theorem check of a GIVEN approximate design.

lin_info <- function(x) tcrossprod(c(1, x[1]))

test_that("confirms a D-optimal linear design (box path)", {
  res <- suppressWarnings(optimal_design(info_matrix = lin_info,
    design_box = list(dose = c(-1, 1)), step_sequence = c(0.1, 0.05), p = 0))
  v <- verify_optimality(res$support, res$weights, info_matrix = lin_info,
                         design_box = list(dose = c(-1, 1)), step = 0.05, p = 0)
  expect_true(v$is_optimal$value)
  expect_match(v$is_optimal$note, "tol")            # note lives under is_optimal
  expect_lt(v$max_sensitivity, 1e-4)
  expect_equal(v$criterion, res$criterion, tolerance = 1e-5)
  expect_equal(unname(v$information), unname(res$information), tolerance = 1e-5)
  expect_equal(dim(v$information), c(2L, 2L))
})

test_that("A-optimal design and design-result input both verify", {
  res <- suppressWarnings(optimal_design(info_matrix = lin_info,
    design_box = list(c(-1, 1)), step_sequence = c(0.1, 0.05), p = 1))
  v <- verify_optimality(res, info_matrix = lin_info,      # pass the result object
                         design_box = list(c(-1, 1)), step = 0.05, p = 1)
  expect_true(v$is_optimal$value)
  expect_equal(v$criterion, res$criterion, tolerance = 1e-5)
})

test_that("flags a non-optimal design (perturbed weights)", {
  v <- verify_optimality(matrix(c(-1, 1), 2, 1), c(0.7, 0.3),
                         info_matrix = lin_info, design_box = list(c(-1, 1)),
                         step = 0.05, p = 0)
  expect_false(v$is_optimal$value)
  expect_gt(v$max_sensitivity, 1e-3)
})

test_that("tol controls is_optimal$value and is echoed in the note", {
  v <- verify_optimality(matrix(c(-1, 1), 2, 1), c(0.5, 0.5), info_matrix = lin_info,
                         design_box = list(c(-1, 1)), step = 0.05, p = 0, tol = 1e-3)
  expect_true(v$is_optimal$value)
  expect_match(v$is_optimal$note, "0.001")
})

test_that("weights are validated as hard errors", {
  expect_error(verify_optimality(matrix(c(-1, 1), 2, 1), c(0.6, 0.6),
    info_matrix = lin_info, design_box = list(c(-1, 1)), step = 0.05, p = 0),
    "sum to 1")
  expect_error(verify_optimality(matrix(c(-1, 1), 2, 1), c(1.2, -0.2),
    info_matrix = lin_info, design_box = list(c(-1, 1)), step = 0.05, p = 0),
    "nonnegative")
})

test_that("support points off the grid warn but still return a result", {
  # -0.97 is not a multiple-of-0.05 grid point on [-1, 1]
  expect_warning(
    v <- verify_optimality(matrix(c(-0.97, 1), 2, 1), c(0.5, 0.5),
      info_matrix = lin_info, design_box = list(c(-1, 1)), step = 0.05, p = 0),
    "not among the design points")
  expect_type(v$is_optimal$value, "logical")
  expect_true(is.finite(v$max_sensitivity))
})

test_that(".offgrid_dmax matches the full grid scan, coordinate by coordinate", {
  box <- list(A = 3, x = c(-1, 1), z = c(0, 10))
  m   <- owea:::.parse_design_box(box)
  by  <- owea:::.expand_stage_step(c(0.25, 2), m$is_factor,
                                   sum(!m$is_factor), length(box))
  pts <- rbind(c(1, -0.25, 4),     # exactly on the grid
               c(3, 0.13, 4),      # x between grid points (nearest 0.25)
               c(2, 1.4, 11))      # outside the box in x and z
  fast <- owea:::.offgrid_dmax(pts, m$lo, m$hi, by, m$is_factor, m$nlevels)

  # reference: infinity-norm distance to the nearest row of the BUILT grid
  X    <- candidate_grid(box, c(0.25, 2))
  idx  <- owea:::.nearest_idx(X, pts)
  slow <- vapply(seq_len(nrow(pts)),
                 function(i) max(abs(pts[i, ] - X[idx[i], ])), numeric(1))
  expect_equal(fast, slow)
  expect_equal(fast[1], 0)                       # on-grid point has zero deviation
  expect_gt(fast[2], 1e-6)                       # off-grid point is flagged
  expect_gt(fast[3], 1e-6)                       # out-of-box point is flagged
})

test_that("candidate_set path verifies and warns on non-members", {
  X <- matrix(seq(-1, 1, by = 0.1), ncol = 1)
  res <- suppressWarnings(optimal_design(info_matrix = lin_info,
                                         candidate_set = X, p = 0))
  v <- verify_optimality(res$support, res$weights, info_matrix = lin_info,
                         candidate_set = X, p = 0)
  expect_true(v$is_optimal$value)
  # a support point that is not one of the candidate points -> warning, not error
  expect_warning(
    v2 <- verify_optimality(matrix(c(-0.95, 1), 2, 1), c(0.5, 0.5),
      info_matrix = lin_info, candidate_set = X, p = 0),
    "not among the candidate_set")
  expect_type(v2$is_optimal$value, "logical")
})

test_that("max_points caps the design_box grid (non-interactive -> error)", {
  # design_box = [0,1] with step 0.001 -> 1001 points > max_points = 100
  expect_error(verify_optimality(matrix(0, 1, 1), 1, info_matrix = lin_info,
    design_box = list(c(0, 1)), step = 0.001, p = 0, max_points = 100),
    "max_points")
  # the message offers the criterion-only escape hatch
  expect_error(verify_optimality(matrix(0, 1, 1), 1, info_matrix = lin_info,
    design_box = list(c(0, 1)), step = 0.001, p = 0, max_points = 100),
    "criterion_only")
  # a modest grid under the default cap runs fine
  expect_error(verify_optimality(matrix(c(-1, 1), 2, 1), c(0.5, 0.5),
    info_matrix = lin_info, design_box = list(c(-1, 1)), step = 0.05, p = 0), NA)
})

test_that("criterion_only returns the criterion without the sensitivity scan", {
  sup <- matrix(c(-1, 1), 2, 1); w <- c(0.5, 0.5)
  full <- verify_optimality(sup, w, info_matrix = lin_info,
                            design_box = list(c(-1, 1)), step = 0.05, p = 0)
  co   <- verify_optimality(sup, w, info_matrix = lin_info,
                            design_box = list(c(-1, 1)), step = 0.05, p = 0,
                            criterion_only = TRUE)
  expect_equal(co$criterion, full$criterion)           # same criterion value
  expect_equal(co$information, full$information)       # same information matrix
  expect_identical(co$max_sensitivity, NA_real_)       # no sensitivity computed
  expect_identical(co$is_optimal$value, NA)            # optimality NOT assessed
  expect_match(co$is_optimal$note, "criterion_only")
  expect_null(co$maximiser)
})

test_that("criterion_only never builds the grid, so any grid size works", {
  # ~1e9 grid points with max_points = 100: full verify would refuse, but
  # criterion_only runs instantly and silently
  expect_warning(
    v <- verify_optimality(matrix(c(0, 1), 2, 1), c(0.5, 0.5),
      info_matrix = lin_info, design_box = list(c(0, 1)), step = 1e-9, p = 0,
      max_points = 100, criterion_only = TRUE),
    NA)
  expect_true(is.finite(v$criterion))

  # the design is STILL validated: bad weights are hard errors ...
  expect_error(verify_optimality(matrix(c(0, 1), 2, 1), c(0.6, 0.6),
    info_matrix = lin_info, design_box = list(c(0, 1)), step = 1e-9, p = 0,
    criterion_only = TRUE), "sum to 1")
  # ... and off-grid support points still warn (checked without the grid)
  expect_warning(verify_optimality(matrix(c(0.25, 1), 2, 1), c(0.5, 0.5),
    info_matrix = lin_info, design_box = list(c(0, 1)), step = 0.4, p = 0,
    criterion_only = TRUE), "not among the design points")
})

test_that("works with a formula-style model spec (logit)", {
  res <- suppressWarnings(optimal_design(design_box = list(dose = c(-1, 1)),
    step_sequence = c(0.2, 0.1), link = "logit", x = 1,
    theta = c(0.2, 1.0), p = 0))
  v <- suppressWarnings(verify_optimality(res$support, res$weights,
    design_box = list(dose = c(-1, 1)), step = 0.05, link = "logit", x = 1,
    theta = c(0.2, 1.0), p = 0))
  expect_true(v$is_optimal$value)
  expect_equal(rownames(v$information), c("(Intercept)", "dose"))
})
