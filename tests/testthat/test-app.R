# Tests for the web-app translator helpers (R/app.R): the friendly UI state
# must map to exactly the optimal_design() arguments the demos use.

test_that(".ui_design_box encodes factors and a continuous step sequence", {
  covs <- list(
    list(name = "A", type = "factor", nlevels = 3),
    list(name = "B", type = "continuous", lo = 0, hi = 1, steps = c(0.2, 0.1)))
  db <- owea:::.ui_design_box(covs)
  expect_equal(db$design_box, list(A = 3L, B = c(0, 1)))
  expect_true(is.list(db$step_sequence))          # one continuous covariate
  expect_equal(length(db$step_sequence), 2L)      # two stages (coarsest first)
  expect_equal(db$step_sequence[[1]], 0.2)
  expect_equal(db$step_sequence[[2]], 0.1)        # finest stage
  expect_equal(db$finest, 0.1)
})

test_that(".ui_design_box sorts coarsest-first and pads unequal sequences", {
  covs <- list(
    list(name = "u", type = "continuous", lo = -1, hi = 1, steps = c(0.05, 0.2, 0.1)),
    list(name = "v", type = "continuous", lo = 0, hi = 10, steps = c(2, 1)))
  db <- owea:::.ui_design_box(covs)
  expect_equal(length(db$step_sequence), 3L)      # padded to K = 3
  expect_equal(db$step_sequence[[1]], c(0.2, 2))  # coarsest; v padded with its coarsest
  expect_equal(db$step_sequence[[2]], c(0.1, 2))
  expect_equal(db$step_sequence[[3]], c(0.05, 1)) # finest
  expect_equal(db$finest, c(0.05, 1))
})

test_that(".ui_design_box accepts a scalar step (back-compat) and all-factor box", {
  db <- owea:::.ui_design_box(list(
    list(name = "x", type = "continuous", lo = -1, hi = 1, step = 0.05)))
  expect_equal(db$step_sequence, list(0.05))      # single stage
  expect_equal(db$finest, 0.05)
  db2 <- owea:::.ui_design_box(list(
    list(name = "A", type = "factor", nlevels = 2),
    list(name = "B", type = "factor", nlevels = 3)))
  expect_identical(db2$step_sequence, numeric(0))
  expect_identical(db2$finest, numeric(0))
})

test_that(".ui_grid_sizes counts the grid at each step and flags huge grids", {
  covs <- list(list(name = "A", type = "factor", nlevels = 3),
               list(name = "B", type = "continuous", lo = 0, hi = 10, steps = c(2, 1)))
  gs <- owea:::.ui_grid_sizes(covs)               # 3 * (10/2+1)=18 ; 3 * (10/1+1)=33
  expect_equal(gs, c(18, 33))
  big <- owea:::.ui_grid_sizes(list(
    list(name = "x", type = "continuous", lo = 0, hi = 100, steps = c(1, 0.001))))
  expect_gt(max(big), 1e5)
})

test_that(".ui_design_box rejects bad ranges / levels", {
  expect_error(owea:::.ui_design_box(list(
    list(name = "A", type = "continuous", lo = 1, hi = 0, steps = 0.1))),
    "low < high")
  expect_error(owea:::.ui_design_box(list(
    list(name = "A", type = "factor", nlevels = 1))),
    "at least 2 levels")
})

test_that(".ui_model_spec computes within-kind indices and interaction codes", {
  covs <- list(
    list(name = "dose", type = "factor", nlevels = 2),
    list(name = "conc", type = "continuous", lo = -1, hi = 1, steps = c(0.05)))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)), link = "logit")
  expect_equal(sp$f, 1L)                 # first (only) factor
  expect_equal(sp$x, 1L)                 # first (only) continuous
  expect_equal(sp$fx, list(c(1L, 1L)))   # factor 1 x continuous 1
  expect_null(sp$xx); expect_null(sp$ff)
})

test_that(".ui_model_spec: factor:factor and quadratic terms", {
  covs <- list(list(name = "A", type = "factor", nlevels = 2),
               list(name = "B", type = "factor", nlevels = 3),
               list(name = "z", type = "continuous", lo = -1, hi = 1, steps = c(0.1)))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)),
                              quadratics = 3L, link = "identity")
  expect_equal(sp$f, c(1L, 2L))
  expect_equal(sp$x, 1L)
  expect_equal(sp$ff, list(c(1L, 2L)))
  expect_equal(sp$xx, list(c(1L, 1L)))   # quadratic in the single continuous cov
})

test_that(".ui_coef_names returns labels of the right length/order", {
  covs <- list(
    list(name = "dose", type = "factor", nlevels = 2),
    list(name = "conc", type = "continuous", lo = -1, hi = 1, steps = c(0.05)))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)), link = "logit")
  cn <- owea:::.ui_coef_names(sp)
  expect_equal(cn, c("(Intercept)", "dose", "conc", "dose:conc"))
})

test_that("translator spec reproduces a direct optimal_design() call", {
  covs <- list(
    list(name = "dose", type = "factor", nlevels = 2),
    list(name = "conc", type = "continuous", lo = -1, hi = 1, steps = c(0.2, 0.1, 0.05)))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)), link = "logit")
  th <- c(0.5, -0.8, 1.2, 0.4)
  a <- optimal_design(design_box = sp$design_box, step_sequence = sp$step_sequence,
                      link = sp$link, f = sp$f, x = sp$x, fx = sp$fx,
                      theta = th, p = 1)
  b <- suppressWarnings(optimal_design(
    design_box = list(dose = 2, conc = c(-1, 1)),
    step_sequence = c(0.2, 0.1, 0.05),
    link = "logit", f = 1, x = 1, fx = c(11), theta = th, p = 1))
  ord <- function(r) r$support[do.call(order, as.data.frame(r$support)), , drop = FALSE]
  expect_equal(unname(ord(a)), unname(ord(b)), tolerance = 1e-5)
  expect_equal(a$criterion, b$criterion, tolerance = 1e-5)
})

test_that("plot_design runs for 1, 2 and 3 covariates without error", {
  lin <- function(x) tcrossprod(c(1, x[1]))
  r1 <- suppressWarnings(optimal_design(info_matrix = lin,
                                        design_box = list(d = c(-1, 1)),
                                        step_sequence = c(0.2, 0.1), p = 0))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp); on.exit(unlink(tmp))
  expect_silent(plot_design(r1))
  grDevices::dev.off()
})
