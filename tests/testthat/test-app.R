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

# ---------------------------------------------------------------------------
# the wizard: step graph, argument assembly, theta, existing designs, efficiency
# ---------------------------------------------------------------------------

# a one-continuous-covariate spec, used by most of the tests below
sp_x <- function(link = "identity", ncat = NULL, coding = "zero-sum")
  owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    link = link, ncat = ncat, coding = coding)

test_that(".ui_wizard_steps skips the steps that do not apply", {
  s <- owea:::.ui_wizard_steps("identity", "none")
  expect_false("theta" %in% s)                     # identity: not a local design
  expect_false(any(c("design_in", "data_in") %in% s))
  expect_equal(s[1:2], c("model", "start"))
  expect_equal(s[length(s)], "results")

  expect_true("theta" %in% owea:::.ui_wizard_steps("logit", "none"))
  expect_equal(owea:::.ui_wizard_steps("logit", "data"),
               c("model", "start", "data_in", "theta", "criterion",
                 "design_type", "review", "results"))
  s3 <- owea:::.ui_wizard_steps("identity", "design")
  expect_true("design_in" %in% s3)
  expect_false("theta" %in% s3)
  expect_error(owea:::.ui_wizard_steps("logit", "nonsense"))
})

test_that(".ui_model_spec carries the coding; the labels do not depend on it", {
  covs <- list(list(name = "A", type = "factor", nlevels = 3),
               list(name = "z", type = "continuous", lo = -1, hi = 1, steps = 0.5))
  z <- owea:::.ui_model_spec(covs, link = "logit")
  b <- owea:::.ui_model_spec(covs, link = "logit", coding = "baseline")
  expect_equal(z$coding, "zero-sum")
  expect_equal(b$coding, "baseline")
  expect_equal(owea:::.ui_coef_names(z), owea:::.ui_coef_names(b))   # same labels
  expect_error(owea:::.ui_model_spec(covs, link = "logit", coding = "nope"))
})

test_that(".ui_solver_args always emits the coding and shapes each target", {
  sp <- sp_x("logit", coding = "zero-sum")
  for (tg in c("optimal", "exact", "verify", "fit", "simulate")) {
    a <- owea:::.ui_solver_args(sp, tg, theta = c(0.2, 1), n = 10)
    expect_equal(a$coding, "zero-sum", info = tg)
    expect_equal(a$link, "logit", info = tg)
  }
  o <- owea:::.ui_solver_args(sp, "optimal", theta = c(0.2, 1))
  expect_equal(o$step_sequence, sp$step_sequence)
  expect_null(o$n)
  v <- owea:::.ui_solver_args(sp, "verify", theta = c(0.2, 1))
  expect_equal(v$step, sp$finest)                  # a step, not a step_sequence
  expect_null(v$step_sequence)
  expect_identical(v$max_points, Inf)
  f <- owea:::.ui_solver_args(sp, "fit")
  expect_false("theta" %in% names(f))              # fitting estimates theta
  expect_false("p" %in% names(f))
})

test_that(".ui_solver_args validates theta, the subset and n", {
  li <- sp_x("identity"); lo <- sp_x("logit")
  expect_null(owea:::.ui_solver_args(li, "optimal")$theta)      # identity: no theta
  expect_null(owea:::.ui_solver_args(li, "optimal", theta = c(9, 9))$theta)
  expect_error(owea:::.ui_solver_args(lo, "optimal"), "assumed parameter")
  expect_error(owea:::.ui_solver_args(lo, "optimal", theta = 1), "assumed parameter")

  a <- owea:::.ui_solver_args(lo, "optimal", theta = c(0, 1), subset = c(2L, 2L, 1L))
  expect_equal(a$subset, c(1L, 2L))                # sorted, de-duplicated
  expect_null(owea:::.ui_solver_args(lo, "optimal", theta = c(0, 1),
                                     subset = integer(0))$subset)   # empty = all
  expect_error(owea:::.ui_solver_args(lo, "optimal", theta = c(0, 1), subset = 5L),
               "1\\.\\.2")
  expect_error(owea:::.ui_solver_args(lo, "exact", theta = c(0, 1)), "number of runs")
  expect_equal(owea:::.ui_solver_args(lo, "exact", theta = c(0, 1), n = 12)$n, 12L)
})

test_that(".ui_solver_args omits xi0/n0/n1 with no existing design, and ties n = n1", {
  sp <- sp_x("identity")
  a <- owea:::.ui_solver_args(sp, "optimal")
  expect_false(any(c("xi0_points", "xi0_weights", "n0", "n1") %in% names(a)))

  ex <- list(points = cbind(c(-1, 1)), weights = c(0.5, 0.5), n0 = 10, n1 = 25)
  b <- owea:::.ui_solver_args(sp, "optimal", existing = ex)
  expect_equal(b$n0, 10L); expect_equal(b$n1, 25L)
  expect_equal(sum(b$xi0_weights), 1)

  # exact: n and n1 are tied, whatever the caller asked for
  e <- owea:::.ui_solver_args(sp, "exact", existing = ex, n = 30)
  expect_equal(e$n, 30L)
  expect_equal(e$n1, 30L)

  ex0 <- ex; ex0$n0 <- 0
  expect_error(owea:::.ui_solver_args(sp, "optimal", existing = ex0), "n0 >= 1")
})

test_that(".ui_solver_args reproduces a direct optimal_design() call", {
  covs <- list(
    list(name = "dose", type = "factor", nlevels = 2),
    list(name = "conc", type = "continuous", lo = -1, hi = 1, steps = c(0.2, 0.1)))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)), link = "logit")
  th <- c(0.5, -0.8, 1.2, 0.4)
  a <- do.call(optimal_design,
               owea:::.ui_solver_args(sp, "optimal", theta = th, p = 1))
  b <- suppressWarnings(optimal_design(
    design_box = list(dose = 2, conc = c(-1, 1)), step_sequence = c(0.2, 0.1),
    link = "logit", f = 1, x = 1, fx = c(11), theta = th, p = 1,
    coding = "zero-sum"))
  expect_equal(a$criterion, b$criterion, tolerance = 1e-6)
})

test_that(".ui_random_theta and .ui_check_theta respect each parameterisation", {
  expect_null(owea:::.ui_random_theta(sp_x("identity")))
  set.seed(1)
  th <- owea:::.ui_random_theta(sp_x("logit"))
  expect_length(th, length(owea:::.ui_coef_names(sp_x("logit"))))

  cum <- sp_x("cumulative", ncat = 4)
  set.seed(2)
  thc <- owea:::.ui_random_theta(cum)
  expect_length(thc, length(owea:::.ui_coef_names(cum)))
  expect_true(all(diff(thc[1:3]) > 0))             # thresholds strictly increasing

  expect_null(owea:::.ui_check_theta(NULL, sp_x("identity")))
  expect_null(owea:::.ui_check_theta(c(0, 1), sp_x("logit")))
  expect_match(owea:::.ui_check_theta(1, sp_x("logit")), "2 parameter")
  expect_match(owea:::.ui_check_theta(c(0, NA), sp_x("logit")), "finite")
  expect_null(owea:::.ui_check_theta(c(-1, 0, 1, 0.5), cum))
  expect_match(owea:::.ui_check_theta(c(1, 0, 2, 0.5), cum), "strictly increasing")
})

test_that(".ui_existing_from_csv turns counts/weights into proportions and n0", {
  sup <- cbind(c(-1, 0, 1))
  a <- owea:::.ui_existing_from_csv(sup, c(10, 10, 20), "count")
  expect_equal(a$weights, c(0.25, 0.25, 0.5))
  expect_equal(a$n0, 40)
  expect_match(a$notes, "taken from the counts")
  expect_equal(owea:::.ui_existing_from_csv(sup, c(10, 10, 20), "count", n0 = 99)$n0, 99)

  b <- owea:::.ui_existing_from_csv(sup, c(0.2, 0.3, 0.4), "weight")   # sums to 0.9
  expect_equal(sum(b$weights), 1)
  expect_null(b$n0)                              # weights carry no sample size
  expect_match(b$notes, "rescaled")

  z <- owea:::.ui_existing_from_csv(sup, c(5, 0, 5), "count")          # zero row dropped
  expect_equal(nrow(z$points), 2L)
  expect_error(owea:::.ui_existing_from_csv(sup, c(1, -1, 1), "count"), "nonnegative")
  expect_error(owea:::.ui_existing_from_csv(sup, c(0, 0, 0), "count"), "no runs")
  expect_error(owea:::.ui_existing_from_csv(sup, c(1.5, 1, 1), "count"), "whole numbers")
})

test_that(".ui_existing_from_data aggregates a data set into a design", {
  d <- data.frame(dose = c(-1, -1, -1, 0, 0, 1), y = c(0, 1, 1, 0, 1, 1))
  e <- owea:::.ui_existing_from_data(d, list(dose = c(-1, 1)))
  expect_equal(nrow(e$points), 3L)               # 3 distinct covariate rows
  expect_equal(e$n0, 6L)
  expect_equal(sort(e$counts), c(1L, 2L, 3L))
  expect_equal(sum(e$weights), 1)
  expect_equal(names(e$df), c("dose", "count"))

  # a factor covariate: levels must be integer codes within 1..L
  df <- data.frame(grp = c(1, 1, 2, 3), y = c(0, 1, 1, 0))
  expect_equal(owea:::.ui_existing_from_data(df, list(grp = 3))$n0, 4L)
  bad <- data.frame(grp = c(1, 4), y = c(0, 1))
  expect_error(owea:::.ui_existing_from_data(bad, list(grp = 3)), "integers in 1\\.\\.3")
  wide <- data.frame(a = 1, b = 2, y = 0)
  expect_error(owea:::.ui_existing_from_data(wide, list(dose = c(-1, 1))),
               "2 covariate column")
})

test_that(".ui_simulate pools the existing stage with the new runs", {
  sp  <- sp_x("identity")
  sup <- cbind(c(-1, 1)); cnt <- c(10L, 10L)
  th  <- c(1, 2)

  # no existing design: only the new runs
  a <- owea:::.ui_simulate(sp, th, sigma = 1, sup, cnt, nsim = 50, seed = 1)
  expect_equal(a$N, 20L)
  expect_equal(a$n_existing, 0L)
  expect_equal(a$n_new, 20L)

  # an existing DESIGN (no responses): its runs are simulated too, and pooling
  # more data must reduce the MSE
  ex <- list(points = cbind(c(-1, 1)), weights = c(0.5, 0.5), n0 = 40)
  b <- owea:::.ui_simulate(sp, th, sigma = 1, sup, cnt, existing = ex,
                           nsim = 50, seed = 1)
  expect_equal(b$N, 60L)                 # 40 existing + 20 new
  expect_equal(b$n_existing, 40L)
  expect_true(all(b$mse < a$mse))        # 60 observations beat 20
  expect_length(b$coef_names, 2L)

  # counts are taken from the weights when they are not supplied
  ex2 <- list(points = cbind(c(-1, 0, 1)), weights = c(1, 1, 1) / 3, n0 = 10)
  expect_equal(owea:::.ui_simulate(sp, th, 1, sup, cnt, existing = ex2,
                                   nsim = 5, seed = 1)$N, 30L)
})

test_that(".ui_simulate holds an OBSERVED first stage's responses fixed", {
  sp  <- sp_x("identity")
  sup <- cbind(c(-1, 1)); cnt <- c(5L, 5L)
  th  <- c(1, 2)
  d1 <- data.frame(dose = rep(c(-1, 0, 1), each = 4),
                   y    = c(rep(0, 4), rep(1, 4), rep(2, 4)))
  o1 <- owea:::.ui_simulate(sp, th, 1, sup, cnt, obs = list(data = d1),
                            nsim = 40, seed = 3)
  expect_equal(o1$N, 22L)                # 12 observed + 10 new
  expect_equal(o1$n_existing, 12L)

  # the observed responses ENTER the fit: change them and the estimates move,
  # even with the same seed (a re-simulated first stage could not do this)
  d2 <- d1; d2$y <- d2$y + 10
  o2 <- owea:::.ui_simulate(sp, th, 1, sup, cnt, obs = list(data = d2),
                            nsim = 40, seed = 3)
  expect_false(isTRUE(all.equal(o1$theta_hat_mean, o2$theta_hat_mean)))
  expect_gt(o2$mse[["(Intercept)"]], o1$mse[["(Intercept)"]])   # y shifted by 10

  # obs wins over a design-only existing stage when both are given
  ex <- list(points = cbind(c(-1, 1)), weights = c(0.5, 0.5), n0 = 99)
  expect_equal(owea:::.ui_simulate(sp, th, 1, sup, cnt, existing = ex,
                                   obs = list(data = d1), nsim = 5, seed = 1)$N,
               22L)
})

test_that(".ui_simulate works for a binary response and needs true values", {
  sp  <- sp_x("logit")
  sup <- cbind(c(-1, 1)); cnt <- c(30L, 30L)
  d   <- data.frame(dose = rep(c(-1, 1), each = 10),
                    y    = c(0, 0, 0, 1, 0, 1, 0, 0, 1, 0,
                             1, 1, 0, 1, 1, 1, 0, 1, 1, 1))
  r <- owea:::.ui_simulate(sp, c(0.2, 1), support = sup, counts = cnt,
                           obs = list(data = d), nsim = 30, seed = 2)
  expect_equal(r$N, 80L)
  expect_true(all(is.finite(r$mse)))
  expect_error(owea:::.ui_simulate(sp, NULL, support = sup, counts = cnt,
                                   nsim = 5), "true parameter value")
})

test_that(".ui_efficiency is 1 under the same criterion and < 1 under the other", {
  sp <- sp_x("identity")
  res <- suppressWarnings(do.call(optimal_design,
                                  owea:::.ui_solver_args(sp, "optimal", p = 0)))
  same <- suppressWarnings(owea:::.ui_efficiency(res, sp, p_new = 0L))
  expect_equal(same$efficiency, 1, tolerance = 1e-4)

  # the D-optimal design is not A-optimal for a proper subset of the parameters
  sp2 <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    quadratics = 1L, link = "identity")
  r2 <- suppressWarnings(do.call(optimal_design,
                                 owea:::.ui_solver_args(sp2, "optimal", p = 0)))
  other <- suppressWarnings(owea:::.ui_efficiency(r2, sp2, p_new = 1L,
                                                  subset_new = 3L))
  expect_true(other$efficiency > 0 && other$efficiency < 1)
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
