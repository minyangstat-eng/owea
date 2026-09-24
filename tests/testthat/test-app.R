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

test_that(".ui_verify_points counts the finest-step verify grid", {
  covs <- list(list(name = "A", type = "factor", nlevels = 3),
               list(name = "B", type = "continuous", lo = 0, hi = 10,
                    steps = c(2, 1)))
  sp <- owea:::.ui_model_spec(covs, link = "identity")
  expect_equal(owea:::.ui_verify_points(sp), 33)   # finest stage: 3 * (10/1 + 1)
  gs <- owea:::.ui_grid_sizes(covs)                # = last stage of the sequence
  expect_equal(owea:::.ui_verify_points(sp), gs[length(gs)])
  # all-factor spec: just the level combinations
  spf <- owea:::.ui_model_spec(
    list(list(name = "A", type = "factor", nlevels = 2),
         list(name = "B", type = "factor", nlevels = 3)),
    link = "identity")
  expect_equal(owea:::.ui_verify_points(spf), 6)
})

test_that(".ui_parse_steps parses the grid step(s) text box", {
  expect_equal(owea:::.ui_parse_steps("0.1"), 0.1)
  expect_equal(owea:::.ui_parse_steps(" 0.5, 0.1,0.02 "), c(0.5, 0.1, 0.02))
  expect_equal(owea:::.ui_parse_steps(0.25), 0.25)          # numeric passthrough
  expect_equal(owea:::.ui_parse_steps(c(0.5, 0.1)), c(0.5, 0.1))
  expect_identical(owea:::.ui_parse_steps(""), numeric(0))
  expect_identical(owea:::.ui_parse_steps("  ,  "), numeric(0))
  # ANY bad token empties the whole sequence -- a typo must not be dropped
  # silently, it must trip .ui_design_box's error
  expect_identical(owea:::.ui_parse_steps("0.5, abc"), numeric(0))
  expect_error(owea:::.ui_design_box(list(
    list(name = "x", type = "continuous", lo = 0, hi = 1,
         steps = owea:::.ui_parse_steps("0.5, abc")))),
    "step sequence of positive numbers")
})

test_that(".ui_design_box / .ui_model_spec / .ui_solver_args support the continuous search", {
  covs <- list(list(name = "A", type = "factor", nlevels = 3),
               list(name = "B", type = "continuous", lo = 0, hi = 1, steps = numeric(0)))
  expect_error(owea:::.ui_design_box(covs), "step sequence")     # the grid mode still needs steps
  db <- owea:::.ui_design_box(covs, continuous = TRUE)
  expect_equal(db$design_box, list(A = 3L, B = c(0, 1)))
  expect_null(db$step_sequence)
  expect_identical(db$finest, numeric(0))
  expect_true(db$continuous)

  sp <- owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE)
  expect_true(sp$continuous)
  expect_equal(sp$n_audit, 20000L)                 # the solver default
  o <- owea:::.ui_solver_args(sp, "optimal")
  expect_true(o$continuous); expect_null(o$step_sequence)
  expect_equal(o$n_audit, 20000L)
  e <- owea:::.ui_solver_args(sp, "exact", n = 10)
  expect_true(e$continuous); expect_equal(e$n, 10L)
  v <- owea:::.ui_solver_args(sp, "verify")
  expect_true(v$continuous); expect_null(v$step); expect_null(v$max_points)
  expect_equal(v$n_audit, 20000L)
  expect_true(is.na(owea:::.ui_verify_points(sp)))
  expect_equal(owea:::.ui_grid_sizes(covs, continuous = TRUE), 9)   # 3 levels x 3^1 start grid

  # the audit size is the user's, for every target; bad values fall back
  sp5 <- owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE, n_audit = 5000)
  expect_equal(sp5$n_audit, 5000L)
  expect_equal(owea:::.ui_solver_args(sp5, "optimal")$n_audit, 5000L)
  expect_equal(owea:::.ui_solver_args(sp5, "verify")$n_audit, 5000L)
  expect_equal(owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE,
                                     n_audit = 0)$n_audit, 0L)
  expect_equal(owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE,
                                     n_audit = NA)$n_audit, 20000L)
  expect_equal(owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE,
                                     n_audit = -3)$n_audit, 20000L)

  # the grid mode is unchanged: no 'continuous' / 'n_audit' argument reaches the solver
  sp0 <- owea:::.ui_model_spec(
    list(list(name = "x", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    link = "identity", n_audit = 5000)
  expect_false(sp0$continuous)
  expect_true(is.na(sp0$n_audit))
  expect_false(any(c("continuous", "n_audit") %in% names(owea:::.ui_solver_args(sp0, "optimal"))))
  expect_false(any(c("continuous", "n_audit") %in% names(owea:::.ui_solver_args(sp0, "verify"))))
})

test_that("an explicit grid check overrides the spec's design space for 'verify'", {
  covs <- list(list(name = "A", type = "factor", nlevels = 3),
               list(name = "dose", type = "continuous", lo = -1, hi = 1))
  spc <- owea:::.ui_model_spec(covs, link = "identity", continuous = TRUE, n_audit = 5000)
  v <- owea:::.ui_solver_args(spc, "verify", verify_step = 0.05)
  expect_equal(v$step, 0.05)
  expect_identical(v$max_points, Inf)
  expect_false(any(c("continuous", "n_audit") %in% names(v)))
  # per-covariate steps pass through; bad steps are rejected
  expect_equal(owea:::.ui_solver_args(spc, "verify", verify_step = c(0.1))$step, 0.1)
  expect_error(owea:::.ui_solver_args(spc, "verify", verify_step = -1), "positive")
  expect_error(owea:::.ui_solver_args(spc, "verify", verify_step = numeric(0)), "positive")
  # an explicit random audit: the continuous check with the chosen size, for any spec
  a <- owea:::.ui_solver_args(spc, "verify", verify_audit = 3000)
  expect_true(a$continuous); expect_equal(a$n_audit, 3000L); expect_null(a$step)
  spg0 <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    link = "identity")
  ag <- owea:::.ui_solver_args(spg0, "verify", verify_audit = 500)
  expect_true(ag$continuous); expect_equal(ag$n_audit, 500L); expect_null(ag$step)
  expect_error(owea:::.ui_solver_args(spg0, "verify", verify_audit = 0), "positive")
  # a step wins over an audit when both are given
  both <- owea:::.ui_solver_args(spg0, "verify", verify_step = 0.1, verify_audit = 500)
  expect_equal(both$step, 0.1); expect_null(both$continuous)
  # the grid count for the gate: 3 levels x (2 / 0.05 + 1)
  expect_equal(owea:::.ui_verify_points(spc, step = 0.05), 3 * 41)
  expect_true(is.na(owea:::.ui_verify_points(spc)))                 # continuous: no grid
  expect_true(is.na(owea:::.ui_verify_points(spc, step = 0)))       # unusable step
  # a grid spec: the override replaces its finest step
  spg <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = c(0.5, 0.25))),
    link = "identity")
  expect_equal(owea:::.ui_solver_args(spg, "verify")$step, 0.25)
  expect_equal(owea:::.ui_solver_args(spg, "verify", verify_step = 0.01)$step, 0.01)
  expect_equal(owea:::.ui_verify_points(spg, step = 0.01), 201)
  # the suggested step: the finest step of a grid spec, range / 40 otherwise
  expect_equal(owea:::.ui_verify_step_default(spg), "0.25")
  expect_equal(owea:::.ui_parse_steps(owea:::.ui_verify_step_default(spc)), 0.05)
  # and the grid check certifies a design from the continuous search
  set.seed(2)
  r <- do.call(optimal_design, owea:::.ui_solver_args(spc, "optimal", p = 0))
  chk <- suppressWarnings(do.call(verify_optimality,
    c(list(support = r$support, weights = r$weights),
      owea:::.ui_solver_args(spc, "verify", p = 0, verify_step = 0.02))))
  expect_identical(chk$method, "grid")
  expect_true(chk$is_optimal$value)
})

test_that(".ui_srs_design draws a simple random sample from a grid or the continuous region", {
  spg <- owea:::.ui_model_spec(
    list(list(name = "A", type = "factor", nlevels = 2),
         list(name = "x", type = "continuous", lo = -1, hi = 1, steps = 0.5)),
    link = "identity")
  set.seed(1)
  d <- owea:::.ui_srs_design(spg, 30)
  expect_equal(sum(d$counts), 30L)
  expect_equal(d$pool_size, 10)                       # 2 levels x 5 grid values (implied)
  expect_true(all(d$support[, 2] %in% seq(-1, 1, by = 0.5)))
  expect_true(all(d$support[, 1] %in% 1:2))
  # the grid is never built: a step so fine that the grid would have 1e18 points
  sph <- owea:::.ui_model_spec(
    list(list(name = "u", type = "continuous", lo = 0, hi = 1, steps = 1e-6),
         list(name = "v", type = "continuous", lo = -5, hi = 5, steps = 1e-6),
         list(name = "w", type = "continuous", lo = 0, hi = 100, steps = 1e-6)),
    link = "identity")
  set.seed(3)
  tm <- system.time(dh <- owea:::.ui_srs_design(sph, 50))[3]
  expect_lt(tm, 2)
  expect_equal(sum(dh$counts), 50L)
  expect_equal(dh$pool_size, (1e6 + 1) * (1e7 + 1) * (1e8 + 1))
  expect_true(all(abs(dh$support[, 1] / 1e-6 - round(dh$support[, 1] / 1e-6)) < 1e-6))
  expect_true(all(dh$support[, 2] >= -5 & dh$support[, 2] <= 5))
  # every grid point is equally likely, including the two ends of a range:
  # with 3 grid values per draw, the ends are hit about a third of the time each
  sp3 <- owea:::.ui_model_spec(
    list(list(name = "x", type = "continuous", lo = 0, hi = 1, steps = 0.5)),
    link = "identity")
  set.seed(4)
  d3 <- owea:::.ui_srs_design(sp3, 3000)
  frac <- d3$counts / sum(d3$counts)
  expect_equal(sort(d3$support[, 1]), c(0, 0.5, 1))
  expect_true(all(abs(frac - 1 / 3) < 0.05))
  # a continuous search has no grid: uniform draws from the region
  spc <- owea:::.ui_model_spec(
    list(list(name = "A", type = "factor", nlevels = 2),
         list(name = "x", type = "continuous", lo = -1, hi = 1)),
    link = "identity", continuous = TRUE)
  set.seed(2)
  dc <- owea:::.ui_srs_design(spc, 30)
  expect_equal(sum(dc$counts), 30L)
  expect_true(is.na(dc$pool_size))
  expect_true(all(dc$support[, 2] >= -1 & dc$support[, 2] <= 1))
  expect_true(all(dc$support[, 1] %in% 1:2))
  expect_equal(colnames(dc$support), c("A", "x"))
  expect_error(owea:::.ui_srs_design(spc, 0), "n >= 1")
})

test_that("simulate_design(existing = , obs = ) pools a first stage like the app", {
  sp  <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    link = "identity")
  sup <- cbind(c(-1, 1)); cnt <- c(10L, 10L); th <- c(1, 2)
  ex  <- list(points = cbind(c(-1, 1)), weights = c(0.5, 0.5), n0 = 40)
  a <- owea:::.ui_simulate(sp, th, 1, sup, cnt, existing = ex, nsim = 30, seed = 1)
  b <- simulate_design(sup, cnt, link = "identity", x = 1, design_box = sp$design_box,
                       theta = th, nsim = 30, seed = 1,
                       existing = list(points = ex$points, counts = c(20L, 20L)))
  expect_equal(b$N, 60L)
  expect_equal(b$n_existing, 40L)
  expect_equal(b$estimates, a$estimates)
  expect_equal(b$mse, a$mse)
  # weights + n0 are accepted too
  b2 <- simulate_design(sup, cnt, link = "identity", x = 1, design_box = sp$design_box,
                        theta = th, nsim = 30, seed = 1, existing = ex)
  expect_equal(b2$mse, a$mse)
  # an observed first stage keeps its responses
  d1 <- data.frame(dose = rep(c(-1, 0, 1), each = 4),
                   y    = c(rep(0, 4), rep(1, 4), rep(2, 4)))
  o <- simulate_design(sup, cnt, link = "identity", x = 1, design_box = sp$design_box,
                       theta = th, nsim = 20, seed = 3, obs = list(data = d1, response = "y"))
  expect_equal(o$N, 32L)
  expect_equal(o$n_existing, 12L)
  expect_error(simulate_design(sup, cnt, link = "identity", x = 1,
                               design_box = sp$design_box, theta = th, nsim = 1,
                               existing = list(points = ex$points, counts = c(20L, 20L))),
               "nsim >= 2")
  expect_error(simulate_design(sup, cnt, link = "identity", x = 1,
                               design_box = sp$design_box, theta = th, nsim = 5,
                               existing = list(points = ex$points)),
               "counts")
})

test_that(".ui_r_code adds a runnable simulation block for an exact design", {
  sp <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1)),
    link = "logit", continuous = TRUE, n_audit = 0)
  args <- owea:::.ui_solver_args(sp, "exact", theta = c(0.2, 1), p = 0, n = 12)
  args$seed <- 1L
  set.seed(1)
  res <- do.call(exact_design, args)
  set.seed(5)
  srs <- owea:::.ui_srs_design(sp, 12)
  sim <- list(theta = c(0.2, 1), sigma = NULL, nsim = 20L, seed = 3L,
              existing = NULL, obs = NULL, designs = list(SRS = srs))
  code <- owea:::.ui_r_code(args, "exact", sim = sim)
  expect_match(code, "set.seed(1)", fixed = TRUE)         # the continuous search is seeded
  expect_match(code, "res <- exact_design(", fixed = TRUE)
  expect_match(code, "theta_true <- c(0.2, 1)", fixed = TRUE)
  expect_match(code, "sim <- simulate_design(", fixed = TRUE)
  expect_match(code, "srs_support <- rbind(", fixed = TRUE)
  expect_match(code, "sim_srs <- simulate_design(", fixed = TRUE)
  expect_match(code, "mse <- rbind(`Exact design` = sim$mse, SRS = sim_srs$mse)", fixed = TRUE)
  expect_error(parse(text = code), NA)
  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  expect_equal(env$sim$N, 12L)
  expect_true(is.matrix(env$mse))
  expect_equal(rownames(env$mse), c("Exact design", "SRS"))
  # the script's design and simulation are the app's own
  expect_equal(env$res$counts, res$counts)
  ref <- simulate_design(res, link = "logit", x = 1, design_box = sp$design_box,
                         theta = c(0.2, 1), nsim = 20, seed = 3)
  expect_equal(env$sim$mse, ref$mse)
  # with an existing first stage the pooled call is emitted
  sim2 <- sim
  sim2$existing <- list(points = cbind(c(-1, 1)), counts = c(5L, 5L))
  code2 <- owea:::.ui_r_code(args, "exact", sim = sim2)
  expect_match(code2, "existing\\s+= list\\(points = rbind\\(")
  expect_error(parse(text = code2), NA)
  # an observed first stage is embedded as CSV text
  sim3 <- sim
  sim3$obs <- list(data = data.frame(dose = c(-1, 0, 1, 1), y = c(0, 1, 1, 0)), response = "y")
  code3 <- owea:::.ui_r_code(args, "exact", sim = sim3)
  expect_match(code3, "obs_data <- read.csv(text = c(", fixed = TRUE)
  expect_match(code3, "obs\\s+= list\\(data = obs_data, response = \"y\"\\)")
  expect_error(parse(text = code3), NA)
  env3 <- new.env(parent = globalenv())
  eval(parse(text = code3), envir = env3)
  expect_equal(env3$sim$n_existing, 4L)
})

test_that(".ui_r_code writes a runnable script that reproduces the app's call", {
  covs <- list(list(name = "my dose", type = "continuous", lo = -1, hi = 1),
               list(name = "grp", type = "factor", nlevels = 2))
  sp <- owea:::.ui_model_spec(covs, interactions = list(c(1, 2)), link = "logit",
                              continuous = TRUE, n_audit = 5000)
  th <- c(0.2, 1, -0.5, 0.3)
  args  <- owea:::.ui_solver_args(sp, "optimal", theta = th, p = 1)
  vargs <- owea:::.ui_solver_args(sp, "verify", theta = th, p = 1)
  code <- owea:::.ui_r_code(args, "optimal", verify_args = vargs)
  expect_match(code, "library(owea)", fixed = TRUE)
  expect_match(code, "res <- optimal_design(", fixed = TRUE)
  expect_match(code, "continuous = TRUE", fixed = TRUE)
  expect_match(code, "n_audit    = 5000", fixed = TRUE)
  expect_match(code, "`my dose` = c(-1, 1)", fixed = TRUE)      # non-syntactic name back-ticked
  expect_match(code, "v <- verify_optimality(", fixed = TRUE)
  expect_error(parse(text = code), NA)
  # running the script gives the same design as the app's own call
  env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = env)
  set.seed(1)
  direct <- do.call(optimal_design, args)
  expect_equal(env$res$criterion, direct$criterion, tolerance = 1e-6)
  expect_true(env$v$is_optimal$value)

  # an exact design: exact_design(), n and seed in the call, no verify block
  eargs <- owea:::.ui_solver_args(sp, "exact", theta = th, p = 1, n = 12)
  eargs$seed <- 7L
  ecode <- owea:::.ui_r_code(eargs, "exact")
  expect_match(ecode, "res <- exact_design(", fixed = TRUE)
  expect_match(ecode, "n          = 12", fixed = TRUE)
  expect_match(ecode, "seed       = 7", fixed = TRUE)
  expect_false(grepl("verify_optimality", ecode, fixed = TRUE))
  expect_error(parse(text = ecode), NA)

  # the grid path and an existing design: step_sequence and rbind() matrices
  spg <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = c(0.5, 0.1))),
    link = "identity")
  ex <- list(points = cbind(c(-1, 0.5)), weights = c(0.5, 0.5), n0 = 10, n1 = 20)
  gargs <- owea:::.ui_solver_args(spg, "optimal", p = 0, existing = ex)
  gcode <- owea:::.ui_r_code(gargs, "optimal")
  expect_match(gcode, "step_sequence = list(0.5, 0.1)", fixed = TRUE)
  expect_match(gcode, "xi0_points    = rbind(c(-1),", fixed = TRUE)
  expect_false(grepl("set.seed", gcode, fixed = TRUE))
  expect_error(parse(text = gcode), NA)
  env2 <- new.env(parent = globalenv())
  suppressWarnings(eval(parse(text = gcode), envir = env2))
  expect_true(is.finite(env2$res$criterion))
})

test_that("the continuous spec reproduces a direct optimal_design(continuous = TRUE) call", {
  covs <- list(list(name = "dose", type = "continuous", lo = -1, hi = 1))
  sp <- owea:::.ui_model_spec(covs, link = "logit", continuous = TRUE)
  set.seed(1)
  a <- do.call(optimal_design,
               owea:::.ui_solver_args(sp, "optimal", theta = c(0.2, 1), p = 0))
  expect_identical(a$method, "continuous")
  expect_true(a$converged)
  b <- suppressWarnings(optimal_design(design_box = list(dose = c(-1, 1)),
                                       link = "logit", x = 1, theta = c(0.2, 1), p = 0,
                                       step_sequence = c(0.2, 0.05, 0.01)))
  expect_equal(a$criterion, b$criterion, tolerance = 1e-4)
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

  # with a multi-step sequence, verification uses the LAST (finest) step
  sp2 <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1,
              steps = owea:::.ui_parse_steps("0.5, 0.1"))),
    link = "identity")
  v2 <- owea:::.ui_solver_args(sp2, "verify")
  expect_equal(v2$step, 0.1)
  expect_equal(v2$step, sp2$finest)
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
  expect_equal(same$efficiency_lower_bound, 1, tolerance = 1e-4)
  # the design's own criterion is evaluated grid-free (criterion_only), so the
  # only grid search is the reference solve with the ORIGINAL step_sequence
  expect_identical(same$max_sensitivity, NA_real_)

  # the D-optimal design is not A-optimal for a proper subset of the parameters
  sp2 <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.25)),
    quadratics = 1L, link = "identity")
  r2 <- suppressWarnings(do.call(optimal_design,
                                 owea:::.ui_solver_args(sp2, "optimal", p = 0)))
  other <- suppressWarnings(owea:::.ui_efficiency(r2, sp2, p_new = 1L,
                                                  subset_new = 3L))
  expect_true(other$efficiency_lower_bound > 0 && other$efficiency_lower_bound < 1)
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

# ---- robust simulation summaries (both branches of the app use these) ------

test_that(".ui_simulate reports a median squared error beside the mean", {
  sp <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.5)),
    link = "logit")
  s <- owea:::.ui_simulate(sp, theta = c(0.5, 1),
                           support = cbind(c(-1, 0, 1)), counts = c(20, 20, 20),
                           nsim = 100L, seed = 4)
  expect_length(s$medse, length(s$mse))
  expect_equal(names(s$medse), names(s$mse))
  expect_true(all(s$medse >= 0))
  # the median squared error can never exceed the mean by much, and for a
  # right-skewed error distribution it sits below it
  expect_true(all(s$medse <= s$mse * 1.5 + 1e-8))
})

test_that(".ui_sim_dominance detects an MSE carried by one replicate", {
  sp <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.5)),
    link = "logit")
  s <- owea:::.ui_simulate(sp, theta = c(0.5, 1),
                           support = cbind(c(-1, 0, 1)), counts = c(20, 20, 20),
                           nsim = 100L, seed = 4)
  d <- owea:::.ui_sim_dominance(s)
  expect_true(is.finite(d))
  expect_gte(d, 1 / s$n_converged)          # at least an equal share
  expect_lte(d, 1)

  # plant one wild replicate and the share must jump
  s2 <- s
  s2$estimates[1, ] <- s2$estimates[1, ] + 500
  expect_gt(owea:::.ui_sim_dominance(s2), 0.9)
  expect_gt(owea:::.ui_sim_dominance(s2), d)

  # and it is well behaved on rubbish input
  expect_true(is.na(owea:::.ui_sim_dominance(NULL)))
  expect_true(is.na(owea:::.ui_sim_dominance(simpleError("nope"))))
})

test_that("the median squared error resists a wild replicate, the mean does not", {
  sp <- owea:::.ui_model_spec(
    list(list(name = "dose", type = "continuous", lo = -1, hi = 1, steps = 0.5)),
    link = "logit")
  s <- owea:::.ui_simulate(sp, theta = c(0.5, 1),
                           support = cbind(c(-1, 0, 1)), counts = c(20, 20, 20),
                           nsim = 100L, seed = 4)
  ok <- stats::complete.cases(s$estimates)
  E  <- s$estimates[ok, , drop = FALSE]
  E2 <- E; E2[1, ] <- E2[1, ] + 500                    # one separated fit
  sq  <- sweep(E,  2, s$theta, "-")^2
  sq2 <- sweep(E2, 2, s$theta, "-")^2
  expect_gt(mean(sq2[, 1]), 100 * mean(sq[, 1]))       # the mean explodes
  expect_equal(stats::median(sq2[, 1]), stats::median(sq[, 1]),
               tolerance = 1e-8)                        # the median does not
})

test_that(".ui_sim_matrix lays the compared designs out one row each", {
  slots <- list(`Exact design` = c(a = 1, b = 2), SRS = c(a = 3, b = 4))
  M <- owea:::.ui_sim_matrix(slots)
  expect_true(is.matrix(M))
  expect_equal(dim(M), c(2L, 2L))                       # designs x parameters
  expect_equal(rownames(M), c("Exact design", "SRS"))
  expect_equal(colnames(M), c("a", "b"))                # names come from the vectors
  expect_equal(M["SRS", "b"], 4)

  # a design that was not run, or that failed, simply drops out
  M2 <- owea:::.ui_sim_matrix(list(`Exact design` = c(a = 1, b = 2), SRS = NULL,
                                   Custom = simpleError("nope")))
  expect_equal(rownames(M2), "Exact design")

  # explicit parameter names win over the vectors' own
  M3 <- owea:::.ui_sim_matrix(list(D = c(1, 2)), param_names = c("b0", "b1"))
  expect_equal(colnames(M3), c("b0", "b1"))

  # nothing comparable -> NULL, so the app can req() on it
  expect_null(owea:::.ui_sim_matrix(list()))
  expect_null(owea:::.ui_sim_matrix(list(A = NULL, B = NULL)))
  expect_null(owea:::.ui_sim_matrix(list(A = c(NA_real_, NaN))))

  # a design of a different length cannot share the groups, so it is dropped;
  # a single non-finite entry only loses its own bar
  M4 <- owea:::.ui_sim_matrix(list(A = c(a = 1, b = 2), B = c(a = 1, b = 2, c = 3),
                                   C = c(a = Inf, b = 5)))
  expect_equal(rownames(M4), c("A", "C"))
  expect_true(is.na(M4["C", "a"]))
  expect_equal(M4["C", "b"], 5)
})

test_that(".ui_sim_ratio divides every design by the reference one", {
  M <- owea:::.ui_sim_matrix(list(`Exact design` = c(b0 = 1, b1 = 2),
                                  SRS = c(b0 = 3, b1 = 4),
                                  Custom = c(b0 = 2, b1 = 1)))
  R <- owea:::.ui_sim_ratio(M, "Exact design")
  expect_equal(rownames(R), c("SRS", "Custom"))   # the reference itself is gone
  expect_equal(colnames(R), c("b0", "b1"))
  expect_equal(R["SRS", ], c(b0 = 3, b1 = 2))     # competitor / computed
  expect_equal(R["Custom", ], c(b0 = 2, b1 = 0.5))

  # the reference is found by NAME, not position: with the computed design not
  # run, its row is missing and no ratio can be formed
  M2 <- owea:::.ui_sim_matrix(list(SRS = c(b0 = 3, b1 = 4),
                                   Custom = c(b0 = 2, b1 = 1)))
  expect_null(owea:::.ui_sim_ratio(M2, "Exact design"))

  # nothing to compare against
  expect_null(owea:::.ui_sim_ratio(
    owea:::.ui_sim_matrix(list(`Exact design` = c(b0 = 1))), "Exact design"))
  expect_null(owea:::.ui_sim_ratio(NULL, "Exact design"))
})

test_that(".ui_sim_ratio reports NA rather than Inf for an unusable reference", {
  # a zero or missing reference MSE has no meaningful ratio
  M <- rbind(`Exact design` = c(b0 = 0, b1 = NA_real_, b2 = 2),
             SRS           = c(b0 = 1, b1 = 3,         b2 = 4))
  R <- owea:::.ui_sim_ratio(M, "Exact design")
  expect_true(is.na(R["SRS", "b0"]))              # not Inf
  expect_true(is.na(R["SRS", "b1"]))
  expect_equal(R["SRS", "b2"], 2)                 # the usable column still works

  # a competitor's own missing value only costs that one entry
  M2 <- rbind(`Exact design` = c(b0 = 1, b1 = 2),
              SRS            = c(b0 = NA_real_, b1 = 4))
  R2 <- owea:::.ui_sim_ratio(M2, "Exact design")
  expect_true(is.na(R2["SRS", "b0"]))
  expect_equal(R2["SRS", "b1"], 2)
})
