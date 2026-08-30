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
  crit <- vapply(c("minmax", "minmaxmedian", "random", "MA"), function(m) {
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

test_that("init_method = 'MA' warm-starts a higher-dimensional quadratic model", {
  # the second-order linear model from the user report (identity link).
  args <- list(x = c(1, 2, 3, 4, 5),
               xx = c(11,12,13,14,15,22,23,24,25,33,34,35,44,45,55),
               link = "identity",
               design_box = replicate(5, c(-1, 1), simplify = FALSE),
               step_sequence = c(1), p = 0, check_global = TRUE)
  res_ma <- suppressWarnings(do.call(optimal_design,
                                     c(args, list(init_method = "MA"))))
  res_mm <- suppressWarnings(do.call(optimal_design,
                                     c(args, list(init_method = "minmax"))))
  expect_true(res_ma$converged)
  expect_lt(abs(res_ma$criterion - res_mm$criterion), 1e-6)
  # 'MA' accepts a lower-case alias too.
  res_lc <- suppressWarnings(do.call(optimal_design,
                                     c(args, list(init_method = "ma"))))
  expect_lt(abs(res_lc$criterion - res_mm$criterion), 1e-6)
  # ma_max_iter is accepted and reaches the same optimum with a tiny cap.
  res_it <- suppressWarnings(do.call(optimal_design,
                              c(args, list(init_method = "MA", ma_max_iter = 5))))
  expect_true(res_it$converged)
  expect_lt(abs(res_it$criterion - res_mm$criterion), 1e-6)
})

test_that("solver = 'MA' solves D-optimality directly and matches OWEA", {
  args <- list(x = c(1, 2, 3, 4, 5),
               xx = c(11,12,13,14,15,22,23,24,25,33,34,35,44,45,55),
               link = "identity",
               design_box = replicate(5, c(-1, 1), simplify = FALSE),
               step_sequence = c(1), p = 0)
  res_ma <- suppressWarnings(do.call(optimal_design,
                                     c(args, list(solver = "MA"))))
  res_ow <- suppressWarnings(do.call(optimal_design,
                                     c(args, list(solver = "owea"))))
  expect_true(res_ma$converged)
  expect_lte(res_ma$max_d, 1e-6)
  expect_lt(abs(sum(res_ma$weights) - 1), 1e-8)
  expect_lt(abs(res_ma$criterion - res_ow$criterion), 1e-6)
  # "multiplicative" is an accepted alias.
  res_al <- suppressWarnings(do.call(optimal_design,
                              c(args, list(solver = "multiplicative"))))
  expect_lt(abs(res_al$criterion - res_ow$criterion), 1e-6)
  # unknown solver errors.
  expect_error(do.call(optimal_design, c(args, list(solver = "bogus"))))
})

test_that("solver = 'MA' certifies A-optimality, subset, grad_g (>= OWEA)", {
  # Fixed candidate grid keeps it a single, deterministic solve.
  Xg <- make_grid(c(-2, -1, -3), c(2, 1, 3), 0.5)
  base <- list(info_matrix = logistic_info, candidate_set = Xg)
  gg <- function(th) matrix(c(0,1,0,0, 0,0,0,1), nrow = 2, byrow = TRUE)
  cases <- list(
    "A-opt full" = list(p = 1),
    "A subset"   = list(p = 1, subset = c(2, 4)),
    "D grad_g"   = list(p = 0, grad_g = gg))
  for (nm in names(cases)) {
    a <- c(base, cases[[nm]])
    r_ma <- suppressWarnings(do.call(optimal_design, c(a, list(solver = "MA"))))
    r_ow <- suppressWarnings(do.call(optimal_design, c(a, list(solver = "owea"))))
    expect_true(r_ma$converged, info = nm)
    expect_lte(r_ma$max_d, 1e-6)
    # a certified MA design is optimal, so no worse than OWEA (which can stall).
    expect_lte(r_ma$criterion, r_ow$criterion + 1e-5)
  }
})

test_that("solver = 'MA' matches OWEA with an existing design (n0 > 0)", {
  args <- list(x = c(1, 2, 3), xx = c(11, 22, 33), link = "identity",
               design_box = list(c(-1, 1), c(-1, 1), c(-1, 1)),
               step_sequence = c(1), p = 0,
               xi0_points = matrix(c(0,0,0, 1,1,1), nrow = 2, byrow = TRUE),
               xi0_weights = c(0.5, 0.5), n0 = 10, n1 = 20)
  r_ma <- suppressWarnings(do.call(optimal_design, c(args, list(solver = "MA"))))
  r_ow <- suppressWarnings(do.call(optimal_design, c(args, list(solver = "owea"))))
  expect_true(r_ma$converged)
  expect_lt(abs(r_ma$criterion - r_ow$criterion), 1e-5)
})

test_that("solver = 'MA' handles c-optimality (rank-1 wb) gracefully", {
  # a single linear combination of parameters -> rank-1 wb.  c-optimal designs
  # are typically singular, which is outside MA's guaranteed-convergent class;
  # the solver must still run without error and return a finite design (falling
  # back to OWEA per grid when it cannot certify).
  args <- list(info_matrix = logistic_info,
               design_box = list(c(-2, 2), c(-1, 1), c(-3, 3)),
               step_sequence = c(0.2, 0.1), p = 0,
               grad_g = function(th) matrix(c(0, 1, -1, 0), nrow = 1))
  r_ma <- suppressWarnings(do.call(optimal_design, c(args, list(solver = "MA"))))
  expect_true(is.finite(r_ma$criterion))
  expect_true(nrow(as.matrix(r_ma$support)) >= 1L)
})

test_that("total_time is the whole call, times is the solver only", {
  # A grid big enough that building it costs more than solving on it: the two
  # numbers must then differ, and total_time must be the larger.
  f <- function(x, th) { q <- c(1, x[1], x[2], x[3]); e <- sum(q * th)
                         (exp(e / 2) / (1 + exp(e))) * q }
  th  <- c(0.2, 0.8, -0.6, 0.5)
  box <- list(x1 = c(-1, 1), x2 = c(-1, 1), x3 = c(-1, 1))

  wall <- system.time(
    r <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = 0,
                                         design_box = box,
                                         step_sequence = 0.02)))[3]
  expect_true(is.finite(r$total_time))
  # it tracks the call's wall clock, not the solver alone
  expect_lte(r$total_time, wall + 0.05)
  expect_gte(r$total_time, sum(r$times))
  # and on this grid the untimed work (grid build, model evaluation) is real
  expect_gt(r$total_time, sum(r$times))

  # the per-stage vector is untouched: one entry per stage of the sequence
  m <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = 0,
                                       design_box = box,
                                       step_sequence = c(0.5, 0.2, 0.05)))
  expect_length(m$times, 3L)
  expect_length(m$grid_sizes, 3L)
  expect_gte(m$total_time, sum(m$times))

  # check_global's verification is inside the call, so it counts too
  g <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = 0,
                                       design_box = box,
                                       step_sequence = c(0.5, 0.2, 0.05),
                                       check_global = TRUE))
  expect_gte(g$total_time, sum(g$times))

  # the candidate_set path reports a whole-call time as well
  Xg <- candidate_grid(box, 0.1)
  s <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = 0,
                                       candidate_set = Xg))
  expect_true(is.finite(s$total_time))
  expect_gte(s$total_time, sum(s$times))
})
