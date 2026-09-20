# optimal_design(continuous = TRUE): the grid-free search over continuous
# covariates (R/continuous.R), its analytic Jacobians, verification and the
# exact-design and app pass-throughs.

box3 <- list(c(-2, 2), c(-1, 1), c(-3, 3))
th3  <- c(1, -0.5, 0.5, 1)

test_that("continuous A-optimal logistic design matches the known optimum", {
  set.seed(1)
  r <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                      p = 1, continuous = TRUE)
  expect_true(r$converged)
  expect_lte(r$max_d, 1e-6)
  # the 7-point A-optimal design: tr(M^-1) = 19.828326, criterion = tr / v
  expect_equal(r$criterion, 19.828326 / 4, tolerance = 1e-6)
  expect_equal(nrow(r$support), 7L)
  expect_lt(abs(sum(r$weights) - 1), 1e-8)
  expect_identical(r$method, "continuous")
  expect_identical(r$jacobian, "analytic")
  expect_true(all(r$support[, 3] >= -3 & r$support[, 3] <= 3))
  # every support point has linear predictor +/- 1.3305
  eta <- drop(cbind(1, r$support) %*% th3)
  expect_true(all(abs(abs(eta) - 1.3305) < 2e-3))
  expect_gte(r$efficiency_lower_bound, 0.99999)
  expect_true(is.data.frame(r$history))
  expect_equal(r$grid_sizes, 27)                    # 3^3 coarse start grid
  expect_length(r$times, 2L)
  expect_true(r$global_check)
})

test_that("continuous D-optimal design is at least as good as the grid path", {
  set.seed(2)
  rc <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                       p = 0, continuous = TRUE)
  rg <- suppressWarnings(optimal_design(design_box = box3, link = "logit", x = 1:3,
                                        theta = th3, p = 0,
                                        step_sequence = c(0.5, 0.1)))
  expect_true(rc$converged)
  expect_lte(rc$criterion, rg$criterion + 1e-6)
})

test_that("user info_vector (finite differences) and info_matrix reproduce the built-in model", {
  set.seed(3)
  r0 <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                       p = 1, continuous = TRUE, n_audit = 0)
  rv <- optimal_design(design_box = box3, info_vector = logistic_vec, theta = th3,
                       p = 1, continuous = TRUE, n_audit = 0)
  rm <- optimal_design(design_box = box3, info_matrix = logistic_info,
                       p = 1, continuous = TRUE, n_audit = 0)
  expect_identical(rv$jacobian, "finite differences")
  expect_identical(rm$jacobian, "finite differences")
  expect_true(rv$converged && rm$converged)
  expect_equal(rv$criterion, r0$criterion, tolerance = 1e-6)
  expect_equal(rm$criterion, r0$criterion, tolerance = 1e-6)
  expect_true(is.na(rv$audit_max_d))              # n_audit = 0: no audit
})

test_that("a user info_jacobian is used, and a wrong shape is rejected", {
  set.seed(4)
  jac <- function(x, th) {
    q <- c(1, x); eta <- sum(q * th); mu <- plogis(eta)
    g <- exp(eta / 2) / (1 + exp(eta)); Jq <- rbind(0, diag(3))
    g * Jq + g * (0.5 - mu) * outer(q, as.numeric(crossprod(Jq, th)))
  }
  r <- optimal_design(design_box = box3, info_vector = logistic_vec, theta = th3,
                      p = 1, continuous = TRUE, info_jacobian = jac, n_audit = 0)
  expect_identical(r$jacobian, "user")
  expect_true(r$converged)
  expect_equal(r$criterion, 19.828326 / 4, tolerance = 1e-6)
  bad <- function(x, th) matrix(0, 2, 2)
  expect_error(optimal_design(design_box = box3, info_vector = logistic_vec,
                              theta = th3, p = 1, continuous = TRUE,
                              info_jacobian = bad), "info_jacobian")
})

test_that("analytic Jacobian and vectorised evaluation of formula-style models are exact", {
  f <- model_info_vector(list(g = 3, x = c(-1, 1), z = c(0, 2), h = 2),
                         link = "logit", f = 1:2, x = 1:2, fx = c(11, 22),
                         xx = c(11, 12), ff = 12)
  th <- seq(0.1, by = 0.1, length.out = length(attr(f, "coef_names")))
  x0 <- c(2, 0.3, 1.2, 1)
  Ja <- attr(f, "jacobian")(x0, th)
  Jn <- sapply(1:2, function(j) {
    h <- 1e-6; xp <- x0; xm <- x0
    xp[1 + j] <- xp[1 + j] + h; xm[1 + j] <- xm[1 + j] - h
    (f(xp, th) - f(xm, th)) / (2 * h)
  })
  expect_equal(Ja, Jn, tolerance = 1e-6)
  set.seed(5)
  X <- cbind(sample(1:3, 6, TRUE), runif(6, -1, 1), runif(6, 0, 2), sample(1:2, 6, TRUE))
  expect_equal(attr(f, "vectorized")(X, th),
               sapply(seq_len(6), function(i) f(X[i, ], th)))
  # loglinear and identity links
  fl <- model_info_vector(list(x = c(-1, 1)), link = "loglinear", x = 1, xx = 11)
  thl <- c(0.2, 0.5, -0.3)
  expect_equal(as.numeric(attr(fl, "jacobian")(0.4, thl)),
               (fl(0.4 + 1e-6, thl) - fl(0.4 - 1e-6, thl)) / 2e-6, tolerance = 1e-6)
  fi <- model_info_vector(list(x = c(-1, 1)), link = "identity", x = 1, xx = 11)
  expect_equal(as.numeric(attr(fi, "jacobian")(0.4, NULL)), c(0, 1, 0.8))
  # the attributes survive the internal normalisation of a one-argument model
  fn <- owea:::.normalize_info_vector(fi, NULL)
  expect_true(is.function(attr(fn, "jacobian")))
  expect_true(is.function(attr(fn, "vectorized")))
})

test_that("factor + continuous covariates: levels stay integers and the design certifies", {
  set.seed(6)
  r <- optimal_design(design_box = list(g = 3, x = c(-1, 1), z = c(0, 2)),
                      link = "logit", f = 1, x = 1:2, fx = 11,
                      theta = c(0.5, -0.3, 0.4, 1, -0.8, 0.2, 0.3), p = 0,
                      continuous = TRUE)
  expect_true(r$converged)
  expect_true(all(r$support[, 1] %in% 1:3))
  expect_true(all(r$support[, 2] >= -1 & r$support[, 2] <= 1))
  expect_true(all(r$support[, 3] >= 0 & r$support[, 3] <= 2))
  expect_equal(r$is_factor, c(TRUE, FALSE, FALSE))
})

test_that("subset, existing design and init_points work on the continuous path", {
  set.seed(7)
  rs <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                       p = 1, subset = 2:4, continuous = TRUE)
  expect_true(rs$converged)
  xi0 <- rbind(c(0, 0, 0), c(1, 1, 1))
  re <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                       p = 0, continuous = TRUE, xi0_points = xi0,
                       xi0_weights = c(0.5, 0.5), n0 = 10, n1 = 20)
  rg <- suppressWarnings(optimal_design(design_box = box3, link = "logit", x = 1:3,
                                        theta = th3, p = 0, step_sequence = c(0.5, 0.1),
                                        xi0_points = xi0, xi0_weights = c(0.5, 0.5),
                                        n0 = 10, n1 = 20))
  expect_true(re$converged)
  expect_lte(re$criterion, rg$criterion + 1e-6)
  # a user-supplied starting support replaces the coarse grid
  verts <- as.matrix(expand.grid(c(-2, 2), c(-1, 1), c(-3, 3)))
  ri <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                       p = 1, continuous = TRUE, init_points = verts)
  expect_true(ri$converged)
  expect_equal(ri$criterion, 19.828326 / 4, tolerance = 1e-6)
  expect_equal(ri$grid_sizes, 8)
  expect_error(optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                              p = 1, continuous = TRUE,
                              init_points = rbind(c(0, 0, 9), verts)),
               "inside the design box")
})

test_that("verify_optimality(continuous = TRUE) confirms the design and flags a perturbed one", {
  set.seed(8)
  r <- optimal_design(design_box = box3, link = "logit", x = 1:3, theta = th3,
                      p = 1, continuous = TRUE)
  v <- verify_optimality(r$support, r$weights, design_box = box3, link = "logit",
                         x = 1:3, theta = th3, p = 1, continuous = TRUE)
  expect_true(v$is_optimal$value)
  expect_identical(v$method, "continuous")
  expect_match(v$is_optimal$note, "audit")
  expect_equal(v$criterion, r$criterion, tolerance = 1e-8)
  expect_length(v$maximiser, 3L)
  v2 <- verify_optimality(r$support, rev(r$weights), design_box = box3,
                          link = "logit", x = 1:3, theta = th3, p = 1,
                          continuous = TRUE, n_audit = 0)
  expect_false(v2$is_optimal$value)
  expect_gt(v2$max_sensitivity, 1e-2)
  expect_warning(
    verify_optimality(rbind(c(0, 0, 5)), 1, design_box = box3, link = "logit",
                      x = 1:3, theta = th3, p = 1, continuous = TRUE, n_audit = 0),
    "outside the design box")
  # a criterion-only call needs no search
  co <- verify_optimality(r$support, r$weights, design_box = box3, link = "logit",
                          x = 1:3, theta = th3, p = 1, continuous = TRUE,
                          criterion_only = TRUE)
  expect_identical(co$max_sensitivity, NA_real_)
  expect_equal(co$criterion, r$criterion, tolerance = 1e-8)
})

test_that("exact_design(continuous = TRUE) rounds the continuous design", {
  set.seed(9)
  ex <- exact_design(30, design_box = box3, link = "logit", x = 1:3, theta = th3,
                     p = 1, continuous = TRUE, seed = 1)
  expect_equal(sum(ex$counts), 30L)
  expect_true(is.finite(ex$criterion))
  expect_gt(ex$efficiency_lower_bound, 0.95)
  expect_identical(ex$approx$method, "continuous")
})

test_that("continuous = TRUE on an all-factor box falls back to the enumeration", {
  expect_message(
    r <- suppressWarnings(optimal_design(design_box = list(2, 3), link = "logit",
                                         f = 1:2, theta = c(0.2, 0.5, -0.3, 0.4),
                                         p = 0, continuous = TRUE)),
    "every covariate is a factor")
  expect_true(r$converged)
  expect_null(r$method)
})

test_that("check_global with continuous = TRUE needs global_step, then scans a grid", {
  set.seed(10)
  expect_warning(
    r <- optimal_design(design_box = list(c(-1, 1)), link = "logit", x = 1,
                        theta = c(0.2, 1), p = 0, continuous = TRUE,
                        check_global = TRUE),
    "global_step")
  expect_true(r$converged)
  r2 <- optimal_design(design_box = list(c(-1, 1)), link = "logit", x = 1,
                       theta = c(0.2, 1), p = 0, continuous = TRUE,
                       check_global = TRUE, global_step = 0.01)
  expect_true(r2$global_check)
  expect_true(is.finite(r2$grid_check_max_d))
})

test_that("the grid path is untouched by the new arguments' defaults", {
  r <- suppressWarnings(optimal_design(design_box = box3, link = "logit", x = 1:3,
                                       theta = th3, p = 1,
                                       step_sequence = c(0.5, 0.1)))
  expect_null(r$method)
  expect_null(r$history)
  v <- verify_optimality(r$support, r$weights, design_box = box3, step = 0.1,
                         link = "logit", x = 1:3, theta = th3, p = 1)
  expect_identical(v$method, "grid")
})

test_that("six continuous covariates solve in seconds without a grid", {
  set.seed(11)
  box6 <- replicate(6, c(-2, 2), simplify = FALSE)
  th6  <- c(1, -0.5, 0.5, 1, -1, 0.7, -0.7)
  tm <- system.time(
    r <- optimal_design(design_box = box6, link = "logit", x = 1:6, theta = th6,
                        p = 0, continuous = TRUE))[3]
  expect_true(r$converged)
  expect_lt(tm, 60)
  expect_equal(r$grid_sizes, 3^6)
})
