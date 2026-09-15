# Tests for compound_design() / compound_sensitivity() -- the compound
# (multi-criterion, multi-model) code path.  Nothing here touches the
# single-criterion functions except to check that the two agree.

tr <- function(X) sum(diag(X))

# three genuinely different models on a shared design region
th1 <- c(0.5, 1, -1); th2 <- c(0.5, 1, -1, 0.5); th3 <- c(0.2, 0.3, -0.3)
f1 <- function(x, th) { q <- c(1, x[1], x[2]); e <- sum(q * th)
                        (exp(e / 2) / (1 + exp(e))) * q }              # logistic, k=3
f2 <- function(x, th) { q <- c(1, x[1], x[2], x[1] * x[2]); e <- sum(q * th)
                        (exp(e / 2) / (1 + exp(e))) * q }              # + interaction, k=4
f3 <- function(x, th) { q <- c(1, x[1], x[2]); e <- sum(q * th)
                        exp(e / 2) * q }                               # Poisson log, k=3

cmp3 <- list(list(info_vector = f1, theta = th1, p = 0, name = "m1"),
             list(info_vector = f2, theta = th2, p = 0, name = "m2"),
             list(info_vector = f3, theta = th3, p = 1, subset = c(2, 3), name = "m3"))
Xg <- as.matrix(expand.grid(x1 = seq(-2, 2, by = 0.25),
                            x2 = seq(-2, 2, by = 0.25)))

test_that("a single D component reproduces optimal_design()", {
  ro <- suppressWarnings(optimal_design(info_vector = f1, theta = th1,
                                        candidate_set = Xg, p = 0))
  rc <- compound_design(list(list(info_vector = f1, theta = th1, p = 0)),
                        candidate_set = Xg, efficiency = FALSE)
  # optimal_design reports log|S|/v, so Psi_D = exp(-criterion)
  expect_equal(rc$criterion, exp(-ro$criterion), tolerance = 1e-8)
  expect_true(rc$converged)
})

test_that("a single A component reproduces optimal_design()", {
  ro <- suppressWarnings(optimal_design(info_vector = f3, theta = th3,
                                        candidate_set = Xg, p = 1,
                                        subset = c(2, 3)))
  rc <- compound_design(list(list(info_vector = f3, theta = th3, p = 1,
                                  subset = c(2, 3))),
                        candidate_set = Xg, efficiency = FALSE)
  # optimal_design reports tr(S)/v, so Psi_A = 1/criterion
  expect_equal(rc$criterion, 1 / ro$criterion, tolerance = 1e-8)
  expect_true(rc$converged)
})

test_that("the compound equivalence theorem holds at the optimum", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  expect_true(r$converged)
  expect_lt(r$max_d, 1e-8)

  v <- compound_sensitivity(r, cmp3, Xg)
  # sensitivity <= 0 everywhere ...
  expect_lt(max(v$sensitivity), 1e-8)
  # ... and = 0 on the support
  idx <- apply(r$support, 1, function(p)
    which.min(rowSums((Xg - matrix(p, nrow(Xg), 2, byrow = TRUE))^2)))
  expect_lt(max(abs(v$sensitivity[idx])), 1e-8)
  # Euler identity  sum_j tr(K_j M_j) = Psi_alpha
  expect_lt(abs(v$euler), 1e-10)
  # weighted sensitivities sum to zero
  expect_lt(abs(sum(r$weights * v$sensitivity[idx])), 1e-10)
})

test_that("components may have different parameter dimensions", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  ks <- vapply(r$information, nrow, integer(1))
  expect_equal(unname(ks), c(3L, 4L, 3L))
  expect_equal(length(r$psi), 3L)
})

test_that("efficiency weighting bounds Psi_alpha in [max(alpha), 1]", {
  al <- c(0.4, 0.35, 0.25)
  r  <- compound_design(cmp3, alpha = al, candidate_set = Xg)
  expect_true(all(r$efficiency_lower_bound > 0), info = "efficiencies positive")
  expect_true(all(r$efficiency_lower_bound <= 1 + 1e-10), info = "efficiencies at most 1")
  expect_gte(r$criterion, max(al) - 1e-10)
  expect_lte(r$criterion, 1 + 1e-10)
})

test_that("the compound optimum beats every single-component optimum", {
  al <- c(0.4, 0.35, 0.25)
  r  <- compound_design(cmp3, alpha = al, candidate_set = Xg)
  for (j in seq_len(3)) {
    rj <- compound_design(cmp3[j], candidate_set = Xg,
                          psi_star = r$psi_star[j], efficiency = TRUE)
    # score the single-component design under the FULL compound criterion
    psi_j <- vapply(seq_len(3), function(l) {
      cl <- cmp3[[l]]
      rr <- compound_design(cmp3[l], candidate_set = Xg,
                            psi_star = r$psi_star[l], efficiency = TRUE)
      rr$psi[[1]]
    }, numeric(1))
    expect_true(is.finite(rj$criterion))
  }
  # the compound criterion at each pure optimum, computed directly
  score <- function(des) sum(r$at * vapply(seq_len(3), function(l) {
    cl <- cmp3[[l]]
    M <- Reduce(`+`, lapply(seq_len(nrow(des$support)), function(i)
      des$weights[i] * tcrossprod(cl$info_vector(des$support[i, ], cl$theta))))
    G <- if (is.null(cl$subset)) diag(nrow(M)) else diag(nrow(M))[cl$subset, , drop = FALSE]
    S <- G %*% solve(M) %*% t(G)
    if (cl$p == 0) det(solve(S))^(1 / nrow(S)) else nrow(S) / tr(S)
  }, numeric(1)))
  for (j in seq_len(3)) {
    rj <- compound_design(cmp3[j], candidate_set = Xg, efficiency = FALSE)
    expect_lte(score(rj), r$criterion + 1e-8)
  }
})

test_that("the gradient and Hessian match finite differences", {
  # a fixed support, so the derivatives are unambiguous
  idx <- c(1L, 17L, 145L, 289L, 80L)
  w   <- c(0.25, 0.2, 0.2, 0.2, 0.15)
  dat <- lapply(cmp3, function(cc)
    sapply(seq_len(nrow(Xg)), function(i) cc$info_vector(Xg[i, ], cc$theta)))
  G   <- list(diag(3), diag(4), diag(3)[c(2, 3), , drop = FALSE])
  I0  <- lapply(dat, function(d) matrix(0, nrow(d), nrow(d)))
  pp  <- c(0L, 0L, 1L); at <- c(1.4, 1.2, 0.05)

  crit <- function(ww)
    owea:::compound_criterion_cpp(dat, rep(0L, 3), G, I0, pp, at, idx, ww)
  Pw <- function(u) crit(c(u, 1 - sum(u)))
  u0 <- w[1:4]

  gh <- owea:::compound_grad_hess_cpp(dat, rep(0L, 3), G, I0, pp, at, idx, w)
  ng <- vapply(1:4, function(i) { e <- rep(0, 4); e[i] <- 1e-6
                                  (Pw(u0 + e) - Pw(u0 - e)) / 2e-6 }, numeric(1))
  expect_equal(as.numeric(gh$gradient), ng, tolerance = 1e-6)

  nh <- outer(1:4, 1:4, Vectorize(function(i, l) {
    h <- 1e-4; ei <- rep(0, 4); ei[i] <- h; el <- rep(0, 4); el[l] <- h
    (Pw(u0+ei+el) - Pw(u0+ei-el) - Pw(u0-ei+el) + Pw(u0-ei-el)) / (4 * h^2) }))
  expect_equal(gh$hessian, nh, tolerance = 1e-4)

  # structural properties: symmetric and negative semi-definite (concavity)
  expect_equal(gh$hessian, t(gh$hessian), tolerance = 1e-12)
  expect_lte(max(eigen(gh$hessian, symmetric = TRUE, only.values = TRUE)$values), 1e-8)
})

test_that("raw and efficiency weighting give different designs", {
  al <- c(0.4, 0.35, 0.25)
  re <- compound_design(cmp3, alpha = al, candidate_set = Xg, efficiency = TRUE)
  rr <- compound_design(cmp3, alpha = al, candidate_set = Xg, efficiency = FALSE)
  # same components, different effective weights -> different efficiencies
  eff_raw <- rr$psi / re$psi_star
  expect_false(isTRUE(all.equal(as.numeric(re$efficiency_lower_bound), as.numeric(eff_raw),
                                tolerance = 1e-4)))
  expect_true(all(is.na(rr$psi_star)))
})

test_that("the design_box path agrees with the candidate_set path", {
  rb <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25),
                        design_box = list(c(-2, 2), c(-2, 2)),
                        step_sequence = c(0.25))
  rs <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  expect_equal(rb$criterion, rs$criterion, tolerance = 1e-6)
})

test_that("an existing design is combined per model", {
  x0 <- rbind(c(-1, 1), c(1, -1), c(0, 0), c(2, 2))
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg,
                       xi0_points = x0, xi0_weights = rep(0.25, 4),
                       n0 = 30, n1 = 90)
  expect_true(r$converged)
  # efficiencies stay in (0, 1]: the reference solves used the same stage
  # structure, so the existing design cannot push them above 1
  expect_true(all(r$efficiency_lower_bound > 0 & r$efficiency_lower_bound <= 1 + 1e-10))
})

test_that("weights are normalised and recycled", {
  r1 <- compound_design(cmp3, alpha = c(2, 2, 2), candidate_set = Xg)
  r2 <- compound_design(cmp3, alpha = NULL,      candidate_set = Xg)
  expect_equal(unname(r1$alpha), rep(1 / 3, 3))
  expect_equal(r1$criterion, r2$criterion, tolerance = 1e-10)
})

test_that("invalid input is rejected", {
  expect_error(compound_design(list(), candidate_set = Xg), "non-empty list")
  expect_error(compound_design(cmp3, alpha = c(1, 2), candidate_set = Xg),
               "one weight per component")
  expect_error(compound_design(cmp3, alpha = c(-1, 1, 1), candidate_set = Xg),
               "nonnegative")
  expect_error(compound_design(cmp3, candidate_set = Xg,
                               psi_star = c(1, 2)), "one value per component")
  expect_error(compound_design(list(list(info_vector = f1, theta = th1, p = 3)),
                               candidate_set = Xg), "D-optimality")
  expect_error(compound_design(cmp3), "candidate_set")
  expect_error(compound_design(list(list(theta = th1)), candidate_set = Xg),
               "info_vector")
})

test_that("print.compound_design runs", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  expect_output(print(r), "Compound optimal design")
  expect_output(print(r), "Psi_alpha")
})

test_that("compound_criterion scores an arbitrary design", {
  fac <- rbind(c(-2, -2), c(-2, 2), c(2, -2), c(2, 2), c(0, 0))
  a <- compound_criterion(fac, rep(0.2, 5), cmp3, alpha = c(0.4, 0.35, 0.25),
                          candidate_set = Xg)
  expect_true(is.finite(a$criterion))
  expect_length(a$psi, 3L)
  expect_true(all(a$efficiency_lower_bound > 0 & a$efficiency_lower_bound <= 1 + 1e-10))
  expect_false(a$is_optimal)            # a factorial is not compound-optimal
  expect_gt(a$max_d, 0)
  expect_equal(unname(vapply(a$information, nrow, integer(1))), c(3L, 4L, 3L))
})

test_that("compound_criterion needs no design space when psi_star is given", {
  fac <- rbind(c(-2, -2), c(-2, 2), c(2, -2), c(2, 2), c(0, 0))
  ref <- compound_criterion(fac, rep(0.2, 5), cmp3, candidate_set = Xg)
  a   <- compound_criterion(fac, rep(0.2, 5), cmp3, psi_star = ref$psi_star)
  expect_equal(a$criterion, ref$criterion)
  expect_null(a$sensitivity)
  expect_error(compound_criterion(fac, rep(0.2, 5), cmp3), "psi_star")
})

test_that("compound_criterion certifies a compound_design() result", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  d <- compound_criterion(r, components = cmp3, candidate_set = Xg)
  expect_equal(d$criterion, r$criterion, tolerance = 1e-12)
  expect_true(d$is_optimal)
  expect_lt(abs(d$euler), 1e-10)
  expect_equal(unname(d$alpha), unname(r$alpha))   # weights reused
})

test_that("an off-grid support is scored exactly, not snapped", {
  ref <- compound_criterion(rbind(c(-2, -2), c(2, 2), c(-2, 2), c(2, -2)),
                            rep(0.25, 4), cmp3, candidate_set = Xg)
  off <- rbind(c(-1.937, -1.914), c(-1.873, 1.902),
               c(1.888, -1.955), c(1.961, 1.877))   # none on the 0.25 grid
  e1 <- compound_criterion(off, rep(0.25, 4), cmp3, psi_star = ref$psi_star)
  e2 <- compound_criterion(off, rep(0.25, 4), cmp3, psi_star = ref$psi_star,
                           candidate_set = Xg)
  expect_equal(e1$criterion, e2$criterion, tolerance = 1e-14)
})

test_that("compound_sensitivity accepts a bare support matrix", {
  ref <- compound_criterion(rbind(c(-2, -2), c(2, 2), c(-2, 2), c(2, -2)),
                            rep(0.25, 4), cmp3, candidate_set = Xg)
  s <- compound_sensitivity(rbind(c(-2, -2), c(2, 2), c(-2, 2), c(2, -2)),
                            cmp3, Xg, weights = rep(0.25, 4),
                            psi_star = ref$psi_star)
  expect_length(s$sensitivity, nrow(Xg))
  expect_true(is.finite(s$max_d))
})

test_that("compound_criterion validates its input", {
  fac <- rbind(c(-2, -2), c(2, 2))
  expect_error(compound_criterion(fac, components = cmp3), "'weights' is required")
  expect_error(compound_criterion(fac, c(-0.5, 1.5), cmp3, candidate_set = Xg),
               "nonnegative")
  expect_warning(compound_criterion(fac, c(0.3, 0.3), cmp3, candidate_set = Xg),
                 "normalising")
})

test_that("cross_efficiency is returned and has the right structure", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  ce <- r$cross_efficiency
  expect_false(is.null(ce))
  expect_equal(dim(ce), c(4L, 3L))                       # J+1 rows, J columns
  expect_equal(colnames(ce), names(r$psi))
  expect_equal(rownames(ce),
               c(paste("optimal for", names(r$psi)), "THIS DESIGN"))
  # every entry is an efficiency
  expect_true(all(ce > 0 & ce <= 1 + 1e-10))
  # the diagonal of the reference block is 1: each design is optimal for itself
  expect_equal(diag(ce[1:3, , drop = FALSE]), rep(1, 3),
               tolerance = 1e-6, ignore_attr = TRUE)
  # the last row is this design's efficiencies
  expect_equal(as.numeric(ce[4, ]), as.numeric(r$efficiency_lower_bound), tolerance = 1e-10)
  # reference designs are kept too
  expect_length(r$reference_designs, 3L)
  expect_true(all(vapply(r$reference_designs,
                         function(d) !is.null(d$support), logical(1))))
})

test_that("cross_efficiency rows match compound_criterion on those designs", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  for (j in seq_len(3)) {
    d <- r$reference_designs[[j]]
    # both sides are lower bounds: the ratio times the reference designs' bounds
    e <- compound_criterion(d$support, d$weights, cmp3, psi_star = r$psi_star,
                            reference_bound = r$reference_bound)$efficiency_lower_bound
    expect_equal(as.numeric(r$cross_efficiency[j, ]), as.numeric(e),
                 tolerance = 1e-10)
  }
})

test_that("cross_efficiency is absent when no reference solves were run", {
  r1 <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg,
                        efficiency = FALSE)
  expect_null(r1$cross_efficiency)
  ps <- compound_design(cmp3, candidate_set = Xg)$psi_star
  r2 <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg,
                        psi_star = ps)
  expect_null(r2$cross_efficiency)
  expect_output(print(r1), "Compound optimal design")   # still prints
})

test_that("print shows the cross-efficiency table", {
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)
  expect_output(print(r), "efficiency of each design under every component")
  expect_output(print(r), "THIS DESIGN")
})

# ---- exact (integer-run) compound designs ---------------------------------

test_that("compound_exact_design allocates exactly n runs", {
  r <- compound_exact_design(n = 30, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                             candidate_set = Xg, seed = 1)
  expect_s3_class(r, "compound_exact_design")
  expect_equal(sum(r$counts), 30L)
  expect_true(all(r$counts > 0))
  expect_true(all(r$counts == round(r$counts)))
  expect_equal(nrow(r$support), length(r$counts))
  expect_equal(r$weights, r$counts / 30)
  expect_equal(sum(r$weights), 1)
})

test_that("the exact design is bounded by the approximate optimum", {
  r <- compound_exact_design(n = 30, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                             candidate_set = Xg, seed = 1)
  # Psi_alpha is MAXIMISED, so the approximate optimum is an upper bound
  expect_lte(r$criterion, r$criterion_approx + 1e-10)
  expect_gt(r$efficiency_exact_lower_bound, 0)
  expect_lte(r$efficiency_exact_lower_bound, 1 + 1e-10)
  expect_true(all(r$efficiency_lower_bound <= 1 + 1e-10))
})

test_that("the exact criterion agrees with compound_criterion()", {
  r <- compound_exact_design(n = 40, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                             candidate_set = Xg, seed = 2)
  v <- compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                          psi_star = r$psi_star)
  expect_equal(v$criterion, r$criterion, tolerance = 1e-10)
  expect_equal(as.numeric(v$psi), as.numeric(r$psi), tolerance = 1e-10)
})

test_that("a larger sample size gets closer to the approximate optimum", {
  e1 <- compound_exact_design(n = 10, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                              candidate_set = Xg, seed = 1)$efficiency_exact_lower_bound
  e2 <- compound_exact_design(n = 200, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                              candidate_set = Xg, seed = 1)$efficiency_exact_lower_bound
  expect_gt(e2, e1)
  expect_gt(e2, 0.99)
})

test_that("compound_exact_design carries the cross-efficiency table", {
  r <- compound_exact_design(n = 30, components = cmp3, alpha = c(0.4, 0.35, 0.25),
                             candidate_set = Xg, seed = 1)
  ce <- r$cross_efficiency
  expect_equal(dim(ce), c(4L, 3L))
  expect_equal(rownames(ce)[4], "THIS DESIGN")
  # the last row is THIS exact design, not the approximate one
  expect_equal(as.numeric(ce[4, ]), as.numeric(r$efficiency_lower_bound), tolerance = 1e-10)
  expect_output(print(r), "Exact compound optimal design")
  expect_output(print(r), "efficiency of each design under every component")
})

test_that("compound_exact_design rejects an impossible sample size", {
  expect_error(compound_exact_design(n = 2, components = cmp3,
                                     candidate_set = Xg),
               "below the")
  expect_error(compound_exact_design(n = 0, components = cmp3,
                                     candidate_set = Xg),
               "positive integer")
})

test_that("a step sequence does not materialise the finest grid over the box", {
  # Two continuous covariates on [-2, 2] at a final step of 0.05 is an 81 x 81
  # = 6561-point grid; the coarse-plus-neighbourhood construction must stay far
  # below that while still finding a design of comparable quality.
  box <- list(x1 = c(-2, 2), x2 = c(-2, 2))
  r <- compound_exact_design(n = 30, components = cmp3,
                             alpha = c(0.4, 0.35, 0.25),
                             design_box = box, step_sequence = c(1, 0.25, 0.05),
                             seed = 1)
  expect_equal(sum(r$counts), 30L)

  # the same problem searched over the whole fine grid, for comparison
  full <- compound_exact_design(n = 30, components = cmp3,
                                alpha = c(0.4, 0.35, 0.25),
                                candidate_set = as.matrix(expand.grid(
                                  x1 = seq(-2, 2, by = 0.05),
                                  x2 = seq(-2, 2, by = 0.05))),
                                seed = 1)
  expect_equal(full$n_candidates, 6561L)
  expect_lt(r$n_candidates, full$n_candidates / 3)

  # every support point still lies in the box, and on the fine lattice
  sup <- r$support
  expect_true(all(sup >= -2 - 1e-9 & sup <= 2 + 1e-9))
  expect_true(all(abs(sup / 0.05 - round(sup / 0.05)) < 1e-6))

  # and the design is still good: close to what the full fine grid finds, and
  # bounded by its own approximate optimum
  expect_gt(r$criterion, 0.97 * full$criterion)
  expect_lte(r$criterion, r$criterion_approx + 1e-10)
})

test_that("a single-step sequence still searches the whole fine grid", {
  # coarse == fine, so the candidate set is the full grid exactly as before
  box <- list(x1 = c(-2, 2), x2 = c(-2, 2))
  r <- compound_exact_design(n = 20, components = cmp3,
                             alpha = c(0.4, 0.35, 0.25),
                             design_box = box, step_sequence = 0.25, seed = 3)
  expect_equal(r$n_candidates, 17L * 17L)
  expect_equal(sum(r$counts), 20L)
})

test_that("compound_criterion caps an enormous scanning grid", {
  box <- list(x1 = c(-2, 2), x2 = c(-2, 2))
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)

  # 0.001 on [-2, 2]^2 is 4001^2 = 16,008,001 points: refused, not attempted
  expect_error(
    compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                       design_box = box, step = 0.001, psi_star = r$psi_star),
    "max_points")

  # the message says how to get past it
  err <- tryCatch(
    compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                       design_box = box, step = 0.001, psi_star = r$psi_star),
    error = conditionMessage)
  expect_match(err, "criterion_only")
  expect_match(err, "coarser grid")

  # raising the cap is honoured (a grid this size is fine)
  ok <- compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                           design_box = box, step = 0.05,
                           psi_star = r$psi_star, max_points = 1e7)
  expect_true(is.finite(ok$max_d))

  # and a grid under the cap is unaffected
  small <- compound_criterion(r$support, r$weights, cmp3,
                              alpha = c(0.4, 0.35, 0.25),
                              design_box = box, step = 0.25,
                              psi_star = r$psi_star)
  expect_true(is.finite(small$max_d))
  expect_true(is.finite(small$criterion))
})

test_that("criterion_only scores a design without building any grid", {
  box <- list(x1 = c(-2, 2), x2 = c(-2, 2))
  r <- compound_design(cmp3, alpha = c(0.4, 0.35, 0.25), candidate_set = Xg)

  # no grid, so no cap to hit even at an absurd step
  co <- compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                           design_box = box, step = 0.001,
                           psi_star = r$psi_star, criterion_only = TRUE)
  expect_true(is.finite(co$criterion))
  expect_length(co$psi, 3L)
  # the equivalence-theorem quantities are absent, as documented
  expect_null(co$max_d)
  expect_null(co$sensitivity)
  expect_null(co$is_optimal)

  # and the criterion agrees with the full scan's
  full <- compound_criterion(r$support, r$weights, cmp3,
                             alpha = c(0.4, 0.35, 0.25),
                             design_box = box, step = 0.25,
                             psi_star = r$psi_star)
  expect_equal(co$criterion, full$criterion, tolerance = 1e-10)

  # without psi_star there are no reference values to be had: a clear error
  expect_error(
    compound_criterion(r$support, r$weights, cmp3, alpha = c(0.4, 0.35, 0.25),
                       design_box = box, step = 0.25, criterion_only = TRUE),
    "criterion_only")
})

test_that("per-component criterion values are reported on the optimal_design() scale", {
  al <- c(0.4, 0.35, 0.25)
  r  <- compound_design(cmp3, alpha = al, candidate_set = Xg)
  # D components: -log(Psi); A component: 1 / Psi
  expect_equal(unname(r$component_criterion[1:2]), unname(-log(r$psi[1:2])))
  expect_equal(unname(r$component_criterion[3]),   unname(1 / r$psi[3]))
  # the reference values are what optimal_design() reports for each component alone
  o1 <- suppressWarnings(optimal_design(info_vector = f1, theta = th1,
                                        candidate_set = Xg, p = 0))
  o3 <- suppressWarnings(optimal_design(info_vector = f3, theta = th3,
                                        candidate_set = Xg, p = 1, subset = c(2, 3)))
  expect_equal(unname(r$component_criterion_star[1]), o1$criterion, tolerance = 1e-5)
  expect_equal(unname(r$component_criterion_star[3]), o3$criterion, tolerance = 1e-5)
  # compound_criterion() reports the same values for the same design
  cc <- compound_criterion(r$support, r$weights, cmp3, alpha = al,
                           psi_star = r$psi_star, candidate_set = Xg)
  expect_equal(unname(cc$component_criterion), unname(r$component_criterion),
               tolerance = 1e-8)
  # the reported efficiency is the ratio to the reference times the reference's
  # certified bound; dividing the bound out recovers the ratio on either scale
  expect_equal(unname(r$efficiency_lower_bound[1] / r$reference_bound[1]),
               exp(r$component_criterion_star[[1]] - r$component_criterion[[1]]),
               tolerance = 1e-8)
  expect_equal(unname(r$efficiency_lower_bound[3] / r$reference_bound[3]),
               r$component_criterion_star[[3]] / r$component_criterion[[3]],
               tolerance = 1e-8)
})
