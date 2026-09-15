# Every reported `efficiency` is a guaranteed lower bound relative to the TRUE
# optimum (Becker & Yang, "Post Hoc Control Group Selection via Constrained
# Optimal Design", Theorems 4.5 and 4.6), in owea's notation:
#     D: eff >= 1 / (1 + max_d / v),     A: eff >= 1 - max_d / criterion.
# owea's max_d is their optimality gap E(xi) for the unconstrained problem.

th <- c(0.5, 1, -1)
f  <- function(x, theta) { q <- c(1, x[1], x[2]); e <- sum(q * theta)
                           (exp(e / 2) / (1 + exp(e))) * q }          # logistic, k = 3
Xg <- candidate_grid(list(c(-1, 1), c(-1, 1)), 0.25)

true_eff <- function(p, crit_opt, crit) if (p == 0) exp(crit_opt - crit) else crit_opt / crit

test_that("a converged design reports (essentially) full efficiency", {
  for (p in c(0, 1)) {
    r <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = p,
                                         candidate_set = Xg))
    expect_true(r$converged)
    expect_gte(r$efficiency_lower_bound, 1 - 1e-5)
    expect_lte(r$efficiency_lower_bound, 1)
  }
})

test_that("the reported efficiency is a valid lower bound for an unconverged design", {
  for (p in c(0, 1)) for (sub in list(NULL, 2:3)) {
    opt  <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = p,
                                            candidate_set = Xg, subset = sub))
    poor <- suppressWarnings(optimal_design(info_vector = f, theta = th, p = p,
                                            candidate_set = Xg, subset = sub,
                                            max_iter = 2, auto_warm_start = FALSE,
                                            init_method = "random"))
    eff <- true_eff(p, opt$criterion, poor$criterion)   # vs the (certified) optimum
    expect_true(is.finite(poor$efficiency_lower_bound))
    expect_gt(poor$efficiency_lower_bound, 0)
    expect_lte(poor$efficiency_lower_bound, eff + 1e-10)             # a valid lower bound ...
    expect_lte(eff, 1 + 1e-10)                            # ... on a true efficiency
    v <- if (is.null(sub)) 3 else length(sub)
    expected <- if (p == 0) 1 / (1 + max(poor$max_d, 0) / v)
                else        max(0, 1 - max(poor$max_d, 0) / poor$criterion)
    expect_equal(poor$efficiency_lower_bound, expected)
    # verify_optimality() reports the same value for the same design and grid
    vo <- verify_optimality(poor$support, poor$weights, info_vector = f, theta = th,
                            p = p, subset = sub, candidate_set = Xg)
    expect_equal(vo$efficiency_lower_bound, poor$efficiency_lower_bound, tolerance = 1e-8)
  }
})

test_that("exact designs report the certified efficiency (ratio x reference bound)", {
  for (p in c(0, 1)) {
    e <- suppressWarnings(exact_design(n = 12, info_vector = f, theta = th, p = p,
                                       candidate_set = Xg, seed = 1))
    expect_true(is.finite(e$efficiency_lower_bound))
    expect_gt(e$efficiency_lower_bound, 0)
    expect_lte(e$efficiency_lower_bound, 1 + 1e-12)
    expect_gt(e$approx$efficiency_bound, 0)
    expect_lte(e$approx$efficiency_bound, 1)
    expect_equal(e$efficiency_lower_bound,
                 min(true_eff(p, e$approx$criterion, e$criterion), 1) * e$approx$efficiency_bound)
  }
})

test_that("compound designs report certified component efficiencies", {
  th2 <- c(0.5, 1, -1, 0.5)
  f2  <- function(x, theta) { q <- c(1, x[1], x[2], x[1] * x[2]); e <- sum(q * theta)
                              (exp(e / 2) / (1 + exp(e))) * q }
  cmp <- list(list(info_vector = f,  theta = th,  p = 0, name = "D"),
              list(info_vector = f2, theta = th2, p = 1, subset = 2:4, name = "A"))
  r <- compound_design(cmp, alpha = c(0.5, 0.5), candidate_set = Xg)
  expect_true(all(is.finite(r$reference_bound)))
  expect_true(all(r$reference_bound > 0 & r$reference_bound <= 1))
  expect_true(all(r$efficiency_lower_bound > 0 & r$efficiency_lower_bound <= 1 + 1e-12))
  # efficiency = (Psi / Psi*) x reference bound
  expect_equal(unname(r$efficiency_lower_bound), unname(pmin(r$psi / r$psi_star, 1) * r$reference_bound))
  expect_true(is.finite(r$criterion_bound) && r$criterion_bound > 0 && r$criterion_bound <= 1)
  # cross table: each column certified by its component's bound (diagonal = bound)
  expect_equal(unname(diag(r$cross_efficiency[1:2, ])), unname(r$reference_bound), tolerance = 1e-8)
  # compound_criterion() picks the bounds up from the compound_design object
  cc <- compound_criterion(r, components = cmp, candidate_set = Xg)
  expect_equal(unname(cc$reference_bound), unname(r$reference_bound))
  expect_equal(unname(cc$efficiency_lower_bound), unname(r$efficiency_lower_bound), tolerance = 1e-8)
})
