# Tests for simulate_design(): simulation + MLE fit + Monte Carlo.

test_that("seed makes simulate_design reproducible", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  a <- simulate_design(sup, counts = rep(50, 5), link = "logit", x = 1,
                       design_box = list(c(-1, 1)), theta = c(0.3, 1.0), seed = 7)
  b <- simulate_design(sup, counts = rep(50, 5), link = "logit", x = 1,
                       design_box = list(c(-1, 1)), theta = c(0.3, 1.0), seed = 7)
  expect_equal(a$data$y, b$data$y)
  expect_equal(a$theta_hat, b$theta_hat)
})

test_that("single-run structure and identity recovery", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  r <- simulate_design(sup, counts = rep(400, 5), link = "identity",
                       design_box = list(dose = c(-1, 1)), x = 1,
                       theta = c(2, -3), sigma = 0.5, seed = 1)
  expect_named(r$data, c("dose", "y"))
  expect_equal(nrow(r$data), 2000L)
  expect_length(r$theta_hat, 2L)
  expect_equal(names(r$theta_hat), c("(Intercept)", "dose"))
  expect_lt(max(abs(r$theta_hat - c(2, -3))), 0.1)      # well within a few SE
  expect_lt(abs(r$sigma_hat - 0.5), 0.05)
})

test_that("logit MLE recovers theta within a few standard errors", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  th  <- c(0.4, 1.2)
  r <- simulate_design(sup, counts = rep(3000, 5), link = "logit", x = 1,
                       design_box = list(c(-1, 1)), theta = th, seed = 3)
  expect_true(r$converged)
  expect_true(all(abs(r$theta_hat - th) < 4 * r$se))
})

test_that("fit = FALSE returns data only; nsim>1 requires fit", {
  sup <- matrix(c(0, 1), ncol = 1)
  d <- simulate_design(sup, counts = c(5, 5), link = "loglinear", x = 1,
                       design_box = list(c(0, 1)), theta = c(0.2, 0.5),
                       seed = 1, fit = FALSE)
  expect_true(is.null(d$theta_hat))
  expect_true("y" %in% names(d$data))
  expect_error(
    simulate_design(sup, counts = c(5, 5), link = "logit", x = 1,
                    design_box = list(c(0, 1)), theta = c(0, 1),
                    nsim = 5, fit = FALSE),
    "fit = FALSE")
})

test_that("theta length and cumulative threshold order are validated", {
  sup <- matrix(c(-1, 1), ncol = 1)
  expect_error(
    simulate_design(sup, counts = c(5, 5), link = "logit", x = 1,
                    design_box = list(c(-1, 1)), theta = c(1, 2, 3)),
    "must have length")
  expect_error(
    simulate_design(sup, counts = c(5, 5), link = "cumulative", ncat = 3, x = 1,
                    design_box = list(c(-1, 1)), theta = c(0.5, -0.5, 1)),
    "strictly increasing")
})

test_that("Monte Carlo: empirical covariance tracks the design-based M^{-1}/N", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  th  <- c(0.3, 1.0)
  mc <- simulate_design(sup, counts = rep(200, 5), link = "logit", x = 1,
                        design_box = list(c(-1, 1)), theta = th,
                        nsim = 500, seed = 11)
  expect_equal(dim(mc$estimates), c(500L, 2L))
  expect_equal(mc$n_converged, 500L)
  # empirical vs design-predicted standard errors agree to ~15%
  expect_true(all(mc$se_empirical / mc$se_design > 0.8 &
                  mc$se_empirical / mc$se_design < 1.2))
  # estimator is ~unbiased: |bias| small relative to design SE
  expect_true(all(abs(mc$bias) < 0.15 * mc$se_design * sqrt(mc$nsim)))
})

test_that("multinomial fit matches nnet::multinom log-likelihood", {
  skip_if_not_installed("nnet")
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  th  <- c(0.3, 1.0, -0.4, 0.8)
  r <- simulate_design(sup, counts = rep(400, 5), link = "multinomial",
                       ncat = 3, x = 1, design_box = list(c(-1, 1)),
                       theta = th, seed = 5)
  expect_true(all(r$data$y %in% 1:3))
  fit <- nnet::multinom(y ~ x1, data = transform(r$data, x1 = r$data[[1]]),
                        trace = FALSE)
  expect_equal(unname(r$loglik), as.numeric(logLik(fit)), tolerance = 1e-3)
})

test_that("cumulative fit matches MASS::polr log-likelihood", {
  skip_if_not_installed("MASS")
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  th  <- c(-0.8, 0.6, 1.1)                    # (alpha_1 < alpha_2, beta)
  r <- simulate_design(sup, counts = rep(400, 5), link = "cumulative",
                       ncat = 3, x = 1, design_box = list(c(-1, 1)),
                       theta = th, seed = 6)
  expect_true(all(r$data$y %in% 1:3))
  df <- data.frame(y = factor(r$data$y, ordered = TRUE), x1 = r$data[[1]])
  fit <- MASS::polr(y ~ x1, data = df, method = "logistic", Hess = FALSE)
  expect_equal(unname(r$loglik), as.numeric(logLik(fit)), tolerance = 1e-3)
})

test_that("simulate_design accepts an exact_design object", {
  ex <- suppressWarnings(exact_design(n = 30, design_box = list(c(-1, 1)),
                                      step_sequence = c(0.5, 0.25),
                                      link = "logit", x = 1,
                                      theta = c(0.2, 1.0), seed = 1))
  r <- simulate_design(ex, link = "logit", x = 1, design_box = list(c(-1, 1)),
                       theta = c(0.2, 1.0), seed = 2)
  expect_equal(r$N, 30L)
  expect_length(r$theta_hat, 2L)
})
