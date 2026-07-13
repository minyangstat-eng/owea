# Tests for fit_design(): the estimation half of simulate_design(), standalone.

test_that("fit_design reproduces simulate_design's fit exactly (identity)", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  r <- simulate_design(sup, counts = rep(40, 5), link = "identity",
                       design_box = list(dose = c(-1, 1)), x = 1,
                       theta = c(2, -3), sigma = 0.5, seed = 1)
  g <- fit_design(r$data, link = "identity", x = 1,
                  design_box = list(dose = c(-1, 1)))
  expect_equal(g$theta_hat, r$theta_hat)
  expect_equal(g$se, r$se)
  expect_equal(g$loglik, r$loglik)
  expect_equal(g$sigma_hat, r$sigma_hat)
  expect_equal(g$N, nrow(r$data))
})

test_that("fit_design reproduces simulate_design's fit exactly (logit, quadratic)", {
  sup <- matrix(c(-1, -0.5, 0, 0.5, 1), ncol = 1)
  r <- simulate_design(sup, counts = rep(200, 5), link = "logit", x = 1,
                       xx = list(c(1, 1)), design_box = list(c(-1, 1)),
                       theta = c(0.5, 1, -0.5), seed = 3)
  g <- fit_design(r$data, link = "logit", x = 1, xx = list(c(1, 1)),
                  design_box = list(c(-1, 1)))
  expect_true(g$converged)
  expect_equal(g$theta_hat, r$theta_hat)
  expect_equal(unname(g$vcov), unname(as.matrix(r$vcov)))
  expect_equal(dimnames(g$vcov), list(g$coef_names, g$coef_names))
})

test_that("fit_design reproduces simulate_design's fit (loglinear, cumulative)", {
  sup <- matrix(c(0, 0.5, 1), ncol = 1)
  r <- simulate_design(sup, counts = rep(60, 3), link = "loglinear", x = 1,
                       design_box = list(c(0, 1)), theta = c(0.2, 0.5), seed = 2)
  g <- fit_design(r$data, link = "loglinear", x = 1, design_box = list(c(0, 1)))
  expect_equal(g$theta_hat, r$theta_hat)

  r2 <- simulate_design(sup, counts = rep(80, 3), link = "cumulative", x = 1,
                        ncat = 3, design_box = list(c(0, 1)),
                        theta = c(-0.5, 0.8, 1.0), seed = 5)
  g2 <- fit_design(r2$data, link = "cumulative", x = 1, ncat = 3,
                   design_box = list(c(0, 1)))
  expect_equal(g2$theta_hat, r2$theta_hat)
  expect_equal(g2$coef_names, r2$coef_names)
})

test_that("fit_design handles a factor covariate column", {
  sup <- cbind(rep(1:3, each = 2), rep(c(-1, 1), 3))
  r <- simulate_design(sup, counts = rep(50, 6), link = "identity",
                       design_box = list(grp = 3, dose = c(-1, 1)),
                       f = 1, x = 1, theta = c(1, 0.5, -0.2, 2), seed = 4)
  # same data, but with the factor column supplied as an R factor
  dat <- r$data
  dat$grp <- factor(dat$grp, levels = 1:3, labels = c("a", "b", "c"))
  g <- fit_design(dat, link = "identity", f = 1, x = 1)
  expect_equal(unname(g$theta_hat), unname(r$theta_hat))
})

test_that("the response defaults to 'y' / the last column, and can be overridden", {
  sup <- matrix(c(-1, 0, 1), ncol = 1)
  r <- simulate_design(sup, counts = rep(30, 3), link = "identity", x = 1,
                       design_box = list(dose = c(-1, 1)), theta = c(1, 2),
                       sigma = 0.5, seed = 8)
  # last column, not named "y"
  d1 <- r$data
  names(d1) <- c("dose", "resp")
  g1 <- fit_design(d1, link = "identity", x = 1, design_box = list(c(-1, 1)))
  expect_equal(unname(g1$theta_hat), unname(r$theta_hat))

  # response in the FIRST column -- must be pointed at explicitly
  d2 <- r$data[, c("y", "dose")]
  g2 <- fit_design(d2, link = "identity", x = 1, design_box = list(c(-1, 1)))
  expect_equal(unname(g2$theta_hat), unname(r$theta_hat))   # found by name "y"
  names(d2) <- c("resp", "dose")
  g3 <- fit_design(d2, link = "identity", x = 1, response = "resp",
                   design_box = list(c(-1, 1)))
  expect_equal(unname(g3$theta_hat), unname(r$theta_hat))
  g4 <- fit_design(d2, link = "identity", x = 1, response = 1L,
                   design_box = list(c(-1, 1)))
  expect_equal(unname(g4$theta_hat), unname(r$theta_hat))
})

test_that("fit_design validates the response and the model spec", {
  d <- data.frame(dose = c(-1, 0, 1, -1, 0, 1), y = c(0, 1, 1, 0, 1, 0))
  expect_error(fit_design(d, x = 1), "link")
  expect_error(fit_design(d, link = "logit", x = 1, response = "z"), "not found")
  bad <- d; bad$y <- c(0, 1, 2, 3, 1, 0)
  expect_error(fit_design(bad, link = "logit", x = 1,
                          design_box = list(c(-1, 1))), "0/1")
  expect_error(fit_design(d, link = "cumulative", x = 1, ncat = 3,
                          design_box = list(c(-1, 1))), "1\\.\\.3")
})
