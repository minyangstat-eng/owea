# Tests for the formula-style model builder (model.R).

test_that("factor contrast coding matches the spec (zero-sum & baseline)", {
  # 2-level factor, intercept + main effect
  z2 <- model_info_vector(design_box = list(c(2)), link = "identity", f = 1)
  expect_equal(z2(1), c(1,  1))
  expect_equal(z2(2), c(1, -1))
  b2 <- model_info_vector(design_box = list(c(2)), link = "identity", f = 1,
                          coding = "baseline")
  expect_equal(b2(1), c(1, 1))
  expect_equal(b2(2), c(1, 0))

  # 3-level factor
  z3 <- model_info_vector(design_box = list(c(3)), link = "identity", f = 1)
  expect_equal(z3(1), c(1, 1, 0))
  expect_equal(z3(2), c(1, 0, 1))
  expect_equal(z3(3), c(1, -1, -1))
  b3 <- model_info_vector(design_box = list(c(3)), link = "identity", f = 1,
                          coding = "baseline")
  expect_equal(b3(3), c(1, 0, 0))

  # 4-level factor
  z4 <- model_info_vector(design_box = list(c(4)), link = "identity", f = 1)
  expect_equal(z4(1), c(1, 1, 0, 0))
  expect_equal(z4(4), c(1, -1, -1, -1))
  b4 <- model_info_vector(design_box = list(c(4)), link = "identity", f = 1,
                          coding = "baseline")
  expect_equal(b4(4), c(1, 0, 0, 0))
})

test_that("term order and encodings (f, x, fx, xx) are correct", {
  # design: factors 1,2 (levels 2,3); continuous 1,2
  db <- list(c(2), c(3), c(0, 6), c(0, 6))
  iv <- model_info_vector(design_box = db, link = "identity",
                          f = c(1, 2), x = c(1, 2),
                          fx = c(11, 22), xx = c(11, 12))
  x <- c(1, 2, 3, 4)          # factor1=1, factor2=2, cont1=3, cont2=4
  # order: intercept, f(fac1 [1col], fac2 [2col]), x(c1,c2),
  #        fx(fac1 x c1 [1col], fac2 x c2 [2col]), xx(c1^2, c1*c2)
  expect_equal(iv(x), c(1,                # intercept
                        1,                # fac1 level1 zero-sum (L=2)
                        0, 1,             # fac2 level2 zero-sum (L=3)
                        3, 4,             # c1, c2
                        1 * 3,            # fac1 x c1
                        0 * 4, 1 * 4,     # fac2 x c2
                        3 * 3, 3 * 4))    # c1^2, c1*c2
  # explicit pair-list form equals the two-digit form
  iv2 <- model_info_vector(design_box = db, link = "identity",
                           f = c(1, 2), x = c(1, 2),
                           fx = list(c(1, 1), c(2, 2)),
                           xx = list(c(1, 1), c(1, 2)))
  expect_equal(iv2(x), iv(x))
})

test_that("logit / loglinear link formulas match hand-coded versions", {
  db <- list(c(2), c(-1, 1))
  th <- c(0.3, -0.7, 1.1, 0.5)
  hand_logit <- function(x, theta) {
    fc  <- if (x[1] == 1) 1 else -1
    q   <- c(1, fc, x[2], fc * x[2])
    eta <- min(max(sum(q * theta), -500), 500)
    (exp(eta / 2) / (1 + exp(eta))) * q
  }
  hand_pois <- function(x, theta) {
    fc  <- if (x[1] == 1) 1 else -1
    q   <- c(1, fc, x[2], fc * x[2])
    eta <- min(max(sum(q * theta), -500), 500)
    exp(eta / 2) * q
  }
  ivl <- model_info_vector(design_box = db, link = "logit",
                           f = 1, x = 1, fx = c(11))
  ivp <- model_info_vector(design_box = db, link = "loglinear",
                           f = 1, x = 1, fx = c(11))
  pts <- list(c(1, -1), c(1, 0.5), c(2, -0.3), c(2, 1))
  for (x in pts) {
    expect_equal(ivl(x, th), hand_logit(x, th))
    expect_equal(ivp(x, th), hand_pois(x, th))
  }
})

test_that("optimal_design via model spec == hand-coded (identity, mixed factor)", {
  mix_vec <- function(x) {
    contr <- switch(as.character(x[1]),
                    "1" = c(1, 0), "2" = c(0, 1), "3" = c(-1, -1))
    c(1, contr, x[2])
  }
  common <- list(design_box = list(c(3), c(0, 1)),
                 step_sequence = c(0.1, 0.05), p = 0)
  r_hand <- suppressWarnings(do.call(optimal_design,
                c(list(info_vector = mix_vec), common)))
  r_spec <- suppressWarnings(do.call(optimal_design,
                c(list(link = "identity", f = 1, x = 1), common)))
  expect_equal(r_spec$criterion, r_hand$criterion, tolerance = 1e-6)
  ord <- function(r) r$support[do.call(order, as.data.frame(r$support)), ,
                               drop = FALSE]
  expect_equal(unname(ord(r_spec)), unname(ord(r_hand)), tolerance = 1e-4)
})

test_that("optimal_design via model spec == hand-coded (logit, factor x cont)", {
  th <- c(0.5, -0.8, 1.2, 0.4)
  hand_logit <- function(x, theta) {
    fc  <- if (x[1] == 1) 1 else -1
    q   <- c(1, fc, x[2], fc * x[2])
    eta <- min(max(sum(q * theta), -500), 500)
    (exp(eta / 2) / (1 + exp(eta))) * q
  }
  common <- list(design_box = list(c(2), c(-1, 1)),
                 step_sequence = c(0.2, 0.1, 0.05), p = 1)
  r_hand <- suppressWarnings(do.call(optimal_design,
                c(list(info_vector = hand_logit, theta = th), common)))
  r_spec <- suppressWarnings(do.call(optimal_design,
                c(list(link = "logit", f = 1, x = 1, fx = c(11), theta = th),
                  common)))
  expect_equal(r_spec$criterion, r_hand$criterion, tolerance = 1e-6)
})

test_that("missing theta for a GLM draws N(0,1) with a warning and returns it", {
  set.seed(1)
  expect_warning(
    r <- optimal_design(design_box = list(c(2), c(-1, 1)),
                        step_sequence = c(0.2, 0.1),
                        link = "logit", f = 1, x = 1, fx = c(11), p = 0),
    "drawn from N\\(0,1\\)")
  expect_length(r$theta, 4L)          # k = intercept + f + x + fx
  expect_true(r$converged)
})

test_that("identity link needs no theta; spec + explicit model is an error", {
  expect_silent(
    r <- suppressWarnings(optimal_design(design_box = list(c(0, 1)),
              step_sequence = c(0.2, 0.1), link = "identity", x = 1, p = 0)))
  expect_true(r$converged)
  expect_error(
    optimal_design(design_box = list(c(0, 1)), step_sequence = 0.1,
                   link = "identity", x = 1,
                   info_vector = function(x) c(1, x)),
    "not both")
})

test_that("design_information accepts a model spec (== hand-coded matrix)", {
  th <- c(0.5, -0.8, 1.2, 0.4)
  hand_logit <- function(x, theta) {
    fc  <- if (x[1] == 1) 1 else -1
    q   <- c(1, fc, x[2], fc * x[2])
    eta <- min(max(sum(q * theta), -500), 500)
    (exp(eta / 2) / (1 + exp(eta))) * q
  }
  hand_mat <- function(x, theta) tcrossprod(hand_logit(x, theta))
  design  <- rbind(c(1, -1), c(1, 1), c(2, -1), c(2, 1))
  weights <- c(0.25, 0.25, 0.25, 0.25)
  M_hand <- design_information(design, weights, info_matrix = hand_mat, theta = th)
  M_spec <- design_information(design, weights, link = "logit",
                              f = 1, x = 1, fx = c(11), theta = th,
                              factor_levels = c(2, NA))
  expect_equal(M_spec, M_hand, tolerance = 1e-10)

  # infor_matrix (vector-model helper) agrees too
  M_iv <- infor_matrix(design, weights, link = "logit",
                       f = 1, x = 1, fx = c(11), theta = th,
                       factor_levels = c(2, NA))
  expect_equal(M_iv, M_hand, tolerance = 1e-10)
})

test_that("design_information auto-detects factor columns (with a message)", {
  design  <- rbind(c(1, -1), c(1, 1), c(2, -1), c(2, 1))
  weights <- rep(0.25, 4)
  expect_message(
    design_information(design, weights, link = "logit",
                       f = 1, x = 1, fx = c(11), theta = c(0, 0, 0, 0)),
    "auto-detected factor")
})

test_that("factor x factor interaction (ff) matches hand-coded (Kronecker)", {
  # factors 1 (L=2) and 2 (L=3); intercept + both mains + interaction
  db <- list(c(2), c(3))
  iv <- model_info_vector(design_box = db, link = "identity",
                          f = c(1, 2), ff = c(12))
  hand <- function(x) {
    a <- if (x[1] == 1) 1 else -1
    b <- switch(as.character(x[2]),
                "1" = c(1, 0), "2" = c(0, 1), "3" = c(-1, -1))
    c(1, a, b, a * b[1], a * b[2])       # intercept, f(fac1,fac2), ff(fac1:fac2)
  }
  for (l1 in 1:2) for (l2 in 1:3) {
    x <- c(l1, l2)
    expect_equal(iv(x), hand(x))
  }
  # explicit pair-list form agrees
  iv2 <- model_info_vector(design_box = db, link = "identity",
                           f = c(1, 2), ff = list(c(1, 2)))
  expect_equal(iv2(c(2, 3)), iv(c(2, 3)))
})

test_that("ff rejects a factor interacting with itself", {
  expect_error(
    model_info_vector(design_box = list(c(2), c(3)), link = "identity",
                      f = c(1, 2), ff = c(11)),
    "DIFFERENT factors")
})

test_that("optimal_design with ff == hand-coded info_vector", {
  hand <- function(x) {
    a <- if (x[1] == 1) 1 else -1
    b <- switch(as.character(x[2]),
                "1" = c(1, 0), "2" = c(0, 1), "3" = c(-1, -1))
    c(1, a, b, a * b[1], a * b[2])
  }
  common <- list(design_box = list(c(2), c(3)),
                 step_sequence = numeric(0), p = 0)   # all-factor
  r_hand <- suppressWarnings(do.call(optimal_design,
                c(list(info_vector = hand), common)))
  r_spec <- suppressWarnings(do.call(optimal_design,
                c(list(link = "identity", f = c(1, 2), ff = c(12)), common)))
  expect_equal(r_spec$criterion, r_hand$criterion, tolerance = 1e-6)
})

test_that("all-factor design_box needs no step_sequence (omitted or numeric(0))", {
  box <- as.list(rep(2, 7))                 # 7 two-level factors
  r_omit <- suppressWarnings(
    optimal_design(design_box = box, link = "logit", f = 1:7, p = 1))
  expect_true(r_omit$converged)
  r_n0 <- suppressWarnings(
    optimal_design(design_box = box, step_sequence = numeric(0),
                   link = "logit", f = 1:7, p = 1))
  expect_true(r_n0$converged)
  # a continuous covariate without step_sequence is still an error
  expect_error(
    optimal_design(design_box = list(c(2), c(0, 1)), link = "identity", x = 1),
    "step_sequence")
})

test_that("loglinear link builds the Poisson info vector", {
  iv <- model_info_vector(design_box = list(c(0, 1)), link = "loglinear", x = 1)
  q  <- c(1, 0.5); th <- c(0.3, -0.4)
  eta <- sum(q * th)
  expect_equal(iv(c(0.5), th), exp(eta / 2) * q)
})

test_that("coef_names label each term in order (incl. names & ff)", {
  # unnamed -> generic F#/X#
  iv <- model_info_vector(design_box = list(c(2), c(3), c(0, 1)),
                          link = "identity", f = c(1, 2), x = 1,
                          fx = c(21), ff = c(12))
  expect_equal(attr(iv, "coef_names"),
               c("(Intercept)", "F1", "F2.1", "F2.2", "X1",
                 "F2.1:X1", "F2.2:X1", "F1:F2.1", "F1:F2.2"))
  # named design_box -> use names
  ivn <- model_info_vector(design_box = list(A = c(2), B = c(3), z = c(0, 1)),
                           link = "identity", f = c(1, 2), x = 1,
                           xx = c(11))
  expect_equal(attr(ivn, "coef_names"),
               c("(Intercept)", "A", "B.1", "B.2", "z", "z^2"))
})

test_that("model_summary prints labelled terms and returns a data.frame", {
  th <- c(0.5, -0.8, 1.2, 0.4)
  r <- suppressWarnings(optimal_design(
         design_box = list(f = c(2), d = c(-1, 1)),
         step_sequence = c(0.2, 0.1),
         link = "logit", f = 1, x = 1, fx = c(11), theta = th, p = 0))
  expect_output(s <- model_summary(r), "logit link")
  expect_s3_class(s, "data.frame")
  expect_equal(s$term, c("(Intercept)", "f", "d", "f:d"))
  expect_equal(s$theta, th)

  # works on a model function too (via attribute)
  iv <- model_info_vector(design_box = list(c(2), c(-1, 1)), link = "logit",
                          f = 1, x = 1, fx = c(11))
  expect_output(model_summary(iv), "term")

  # errors on a hand-coded model (no labels)
  r2 <- suppressWarnings(optimal_design(info_vector = function(x) c(1, x),
                                        design_box = list(c(0, 1)),
                                        step_sequence = c(0.2, 0.1), p = 0))
  expect_error(model_summary(r2), "no model-spec coefficient labels")
})

# ---- multinomial & cumulative (ordinal) logistic ----------------------------

# Independent per-point Fisher information for ANY multinomial-type model with
# cell probabilities pi(theta): M = sum_c (1/pi_c) (d pi_c/d theta)(.)' , with
# the Jacobian obtained by central finite differences (no reuse of the builder).
.num_info <- function(prob_fn, theta, J, h = 1e-6) {
  k <- length(theta); p0 <- prob_fn(theta)
  Jac <- matrix(0, J, k)
  for (i in seq_len(k)) {
    tp <- theta; tp[i] <- tp[i] + h
    tm <- theta; tm[i] <- tm[i] - h
    Jac[, i] <- (prob_fn(tp) - prob_fn(tm)) / (2 * h)
  }
  M <- matrix(0, k, k)
  for (c in seq_len(J)) M <- M + tcrossprod(Jac[c, ]) / p0[c]
  M
}

test_that("multinomial info matrix == independent finite-difference info (J=3)", {
  db <- list(c(0, 1)); J <- 3; xval <- 0.4
  frow <- c(1, xval)                              # f(x) = (1, x), p = 2
  th <- c(0.5, -0.3, -0.6, 0.9)                   # (beta_1, beta_2)
  im <- model_info_matrix(design_box = db, link = "multinomial", x = 1, ncat = J)
  mn_prob <- function(theta) {
    B <- matrix(theta, 2, J - 1); eta <- as.numeric(crossprod(frow, B))
    ex <- exp(eta); c(ex, 1) / (1 + sum(ex))
  }
  expect_equal(im(xval, th), .num_info(mn_prob, th, J), tolerance = 1e-5)
})

test_that("cumulative info matrix == independent finite-difference info (J=3)", {
  db <- list(c(0, 1)); J <- 3; xval <- 0.4
  ftil <- c(xval)                                 # non-intercept part, q = 1
  th <- c(-0.4, 0.6, 1.1)                         # (alpha_1 < alpha_2, beta)
  im <- model_info_matrix(design_box = db, link = "cumulative", x = 1, ncat = J)
  cum_prob <- function(theta) {
    alpha <- theta[1:2]; beta <- theta[3]
    gamma <- plogis(alpha + sum(ftil * beta))
    diff(c(0, gamma, 1))
  }
  expect_equal(im(xval, th), .num_info(cum_prob, th, J), tolerance = 1e-5)
})

test_that("J=2 multinomial and cumulative reduce exactly to logit", {
  db <- list(c(0, 1)); th <- c(0.3, -0.7)
  lg  <- model_info_matrix(design_box = db, link = "logit", x = 1)
  mn2 <- model_info_matrix(design_box = db, link = "multinomial", x = 1, ncat = 2)
  cm2 <- model_info_matrix(design_box = db, link = "cumulative", x = 1, ncat = 2)
  for (xv in c(-0.2, 0.4, 0.9)) {
    expect_equal(mn2(xv, th), lg(xv, th))
    expect_equal(cm2(xv, th), lg(xv, th))
  }
})

test_that("optimal_design converges for multinomial and cumulative models", {
  r_mn <- suppressWarnings(optimal_design(
    design_box = list(c(-1, 1)), step_sequence = c(0.2, 0.1),
    link = "multinomial", ncat = 3, x = 1,
    theta = c(0.5, -0.3, -0.5, 0.8), p = 0))
  expect_true(r_mn$converged)
  expect_equal(r_mn$coef_names,
               c("cat1:(Intercept)", "cat1:X1", "cat2:(Intercept)", "cat2:X1"))

  r_cm <- suppressWarnings(optimal_design(
    design_box = list(c(-1, 1)), step_sequence = c(0.2, 0.1),
    link = "cumulative", ncat = 3, x = 1,
    theta = c(-0.5, 0.5, 1.0), p = 0))
  expect_true(r_cm$converged)
  expect_equal(r_cm$coef_names, c("(threshold 1)", "(threshold 2)", "X1"))
})

test_that("cumulative thresholds must be strictly increasing; theta length checked", {
  expect_error(
    optimal_design(design_box = list(c(-1, 1)), step_sequence = 0.2,
                   link = "cumulative", ncat = 3, x = 1,
                   theta = c(0.5, -0.5, 1.0)),
    "strictly increasing")
  expect_error(
    optimal_design(design_box = list(c(-1, 1)), step_sequence = 0.2,
                   link = "multinomial", ncat = 3, x = 1, theta = c(1, 2, 3)),
    "has length")
})

test_that("multi-category auto-theta warns and sorts cumulative thresholds", {
  set.seed(1)
  expect_warning(
    r <- optimal_design(design_box = list(c(-1, 1)), step_sequence = 0.2,
                        link = "cumulative", ncat = 4, x = 1, p = 0),
    "thresholds sorted")
  a <- r$theta[1:3]
  expect_true(all(diff(a) > 0))
})

test_that("multi-category guards: ncat required, no info-vector form", {
  expect_error(model_info_matrix(design_box = list(c(0, 1)),
                                 link = "multinomial", x = 1), "ncat")
  expect_error(model_info_vector(design_box = list(c(0, 1)),
                                 link = "cumulative", x = 1, ncat = 3),
               "no information-vector form")
  expect_error(infor_matrix(matrix(0.5, 1, 1), link = "multinomial",
                            x = 1, ncat = 3, theta = c(0, 0, 0, 0)),
               "design_information")
})

test_that("design_information accepts a multinomial model spec", {
  th <- c(0.5, -0.3, -0.6, 0.9)
  design <- matrix(c(-0.5, 0.5), ncol = 1); w <- c(0.5, 0.5)
  im <- model_info_matrix(design_box = list(c(0, 1)), link = "multinomial",
                          x = 1, ncat = 3)
  M_ref <- 0.5 * im(-0.5, th) + 0.5 * im(0.5, th)
  M <- design_information(design, w, link = "multinomial", x = 1, ncat = 3,
                          theta = th)
  expect_equal(M, M_ref, tolerance = 1e-10)
})
