# Wall-clock time limits (R/timeout.R).  Only the Shiny app arms a limit, but
# the machinery lives in the package: the C++ engine's interrupt checks, the
# error handlers that re-raise a time-limit error instead of swallowing it, and
# .with_time_limit(), which turns both the R-level error and the C++-level
# interrupt into one "owea_timeout" error.
#
# Never expect_error() a raw setTimeLimit() around a C++ call: the limit then
# arrives as an interrupt, which is not an error and would abort testthat.

wl    <- owea:::.with_time_limit
clear <- function() setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0
}
# a bounded busy loop (seconds for the default n), so a broken limit fails the
# test instead of hanging it
busy <- function(n = 3e8) { s <- 0; for (i in 1:n) s <- s + 1; s }

covs <- function(k, step) lapply(seq_len(k), function(i)
  list(name = paste0("x", i), type = "continuous", lo = -1, hi = 1, step = step))

test_that(".with_time_limit stops R code and reports an owea_timeout error", {
  on.exit(clear(), add = TRUE)
  expect_identical(wl(1 + 1, Inf), 2)               # no limit: just the value
  expect_identical(wl(1 + 1, NULL), 2)
  expect_identical(wl(1 + 1, 0), 2)
  t <- elapsed(expect_error(wl(busy(), 0.2), class = "owea_timeout"))
  expect_lt(t, 1.5)
  expect_false(is.finite(owea:::.timeout_state$deadline))   # disarmed again
  # the message and the limit travel with the condition
  err <- tryCatch(wl(busy(), 0.2, "custom text"), error = function(e) e)
  expect_s3_class(err, "owea_timeout")
  expect_s3_class(err, "error")
  expect_identical(conditionMessage(err), "custom text")
  expect_equal(err$limit, 0.2)
  # any other error passes through unchanged
  expect_error(wl(stop("boom"), 5), "boom", class = "simpleError")
  # the limit is gone afterwards: a sleep longer than it succeeds
  expect_silent(Sys.sleep(0.3))
})

test_that("a handler that swallows errors does not swallow the time limit", {
  on.exit(clear(), add = TRUE)
  # the pattern used throughout the package: re-raise first, then the fallback.
  # Without it, R would clear the (one-shot) limit when it fires and the loop
  # would run on unlimited.
  guarded <- function(n) {
    s <- 0
    for (i in 1:400)
      s <- s + tryCatch(busy(n),
                        error = function(e) { owea:::.rethrow_if_timeout(e); 0 })
    s
  }
  t <- elapsed(expect_error(wl(guarded(1e6), 0.2), class = "owea_timeout"))
  expect_lt(t, 1.5)
  # with no limit armed the helper is a no-op ...
  expect_null(owea:::.rethrow_if_timeout(simpleError("boom")))
  # ... but it still recognises R's own time-limit error by its message
  expect_error(owea:::.rethrow_if_timeout(simpleError("reached elapsed time limit")),
               "time limit")
  expect_error(owea:::.rethrow_if_timeout(simpleError("reached CPU time limit")),
               "time limit")
})

test_that("the limit interrupts the C++ exchange loop", {
  skip_on_cran()
  on.exit(clear(), add = TRUE)
  # eps0 < 0 can never be met, so the exchange loop runs until max_iter: only
  # the Rcpp::checkUserInterrupt() calls inside the engine can stop it
  sp <- owea:::.ui_model_spec(covs(3, 0.2), link = "logit")
  a  <- owea:::.ui_solver_args(sp, "optimal", theta = c(0.5, 1, -1, 0.5), p = 0L)
  a$eps0 <- -1; a$max_iter <- 10000000L
  t <- elapsed(expect_error(
    wl(suppressMessages(suppressWarnings(do.call(optimal_design, a))), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
  expect_false(is.finite(owea:::.timeout_state$deadline))
  expect_silent(Sys.sleep(0.3))
})

test_that("the limit interrupts the compound engine", {
  skip_on_cran()
  on.exit(clear(), add = TRUE)
  f1 <- function(x, theta) { q <- c(1, x[1], x[2]); eta <- sum(q * theta)
                             (exp(eta / 2) / (1 + exp(eta))) * q }
  f2 <- function(x, theta) { q <- c(1, x[1], x[2], x[1] * x[2]); eta <- sum(q * theta)
                             (exp(eta / 2) / (1 + exp(eta))) * q }
  cmp <- list(list(info_vector = f1, theta = c(0.5, 1, -1), p = 0, name = "main"),
              list(info_vector = f2, theta = c(0.5, 1, -1, 0.5), p = 0, name = "inter"))
  t <- elapsed(expect_error(
    wl(suppressMessages(suppressWarnings(
      compound_design(components = cmp, alpha = c(0.5, 0.5),
                      design_box = list(c(-2, 2), c(-2, 2)), step_sequence = 0.5,
                      eps0 = -1, max_iter = 10000000L))), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
  expect_false(is.finite(owea:::.timeout_state$deadline))
})

test_that("the limit stops a large grid build (R level)", {
  skip_on_cran()
  on.exit(clear(), add = TRUE)
  sp <- owea:::.ui_model_spec(covs(3, 0.05), link = "logit")     # 41^3 points
  a  <- owea:::.ui_solver_args(sp, "optimal", theta = c(0.5, 1, -1, 0.5), p = 0L)
  t <- elapsed(expect_error(
    wl(suppressMessages(suppressWarnings(do.call(optimal_design, a))), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
})

test_that("the limit stops the continuous search (optim handlers re-raise it)", {
  skip_on_cran()
  on.exit(clear(), add = TRUE)
  # ~1 s unlimited on the development machine; the time goes into the R-level
  # polishing / multi-start loops, whose tryCatch(optim()) handlers must not
  # swallow the limit
  th <- c(0.5, 1, -1, 0.5, -0.5, 1, -1, 0.3, -0.3, 0.3)
  sp <- owea:::.ui_model_spec(covs(6, 0.5),
                              interactions = list(c(1, 2), c(3, 4), c(5, 6)),
                              link = "logit", continuous = TRUE, n_audit = 0L)
  a  <- owea:::.ui_solver_args(sp, "optimal", theta = th, p = 0L)
  t <- elapsed(expect_error(
    wl(suppressMessages(suppressWarnings(do.call(optimal_design, a))), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
  expect_false(is.finite(owea:::.timeout_state$deadline))
})

test_that("the limit stops a simulation study (per-replicate handlers re-raise it)", {
  skip_on_cran()
  on.exit(clear(), add = TRUE)
  sp <- owea:::.ui_model_spec(covs(1, 0.5), link = "logit")
  a  <- owea:::.ui_solver_args(sp, "exact", theta = c(0.5, 1), p = 0L, n = 20)
  d  <- suppressWarnings(do.call(exact_design, a))
  # the app's pooled path ...
  t <- elapsed(expect_error(
    wl(owea:::.ui_simulate(sp, theta = c(0.5, 1), sigma = 1, support = d$support,
                           counts = d$counts, nsim = 1000000L, seed = 1), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
  # ... and the exported simulate_design()
  sa <- owea:::.ui_solver_args(sp, "simulate", theta = c(0.5, 1))
  t <- elapsed(expect_error(
    wl(do.call(simulate_design, c(list(support = d$support, counts = d$counts,
                                       nsim = 1000000L, seed = 1), sa)), 0.3),
    class = "owea_timeout"))
  expect_lt(t, 1.5)
  expect_false(is.finite(owea:::.timeout_state$deadline))
  expect_silent(Sys.sleep(0.3))
})
