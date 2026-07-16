# Tests for the Shiny wizard's server logic (inst/shiny/owea-app/app.R): the
# branching, and that each branch reaches the solver with the right arguments.
# The pure translators behind it are tested in test-app.R.

app_env <- function() {
  dir <- system.file("shiny", "owea-app", package = "owea")
  testthat::skip_if(!nzchar(dir) || !file.exists(file.path(dir, "app.R")),
                    "the bundled app is not on disk")
  e <- new.env(parent = globalenv())
  sys.source(file.path(dir, "app.R"), envir = e)
  e
}

# one continuous covariate, so the model step is satisfied by the defaults
set_model <- function(session, link = "identity", cov_step_1 = "0.25", ...) {
  session$setInputs(link = link, ncov = 1, cov_name_1 = "dose",
                    cov_type_1 = "continuous", cov_lo_1 = -1, cov_hi_1 = 1,
                    cov_step_1 = cov_step_1, crit = "0", qoi = "all", ...)
}

test_that("the wizard skips the theta step for the identity link", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none")
    expect_false("theta" %in% steps())
    expect_equal(cur(), "model")

    session$setInputs(next_btn = 1)                  # model -> start
    expect_equal(cur(), "start")
    session$setInputs(next_btn = 2)                  # start -> criterion (no theta)
    expect_equal(cur(), "criterion")
    session$setInputs(back_btn = 1)
    expect_equal(cur(), "start")

    session$setInputs(link = "logit")                # now theta applies
    expect_true("theta" %in% steps())
    session$setInputs(next_btn = 3)
    expect_equal(cur(), "theta")
  })
})

test_that("the existing-design branch reaches the solver with n0, n1 and n = n1", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "design",
              exist_text = "dose,count\n-1,5\n1,5", exist_n0 = 10,
              design_type = "exact", n_new = 20, seed = 1)
    expect_true("design_in" %in% steps())
    ex <- existing()
    expect_equal(nrow(ex$points), 2L)
    expect_equal(ex$weights, c(0.5, 0.5))
    expect_equal(as.integer(ex$n0), 10L)

    session$setInputs(compute = 1)
    cc <- computed()
    expect_null(cc$error)
    expect_equal(cc$args$n0, 10L)
    expect_equal(cc$args$n, 20L)
    expect_equal(cc$args$n1, cc$args$n)            # tied, so no n1 != n warning
    expect_equal(sum(cc$args$xi0_weights), 1)
    expect_equal(sum(cc$res$counts), 20L)
    expect_equal(cur(), "results")                 # a good compute advances
  })
})

test_that("pasting/uploading the data set loads it -- no button click needed", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    dat <- paste("dose,y", "-1,0", "-1,0", "-1,1", "1,1", "1,1", "0,0",
                 sep = "\n")
    set_model(session, link = "logit", start = "data", data_text = dat,
              data_response = "y", use_cov = "yes", use_theta = "yes")
    session$elapse(1000)                    # let the debounce fire

    # data_load was NEVER clicked
    expect_false(is.null(rv$data))
    expect_equal(nrow(rv$data), 6L)
    expect_false(is.null(rv$data_existing))
    expect_false(is.null(rv$data_fit))
    expect_null(step_error("data_in"))      # Next is not blocked
  })
})

test_that("the existing-data branch fits theta and reuses the covariates", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    dat <- paste("dose,y", "-1,0", "-1,0", "-1,1", "1,1", "1,1", "0,0",
                 sep = "\n")
    set_model(session, link = "logit", start = "data", data_text = dat,
              data_response = "y", use_cov = "yes", use_theta = "yes")
    session$setInputs(data_load = 1)

    expect_equal(nrow(rv$data), 6L)
    ex <- rv$data_existing
    expect_equal(ex$n0, 6L)                        # n0 = number of observations
    expect_equal(nrow(ex$points), 3L)              # -1, 1, 0
    fit <- rv$data_fit
    expect_length(fit$theta_hat, 2L)
    # the fitted estimates become the assumed values in the theta boxes
    expect_equal(rv$theta_prefill, as.numeric(fit$theta_hat))
  })
})

test_that("the simulation study pools the existing stage with the new runs", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()

  # (i) an existing DESIGN: its responses are unknown, so they are simulated
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "design",
              exist_text = "dose,count\n-1,15\n1,15", exist_n0 = 30,
              design_type = "exact", n_new = 10, seed = 1)
    session$setInputs(compute = 1)
    expect_null(computed()$error)
    session$setInputs(sim_open = 1, sim_theta_1 = 1, sim_theta_2 = 2,
                      sim_sigma = 1, sim_nsim = 20, sim_seed = 1)
    expect_null(sim_obs())                      # no observed responses to reuse
    session$setInputs(run_sim = 1)
    s <- rv$sim$exact
    expect_false(inherits(s, "error"))
    expect_equal(s$N, 40L)                      # 30 existing + 10 new, pooled
    expect_equal(s$n_existing, 30L)
    expect_equal(s$n_new, 10L)
  })

  # (ii) an existing DATA SET: its real responses are kept, only the new runs
  # are simulated
  shiny::testServer(e$server, {
    dat <- paste("dose,y", "-1,0.4", "-1,1.1", "0,2.0", "0,1.4", "1,3.2",
                 "1,2.6", sep = "\n")
    set_model(session, link = "identity", start = "data", data_text = dat,
              data_response = "y", use_cov = "yes", use_theta = "no",
              design_type = "exact", n_new = 8, seed = 1)
    session$setInputs(data_load = 1)
    session$setInputs(compute = 1)
    expect_null(computed()$error)
    session$setInputs(sim_open = 1, sim_theta_1 = 1, sim_theta_2 = 2,
                      sim_sigma = 1, sim_nsim = 20, sim_seed = 1)
    expect_false(is.null(sim_obs()))            # the observed data IS the stage
    session$setInputs(run_sim = 1)
    s <- rv$sim$exact
    expect_false(inherits(s, "error"))
    expect_equal(s$N, 14L)                      # 6 observed + 8 new, pooled
    expect_equal(s$n_existing, 6L)
    expect_equal(s$n_new, 8L)
  })
})

test_that("a huge grid blocks Next on the model step until acknowledged", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none",
              cov_step_1 = "0.000001")          # (2 / 1e-6) + 1 > 1e6 points
    expect_gt(grid_sizes()[1], 1e6)
    session$setInputs(next_btn = 1)             # blocked: warning, no advance
    expect_equal(cur(), "model")
    session$setInputs(grid_proceed = 1)         # "Proceed anyway"
    expect_equal(cur(), "start")
    session$setInputs(back_btn = 1)             # back to the model step
    expect_equal(cur(), "model")
    session$setInputs(next_btn = 2)             # unchanged grid: no re-prompt
    expect_equal(cur(), "start")
    session$setInputs(cov_step_1 = "0.0000005") # a CHANGED grid re-arms it
    session$setInputs(back_btn = 2)
    session$setInputs(next_btn = 3)
    expect_equal(cur(), "model")
  })
})

test_that("a step sequence gates on its FIRST (coarsest) grid only", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    # fine finest step, but a coarse FIRST step: passes without a warning
    set_model(session, link = "identity", start = "none",
              cov_step_1 = "0.5, 0.001")
    expect_equal(length(spec()$step_sequence), 2L)
    expect_equal(grid_sizes(), c(5, 2001))      # coarsest first
    session$setInputs(next_btn = 1)
    expect_equal(cur(), "start")

    # a sequence whose FIRST grid is already huge is blocked
    session$setInputs(cov_step_1 = "0.000001, 0.0000005")
    session$setInputs(back_btn = 1)
    session$setInputs(next_btn = 2)
    expect_equal(cur(), "model")
    session$setInputs(grid_proceed = 1)         # "Proceed anyway"
    expect_equal(cur(), "start")
  })
})

test_that("the review step reports the candidate-set size", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none")
    expect_match(output$review_out$html,
                 "Candidate set: 9 design points")           # (2/0.25) + 1
    session$setInputs(cov_step_1 = "0.5, 0.25")
    expect_match(output$review_out$html,
                 "Candidate set: 5 design points (first step of the step sequence)",
                 fixed = TRUE)
  })
})

test_that("a comma step sequence flows through compute to the solver", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none",
              cov_step_1 = "0.5, 0.25")
    session$setInputs(compute = 1)
    cc <- computed()
    expect_null(cc$error)
    expect_equal(cc$args$step_sequence, list(0.5, 0.25))
    expect_equal(cur(), "results")
  })
})

test_that("verify panel checks the ORIGINAL criterion at the finest step", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none",
              cov_step_1 = "0.5, 0.25")
    session$setInputs(compute = 1)
    expect_null(computed()$error)

    # the args verify will use: original criterion, finest step of the sequence
    a <- owea:::.ui_solver_args(computed()$sp, "verify", computed()$theta,
                                computed()$p, computed()$subset,
                                computed()$existing)
    expect_equal(a$p, 0L)                        # the criterion chosen at compute
    expect_equal(a$step, 0.25)                   # finest step of "0.5, 0.25"

    session$setInputs(verify_manual = "dose,weight\n-1,0.5\n1,0.5",
                      verify_btn = 1)
    v <- rv$verify$v
    expect_false(inherits(v, "owea_bad"))
    expect_true(is.finite(v$max_sensitivity))    # small grid: full check ran
    expect_true(v$is_optimal$value)              # +/-1 at 1/2 IS D-optimal

    # the modal's "criterion only" option: same criterion, no optimality claim
    session$setInputs(verify_crit_only = 1)
    v2 <- rv$verify$v
    expect_identical(v2$max_sensitivity, NA_real_)
    expect_identical(v2$is_optimal$value, NA)
    expect_equal(v2$criterion, v$criterion)
  })
})

test_that("drawing theta from N(0,1) keeps the ordinal thresholds increasing", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "cumulative", ncat = 4, start = "none")
    session$setInputs(theta_draw = 1)
    th <- rv$theta_prefill             # updateNumericInput is not
    expect_length(th, 4L)                          # echoed back into input
    expect_true(all(diff(th[1:3]) > 0))
  })
})
