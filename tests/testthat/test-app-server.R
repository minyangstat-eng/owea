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

# one continuous covariate, so the model step is satisfied by the defaults.
# mode_choice picks the classical branch: the wizard now opens on the "mode"
# screen, and every test below is about the classical chain that follows it.
set_model <- function(session, link = "identity", cov_step_1 = "0.25", ...) {
  session$setInputs(mode_choice = "classical",
                    link = link, ncov = 1, cov_name_1 = "dose",
                    cov_type_1 = "continuous", cov_lo_1 = -1, cov_hi_1 = 1,
                    cov_step_1 = cov_step_1, crit = "0", qoi = "all", ...)
}

test_that("the wizard skips the theta step for the identity link", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none")
    expect_false("theta" %in% steps())
    expect_equal(cur(), "mode")            # the wizard opens on the mode screen
    goto("model")                          # ... and this test is about what follows
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

test_that("the simulation table gains a ratio column per compared design", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none",
              design_type = "exact", n_new = 12, seed = 1)
    session$setInputs(compute = 1)
    expect_null(computed()$error)

    # nothing run yet: no matrices, so the table stays empty
    session$setInputs(sim_open = 1, sim_theta_1 = 1, sim_theta_2 = 2,
                      sim_sigma = 1, sim_nsim = 20, sim_seed = 1)
    expect_null(sim_mats()[["Mean squared error"]])

    # the computed design alone: nothing to compare it against, so no ratio
    session$setInputs(run_sim = 1)
    m <- sim_mats()
    expect_equal(rownames(m[["Mean squared error"]]), "Exact design")
    expect_equal(colnames(m[["Mean squared error"]]), computed()$coef_names)
    expect_equal(dim(m[["Median squared error"]]),
                 dim(m[["Mean squared error"]]))
    expect_false(grepl("/ Exact design", output$sim_tbl, fixed = TRUE))

    # a compared design becomes a second row, and brings a ratio column with it
    session$setInputs(run_srs = 1)
    m2 <- sim_mats()
    expect_equal(rownames(m2[["Mean squared error"]]), c("Exact design", "SRS"))
    expect_equal(rownames(m2[["Median squared error"]]),
                 c("Exact design", "SRS"))

    tbl <- output$sim_tbl
    expect_true(all(vapply(
      c("Exact design (MSE)", "SRS (MSE)", "SRS / Exact design (MSE)",
        "Exact design (median SE)", "SRS (median SE)",
        "SRS / Exact design (median SE)"),
      function(s) grepl(s, tbl, fixed = TRUE), logical(1))))
    # the reference design never gets a ratio against itself
    expect_false(grepl("Exact design / Exact design", tbl, fixed = TRUE))

    # and the printed ratio really is SRS / Exact design, per parameter
    R <- owea:::.ui_sim_ratio(m2[["Mean squared error"]], "Exact design")
    expect_equal(as.numeric(R["SRS", ]),
                 as.numeric(m2[["Mean squared error"]]["SRS", ] /
                            m2[["Mean squared error"]]["Exact design", ]))
    # the note explains which way round it reads
    expect_match(output$sim_note$html, "ratios", fixed = TRUE)
  })
})

test_that("an exact design's sample size defaults to 200", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none", design_type = "exact")
    expect_match(output$n_ui$html, 'value="200"')

    # the same box means the new stage's size for an approximate design, which
    # keeps its own (much smaller) default
    session$setInputs(design_type = "approx", start = "design",
                      exist_text = "dose,count\n-1,15\n1,15", exist_n0 = 30)
    expect_match(output$n_ui$html, 'value="20"')

    # a value the user typed survives a flip between the two design types
    session$setInputs(n_new = 137, design_type = "exact")
    expect_match(output$n_ui$html, 'value="137"')
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
    goto("model")                               # past the mode screen
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

test_that("the continuous search flows through the app to the solver", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none", search_mode = "continuous",
              n_audit = 5000)
    expect_true(is_continuous())
    expect_true(spec()$continuous)
    expect_null(spec()$step_sequence)
    expect_equal(spec()$n_audit, 5000L)
    expect_match(output$review_out$html, "continuous region", fixed = TRUE)
    expect_match(output$review_out$html, "5,000 points", fixed = TRUE)
    set.seed(1)
    session$setInputs(compute = 1)
    cc <- computed()
    expect_null(cc$error)
    expect_true(cc$args$continuous)
    expect_equal(cc$args$n_audit, 5000L)
    expect_null(cc$args$step_sequence)
    expect_identical(cc$res$method, "continuous")
    expect_equal(cc$res$n_audit, 5000L)
    expect_true(cc$res$converged)
    expect_match(output$status$html, "continuous search", fixed = TRUE)
    expect_match(output$status$html, "5,000 points", fixed = TRUE)
    # the R code behind the result reproduces the call
    code <- output$code_txt
    expect_match(code, "res <- optimal_design(", fixed = TRUE)
    expect_match(code, "continuous = TRUE", fixed = TRUE)
    expect_match(code, "n_audit    = 5000", fixed = TRUE)
    expect_match(code, "verify_optimality(", fixed = TRUE)
    expect_error(parse(text = code), NA)

    # the verify panel checks over a grid at the typed step -- the grid check of
    # a design from the continuous search
    csv <- paste(c("dose,weight",
                   paste(format(cc$res$support[, 1], digits = 15),
                         format(cc$res$weights, digits = 15), sep = ",")),
                 collapse = "\n")
    session$setInputs(verify_manual = csv, verify_step = "0.02", verify_btn = 1)
    v2 <- rv$verify$v
    expect_false(inherits(v2, "owea_bad"))
    expect_identical(v2$method, "grid")
    expect_true(v2$is_optimal$value)
    expect_match(rv$verify$space, "step(s) 0.02 (101 design points)", fixed = TRUE)
    # the R code (rendered last) now ends with that same grid check
    code2 <- output$code_txt
    expect_match(code2, "res <- optimal_design(", fixed = TRUE)
    expect_match(code2, "continuous = TRUE", fixed = TRUE)   # the computation
    expect_match(code2, "verify_optimality(", fixed = TRUE)
    expect_match(code2, "step\\s+= 0\\.02")                  # the check, on the grid
    expect_error(parse(text = code2), NA)
    # an unusable step is refused before anything runs; an empty box falls back
    # to the suggested step
    session$setInputs(verify_step = "abc", verify_btn = 2)
    expect_true(inherits(rv$verify$v, "owea_bad"))
    expect_match(rv$verify$v$msg, "positive grid step")
    session$setInputs(verify_step = "", verify_btn = 3)
    expect_false(inherits(rv$verify$v, "owea_bad"))
    expect_identical(rv$verify$v$method, "grid")
    expect_true(rv$verify$v$is_optimal$value)
    # the verify panel uses the same continuous space
    a <- owea:::.ui_solver_args(cc$sp, "verify", cc$theta, cc$p, cc$subset, cc$existing)
    expect_true(a$continuous)

    # back to the grid: the step box is used again
    session$setInputs(search_mode = "grid")
    expect_false(is_continuous())
    expect_equal(spec()$step_sequence, list(0.25))
  })
})

test_that("the simulation study's simple random sample works after a continuous search", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env()
  shiny::testServer(e$server, {
    set_model(session, link = "identity", start = "none", search_mode = "continuous",
              n_audit = 0, design_type = "exact", n_new = 12, seed = 1)
    session$setInputs(compute = 1)
    expect_null(computed()$error)
    expect_identical(computed()$res$approx$method, "continuous")
    session$setInputs(sim_open = 1, sim_theta_1 = 1, sim_theta_2 = 2,
                      sim_sigma = 1, sim_nsim = 20, sim_seed = 1)
    session$setInputs(run_sim = 1)
    expect_false(inherits(rv$sim$exact, "error"))
    session$setInputs(run_srs = 1)
    expect_false(inherits(rv$sim$srs, "error"))        # used to fail: no grid step
    expect_equal(sum(rv$sim$srs_design$counts), 12L)
    expect_true(all(rv$sim$srs_design$support[, 1] >= -1 &
                    rv$sim$srs_design$support[, 1] <= 1))
    # the R code now includes the simulation study and the SRS comparison
    code <- output$code_txt
    expect_match(code, "res <- exact_design(", fixed = TRUE)
    expect_match(code, "sim <- simulate_design(", fixed = TRUE)
    expect_match(code, "sim_srs <- simulate_design(", fixed = TRUE)
    expect_match(code, "SRS = sim_srs$mse", fixed = TRUE)
    expect_error(parse(text = code), NA)
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
