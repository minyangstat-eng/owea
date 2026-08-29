# Tests for the app's COMPOUND branch: the .uic_* translators (pure, no Shiny)
# and the compound wizard's server logic.  The classical branch is tested in
# test-app.R / test-app-server.R and is untouched by anything here.

# ---- fixtures --------------------------------------------------------------
# two shared continuous covariates, coarse grid so the tests stay quick
covs2 <- function(steps = "0.5")
  list(list(name = "x1", type = "continuous", lo = -2, hi = 2,
            steps = owea:::.ui_parse_steps(steps)),
       list(name = "x2", type = "continuous", lo = -2, hi = 2,
            steps = owea:::.ui_parse_steps(steps)))

# three genuinely different models over those covariates
comps3 <- function()
  list(list(name = "m1", link = "logit", theta = c(0.5, 1, -1), p = 0),
       list(name = "m2", link = "logit", interactions = list(c(1, 2)),
            theta = c(0.5, 1, -1, 0.5), p = 0),
       list(name = "m3", link = "loglinear", theta = c(0.2, 0.3, -0.3),
            p = 1, subset = c(2, 3)))

# ---- the step graph --------------------------------------------------------

test_that(".uic_wizard_steps is the compound chain", {
  s <- owea:::.uic_wizard_steps()
  expect_equal(s, c("cmp_cov", "cmp_models", "cmp_options", "cmp_review",
                    "cmp_results"))
  # and the classical chain is NOT disturbed by the compound branch existing
  expect_equal(owea:::.ui_wizard_steps("identity", "none"),
               c("model", "start", "criterion", "design_type", "review",
                 "results"))
})

# ---- specs -----------------------------------------------------------------

test_that(".uic_specs builds one spec per component over shared covariates", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  expect_length(sp, 3L)
  # every component sees the SAME design region
  for (s in sp) expect_equal(names(s$design_box), c("x1", "x2"))
  # but its own model: 3, 4 and 3 parameters
  expect_equal(vapply(sp, function(s) length(owea:::.ui_coef_names(s)),
                      integer(1)), c(3L, 4L, 3L))
  expect_equal(sp[[1]]$link, "logit")
  expect_equal(sp[[3]]$link, "loglinear")
  expect_null(sp[[1]]$xx)                       # no interaction on component 1
  expect_equal(sp[[2]]$xx, list(c(1L, 2L)))     # x1 x2 on component 2
})

test_that(".uic_specs rejects an empty component list", {
  expect_error(owea:::.uic_specs(covs2(), list()), "at least one component")
})

# ---- validation ------------------------------------------------------------

test_that(".uic_check accepts a good setup and names the bad component", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  expect_null(owea:::.uic_check(sp, comps3(), c(0.4, 0.35, 0.25)))

  bad <- comps3(); bad[[2]]$theta <- c(1, 2)            # needs 4, given 2
  expect_match(owea:::.uic_check(sp, bad, c(1, 1, 1)), "^m2: ")
  expect_match(owea:::.uic_check(sp, bad, c(1, 1, 1)), "4 parameter")

  bad2 <- comps3(); bad2[[3]]$subset <- c(2, 9)         # out of range for k = 3
  expect_match(owea:::.uic_check(sp, bad2, c(1, 1, 1)), "m3.*1\\.\\.3")

  expect_match(owea:::.uic_check(sp, comps3(), c(1, 1)), "one weight per")
  expect_match(owea:::.uic_check(sp, comps3(), c(-1, 1, 1)), "nonnegative")
  expect_match(owea:::.uic_check(sp, comps3(), c(0, 0, 0)), "cannot all be zero")
})

# ---- the argument assembler ------------------------------------------------

test_that(".uic_solver_args shapes the design and verify targets", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  a  <- owea:::.uic_solver_args(sp, comps3(), "design", alpha = c(.4, .35, .25))

  expect_equal(names(a$design_box), c("x1", "x2"))
  expect_equal(a$step_sequence, sp[[1]]$step_sequence)
  expect_true(a$efficiency)
  expect_length(a$components, 3L)
  # the design region is a TOP-LEVEL argument, never per component
  expect_false(any(vapply(a$components,
                          function(z) "design_box" %in% names(z), logical(1))))
  # each component carries its own model and criterion
  expect_equal(a$components[[2]]$xx, list(c(1L, 2L)))
  expect_equal(a$components[[3]]$p, 1L)
  expect_equal(a$components[[3]]$subset, c(2L, 3L))
  expect_null(a$components[[1]]$subset)          # "all" stays absent
  expect_equal(a$components[[1]]$coding, "zero-sum")

  # verify wants a single step, not a sequence, and reuses psi_star
  v <- owea:::.uic_solver_args(sp, comps3(), "verify", alpha = c(.4, .35, .25),
                               psi_star = c(1, 2, 3))
  expect_equal(v$step, sp[[1]]$finest)
  expect_null(v$step_sequence)
  expect_equal(v$psi_star, c(1, 2, 3))
})

test_that(".uic_solver_args omits xi0_* when there is no existing design", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  a  <- owea:::.uic_solver_args(sp, comps3(), "design")
  expect_false(any(c("xi0_points", "xi0_weights", "n0", "n1") %in% names(a)))
  # ... and passes them through when there is one
  ex <- list(points = rbind(c(-2, -2), c(2, 2)), weights = c(0.5, 0.5),
             n0 = 10, n1 = 20)
  b <- owea:::.uic_solver_args(sp, comps3(), "design", existing = ex)
  expect_equal(b$n0, 10L)
  expect_equal(b$n1, 20L)
  expect_equal(sum(b$xi0_weights), 1)
  expect_equal(nrow(b$xi0_points), 2L)
})

test_that(".uic_solver_args demands assumed values for a local model", {
  sp  <- owea:::.uic_specs(covs2(), comps3())
  bad <- comps3(); bad[[1]]$theta <- NULL
  expect_error(owea:::.uic_solver_args(sp, bad, "design"),
               "assumed parameter value")
  # ... but not for the identity link
  ci  <- list(list(name = "lin", link = "identity", p = 0))
  spi <- owea:::.uic_specs(covs2(), ci)
  ai  <- owea:::.uic_solver_args(spi, ci, "design")
  expect_null(ai$components[[1]]$theta)
})

# ---- parity with a direct compound_design() call ---------------------------

test_that(".uic_solver_args reproduces a direct compound_design() call", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  a  <- owea:::.uic_solver_args(sp, comps3(), "design", alpha = c(.4, .35, .25))
  r  <- do.call(compound_design, a)

  d <- compound_design(
    components = list(
      list(link = "logit", x = c(1, 2), theta = c(0.5, 1, -1), p = 0),
      list(link = "logit", x = c(1, 2), xx = list(c(1, 2)),
           theta = c(0.5, 1, -1, 0.5), p = 0),
      list(link = "loglinear", x = c(1, 2), theta = c(0.2, 0.3, -0.3),
           p = 1, subset = c(2, 3))),
    alpha = c(.4, .35, .25),
    design_box = list(x1 = c(-2, 2), x2 = c(-2, 2)),
    step_sequence = list(c(0.5, 0.5)))
  expect_equal(r$criterion, d$criterion, tolerance = 1e-8)
  expect_true(r$converged)
})

# ---- results helpers -------------------------------------------------------

test_that("the results tables have the shape the app renders", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  r  <- do.call(compound_design,
                owea:::.uic_solver_args(sp, comps3(), "design",
                                        alpha = c(.4, .35, .25)))
  s <- owea:::.uic_summary_table(r)
  expect_equal(rownames(s), c("m1", "m2", "m3"))
  expect_true(all(c("alpha", "criterion", "Psi", "Psi_star", "efficiency")
                  %in% names(s)))
  expect_equal(s$criterion, c("D", "D", "A"))

  ct <- owea:::.uic_cross_table(r)
  expect_equal(dim(ct), c(4L, 3L))
  expect_equal(rownames(ct)[4], "THIS DESIGN")

  # raw weighting runs no reference solves, so there is no cross table
  rr <- do.call(compound_design,
                owea:::.uic_solver_args(sp, comps3(), "design",
                                        alpha = c(.4, .35, .25),
                                        efficiency = FALSE))
  expect_null(owea:::.uic_cross_table(rr))
  expect_false("efficiency" %in% names(owea:::.uic_summary_table(rr)))
})

test_that(".uic_verify_points counts the shared grid at the finest step", {
  sp <- owea:::.uic_specs(covs2("0.5, 0.25"), comps3())
  expect_equal(owea:::.uic_verify_points(sp), owea:::.ui_verify_points(sp[[1]]))
  expect_equal(owea:::.uic_verify_points(sp), 17 * 17)     # (4/0.25 + 1)^2
})

# ---- the server ------------------------------------------------------------

app_env_c <- function() {
  dir <- system.file("shiny", "owea-app", package = "owea")
  testthat::skip_if(!nzchar(dir) || !file.exists(file.path(dir, "app.R")),
                    "the bundled app is not on disk")
  e <- new.env(parent = globalenv())
  sys.source(file.path(dir, "app.R"), envir = e)
  e
}

# two shared covariates and two objectives, on a coarse grid.  Anything passed
# through ... REPLACES the corresponding default rather than being appended, so
# a test can override one input without a duplicate-argument error.
set_compound <- function(session, ...) {
  args <- list(
    mode_choice = "compound",
    cmp_ncov = 2,
    ccov_name_1 = "x1", ccov_type_1 = "continuous",
    ccov_lo_1 = -2, ccov_hi_1 = 2, ccov_step_1 = "0.5",
    ccov_name_2 = "x2", ccov_type_2 = "continuous",
    ccov_lo_2 = -2, ccov_hi_2 = 2, ccov_step_2 = "0.5",
    ncomp = 2,
    cmp_name_1 = "main", cmp_link_1 = "logit", cmp_alpha_1 = 0.5,
    cmp_crit_1 = "0", cmp_qoi_1 = "all",
    cmp_theta_1_1 = 0.5, cmp_theta_1_2 = 1, cmp_theta_1_3 = -1,
    cmp_name_2 = "inter", cmp_link_2 = "logit", cmp_alpha_2 = 0.5,
    cmp_crit_2 = "0", cmp_qoi_2 = "all", cmp_int_2 = "1-2",
    cmp_theta_2_1 = 0.5, cmp_theta_2_2 = 1, cmp_theta_2_3 = -1,
    cmp_theta_2_4 = 0.5,
    cmp_efficiency = "eff", cmp_start = "none")
  over <- list(...)
  args[names(over)] <- over
  do.call(session$setInputs, args)
}

test_that("choosing compound switches the wizard to the compound chain", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    session$setInputs(mode_choice = "classical")
    expect_true("model" %in% steps())
    expect_false("cmp_cov" %in% steps())

    session$setInputs(mode_choice = "compound")
    expect_true(is_compound())
    expect_equal(steps(), c("mode", "cmp_cov", "cmp_models", "cmp_options",
                            "cmp_review", "cmp_results"))
    expect_false("model" %in% steps())
  })
})

test_that("the compound branch reaches the solver and converges", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session)
    expect_length(cmp_specs(), 2L)
    expect_equal(vapply(cmp_specs(),
                        function(s) length(owea:::.ui_coef_names(s)), integer(1)),
                 c(3L, 4L))
    expect_null(step_error("cmp_cov"))
    expect_null(step_error("cmp_models"))
    expect_null(step_error("cmp_options"))

    session$setInputs(cmp_compute = 1)
    cc <- cmp_computed()
    expect_null(cc$error)
    expect_true(cc$res$converged)
    expect_lt(cc$res$max_d, 1e-8)
    expect_equal(cur(), "cmp_results")             # a good compute advances

    # the efficiency-weighted value is an average efficiency
    expect_true(cc$res$criterion > 0 && cc$res$criterion <= 1 + 1e-10)
    expect_true(all(cc$res$efficiency <= 1 + 1e-10))
    # ... and the cross-efficiency table is there for the results tab
    expect_equal(dim(cc$res$cross_efficiency), c(3L, 2L))
    # xi0 args are absent with no existing design
    expect_false(any(c("xi0_points", "n0", "n1") %in% names(cc$args)))
  })
})

test_that("the compound branch passes an existing design through", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session, cmp_start = "design",
                 cmp_exist_text = "x1,x2,count\n-2,-2,5\n2,2,5",
                 cmp_exist_n0 = 10, cmp_exist_n1 = 30)
    ex <- cmp_existing()
    expect_equal(nrow(ex$points), 2L)
    expect_equal(as.integer(ex$n0), 10L)
    expect_equal(as.integer(ex$n1), 30L)

    session$setInputs(cmp_compute = 1)
    cc <- cmp_computed()
    expect_null(cc$error)
    expect_equal(cc$args$n0, 10L)
    expect_equal(cc$args$n1, 30L)
    expect_equal(sum(cc$args$xi0_weights), 1)
    # the reference solves used the same stage structure, so efficiencies
    # cannot exceed 1
    expect_true(all(cc$res$efficiency <= 1 + 1e-10))
  })
})

test_that("a bad component blocks Next, naming the component", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    # a subset that reaches past the component's parameter count
    set_compound(session, cmp_qoi_2 = "subset", cmp_subset_2 = c("1", "9"))
    err <- step_error("cmp_models")
    expect_false(is.null(err))
    expect_match(err, "inter")                    # says WHICH objective
    expect_match(err, "1\\.\\.4")                 # and what the valid range is
  })

  # a negative weight is rejected too
  e2 <- app_env_c()
  shiny::testServer(e2$server, {
    set_compound(session, cmp_alpha_1 = -1)
    expect_match(step_error("cmp_models"), "nonnegative")
  })
})

test_that("an empty assumed-value box reads as zero, not as an error", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    # the app's %||% treats NA as missing, so a blank box falls back to 0 --
    # a valid assumed value, not a validation failure
    set_compound(session, cmp_theta_2_4 = NA)
    expect_equal(cmp_comps()[[2]]$theta, c(0.5, 1, -1, 0))
    expect_null(step_error("cmp_models"))
  })
})

test_that("raw weighting drops the reference solves", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session, cmp_efficiency = "raw")
    session$setInputs(cmp_compute = 1)
    cc <- cmp_computed()
    expect_null(cc$error)
    expect_false(cc$args$efficiency)
    expect_true(all(is.na(cc$res$psi_star)))
    expect_null(cc$res$cross_efficiency)
  })
})

# ---- exact designs and the simulation study --------------------------------

test_that(".uic_solver_args shapes the exact target", {
  sp <- owea:::.uic_specs(covs2(), comps3())
  a  <- owea:::.uic_solver_args(sp, comps3(), "exact", alpha = c(.4, .35, .25),
                                n = 25)
  expect_equal(a$n, 25L)
  expect_equal(a$step_sequence, sp[[1]]$step_sequence)   # exact still refines
  expect_error(owea:::.uic_solver_args(sp, comps3(), "exact"), "number of runs")
})

test_that(".uic_simulate runs each component under its own model", {
  sp  <- owea:::.uic_specs(covs2(), comps3())
  sup <- rbind(c(-2, -2), c(2, -2), c(-2, 2), c(2, 2), c(0, 0))
  cnt <- c(8L, 8L, 8L, 8L, 8L)
  s <- owea:::.uic_simulate(sp, comps3(), support = sup, counts = cnt,
                            nsim = 30L, seed = 1)
  expect_length(s, 3L)
  for (j in seq_len(3)) {
    expect_false(inherits(s[[j]], "error"))
    # each component is fitted with ITS OWN model, so its own parameter count
    expect_length(s[[j]]$mse, length(owea:::.ui_coef_names(sp[[j]])))
    expect_equal(s[[j]]$N, sum(cnt))
    expect_gt(s[[j]]$n_converged, 0L)
  }
  expect_equal(vapply(s, function(z) length(z$mse), integer(1)), c(3L, 4L, 3L))
  m <- owea:::.uic_sim_mse(s[[2]])
  expect_length(m, 4L)
  expect_equal(names(m), owea:::.ui_coef_names(sp[[2]]))
  expect_true(all(m >= 0))
})

test_that(".uic_simulate reports a failing component without losing the rest", {
  sp  <- owea:::.uic_specs(covs2(), comps3())
  bad <- comps3(); bad[[1]]$theta <- c(1, 2)          # wrong length for k = 3
  sup <- rbind(c(-2, -2), c(2, 2), c(0, 0)); cnt <- c(5L, 5L, 5L)
  s <- owea:::.uic_simulate(sp, bad, support = sup, counts = cnt,
                            nsim = 10L, seed = 1)
  expect_true(inherits(s[[1]], "error"))
  expect_null(owea:::.uic_sim_mse(s[[1]]))
  expect_false(inherits(s[[2]], "error"))             # the others still ran
})

test_that("the app's exact branch reaches compound_exact_design", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session, cmp_design_type = "exact", cmp_n = 24, cmp_seed = 1)
    expect_true(cmp_is_exact())
    session$setInputs(cmp_compute = 1)
    cc <- cmp_computed()
    expect_null(cc$error)
    expect_true(cc$exact)
    expect_equal(cc$args$n, 24L)
    expect_s3_class(cc$res, "compound_exact_design")
    expect_equal(sum(cc$res$counts), 24L)
    expect_equal(cur(), "cmp_results")
    # the design table switches from weights to counts
    expect_true("count" %in% names(cmp_design_df()))
    expect_false("weight" %in% names(cmp_design_df()))
  })
})

test_that("the app's simulation study runs per component", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session, cmp_design_type = "exact", cmp_n = 24, cmp_seed = 1)
    session$setInputs(cmp_compute = 1)
    expect_null(cmp_computed()$error)
    session$setInputs(cmp_sim_open = 1,
                      cmp_sim_theta_1_1 = 0.5, cmp_sim_theta_1_2 = 1,
                      cmp_sim_theta_1_3 = -1,
                      cmp_sim_theta_2_1 = 0.5, cmp_sim_theta_2_2 = 1,
                      cmp_sim_theta_2_3 = -1, cmp_sim_theta_2_4 = 0.5,
                      cmp_sim_nsim = 30, cmp_sim_seed = 3, cmp_sim_sigma = 1)
    session$setInputs(cmp_run_sim = 1)
    s <- rv$cmp_sim$design
    expect_length(s, 2L)                       # one result per component
    expect_length(owea:::.uic_sim_mse(s[[1]]), 3L)
    expect_length(owea:::.uic_sim_mse(s[[2]]), 4L)

    session$setInputs(cmp_run_srs = 1)
    expect_length(rv$cmp_sim$srs, 2L)
    expect_equal(rv$cmp_sim$srs[[1]]$N, 24L)   # the SRS design has the same n

    # a custom design must have the same number of runs
    session$setInputs(cmp_sim_custom = "x1,x2,count\n-2,-2,5\n2,2,5",
                      cmp_run_custom = 1)
    expect_true(inherits(rv$cmp_sim$custom, "error"))
    expect_match(conditionMessage(rv$cmp_sim$custom), "sum to 24")

    # one matrix per objective: the designs compared are its rows, that
    # objective's parameters its columns
    nms <- vapply(cmp_computed()$comps, function(z) z$name, character(1))
    m1 <- cmp_sim_mats()
    expect_equal(rownames(m1[["Mean squared error"]]), c("This design", "Random"))
    expect_length(colnames(m1[["Mean squared error"]]), 3L)   # component 1
    expect_equal(dim(m1[["Median squared error"]]),
                 dim(m1[["Mean squared error"]]))

    # the compared design brings a ratio column, measured against this design
    tbl <- output$cmp_sim_tbl
    expect_true(all(vapply(
      c("This design (MSE)", "Random (MSE)", "Random / This design (MSE)",
        "Random / This design (median SE)"),
      function(s) grepl(s, tbl, fixed = TRUE), logical(1))))
    expect_false(grepl("This design / This design", tbl, fixed = TRUE))
    R <- owea:::.ui_sim_ratio(m1[["Mean squared error"]], "This design")
    expect_equal(as.numeric(R["Random", ]),
                 as.numeric(m1[["Mean squared error"]]["Random", ] /
                            m1[["Mean squared error"]]["This design", ]))

    # switching objective switches the table, ratio and all
    session$setInputs(cmp_sim_which = nms[2])
    expect_equal(cmp_sim_which()$j, 2L)
    expect_length(colnames(cmp_sim_mats()[["Mean squared error"]]), 4L)
    expect_true(grepl("Random / This design (MSE)", output$cmp_sim_tbl,
                      fixed = TRUE))
  })
})

test_that("the compound exact branch offers 200 runs by default", {
  skip_if_not_installed("shiny")
  dir <- system.file("shiny", "owea-app", package = "owea")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "app.R")),
          "the bundled app is not on disk")
  src <- paste(readLines(file.path(dir, "app.R")), collapse = "\n")
  expect_match(src, 'numericInput\\("cmp_n", "Number of runs \\(n\\)",\\s*\n?\\s*value = N_DEFAULT_EXACT')
  expect_match(src, "N_DEFAULT_EXACT <- 200")
})

# ---- reference-value (Psi*) caching ---------------------------------------

test_that(".uic_psi_key ignores what Psi* does not depend on", {
  sp <- owea:::.uic_specs(covs2(), comps3()[1:2])
  cm <- comps3()[1:2]
  a1 <- owea:::.uic_solver_args(sp, cm, "design", alpha = c(0.5, 0.5))
  a2 <- owea:::.uic_solver_args(sp, cm, "design", alpha = c(0.9, 0.1))
  # alpha does not enter a reference solve (each runs at weight 1), so the
  # stored Psi* is still valid when only the weights move
  expect_equal(owea:::.uic_psi_key(a1), owea:::.uic_psi_key(a2))

  # nor do n, the seed, or the efficiency switch
  a3 <- owea:::.uic_solver_args(sp, cm, "exact", alpha = c(0.5, 0.5), n = 30)
  a3$seed <- 7
  expect_equal(owea:::.uic_psi_key(a1), owea:::.uic_psi_key(a3))

  # but the objectives themselves do: a different theta is a different Psi*
  cm2 <- cm; cm2[[1]]$theta <- c(0.9, 1, -1)
  a4 <- owea:::.uic_solver_args(owea:::.uic_specs(covs2(), cm2), cm2, "design",
                                alpha = c(0.5, 0.5))
  expect_false(identical(owea:::.uic_psi_key(a1), owea:::.uic_psi_key(a4)))

  # and so does the design region / stage structure
  a5 <- owea:::.uic_solver_args(owea:::.uic_specs(covs2("1, 0.5"), cm), cm,
                                "design", alpha = c(0.5, 0.5))
  expect_false(identical(owea:::.uic_psi_key(a1), owea:::.uic_psi_key(a5)))
})

test_that("the app reuses Psi* when only the weights change", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session)
    session$setInputs(cmp_compute = 1)
    c1 <- cmp_computed()
    expect_null(c1$error)
    expect_false(isTRUE(c1$psi_reused))          # nothing cached on the first run
    ps1 <- c1$res$psi_star

    # same objectives, different weights: the reference solves are skipped and
    # the SAME Psi* is used, so the efficiencies stay on the same scale
    expect_false(grepl("were reused", output$cmp_status$html, fixed = TRUE))

    session$setInputs(cmp_alpha_1 = 0.9, cmp_alpha_2 = 0.1, cmp_compute = 2)
    c2 <- cmp_computed()
    expect_null(c2$error)
    expect_true(isTRUE(c2$psi_reused))
    expect_equal(as.numeric(c2$res$psi_star), as.numeric(ps1), tolerance = 1e-12)
    # and the status line says so, so a fast re-run is not a mystery
    expect_true(grepl("were reused", output$cmp_status$html, fixed = TRUE))

    # change an objective and the cache must NOT be used
    session$setInputs(cmp_theta_1_1 = 1.25, cmp_compute = 3)
    c3 <- cmp_computed()
    expect_null(c3$error)
    expect_false(isTRUE(c3$psi_reused))
  })
})

test_that("reusing Psi* gives the same design as recomputing it", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session, cmp_alpha_1 = 0.7, cmp_alpha_2 = 0.3)
    session$setInputs(cmp_compute = 1)
    fresh <- cmp_computed()$res                  # computed its own Psi*

    # same inputs again: this run reuses the cache
    session$setInputs(cmp_compute = 2)
    c2 <- cmp_computed()
    expect_true(isTRUE(c2$psi_reused))
    cached <- c2$res
    expect_equal(cached$criterion, fresh$criterion, tolerance = 1e-10)
    expect_equal(as.numeric(cached$weights), as.numeric(fresh$weights),
                 tolerance = 1e-8)
    expect_equal(as.numeric(cached$efficiency), as.numeric(fresh$efficiency),
                 tolerance = 1e-10)
  })
})

# ---- the verify grid-size gate ---------------------------------------------

# the app's own design_csv_text() is a file-level global that testServer's mask
# does not reach, so the tests below build the pasted design themselves
as_csv <- function(d)
  paste(utils::capture.output(utils::write.csv(d, row.names = FALSE)),
        collapse = "\n")

test_that(".uic_solver_args disables the package cap for the app's verify", {
  sp <- owea:::.uic_specs(covs2("0.05"), comps3())
  v  <- owea:::.uic_solver_args(sp, comps3(), "verify", psi_star = c(1, 2, 3))
  # the app gates the grid itself, with a modal; the package must not ALSO
  # stop (or, inside a Shiny session, open a console menu nobody can answer)
  expect_equal(v$max_points, Inf)
  # the design and exact targets never scan, so they carry no cap
  expect_null(owea:::.uic_solver_args(sp, comps3(), "design")$max_points)
})

test_that("the app gates a huge verify grid instead of building it", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    # a 0.002 step on [-2, 2]^2 is 2001^2 = 4,004,001 points
    set_compound(session, ccov_step_1 = "0.002", ccov_step_2 = "0.002")
    expect_gt(owea:::.uic_verify_points(cmp_specs()), 1e6)
  })

  # a coarse grid is under the cap and scores straight through
  shiny::testServer(e$server, {
    set_compound(session)
    session$setInputs(cmp_compute = 1)
    expect_null(cmp_computed()$error)
    expect_lt(owea:::.uic_verify_points(cmp_specs()), 1e6)
    session$setInputs(cmp_verify_manual = as_csv(cmp_design_df()),
                      cmp_verify_btn = 1)
    v <- rv$cmp_verify
    expect_false(inherits(v, "error"))
    expect_true(is.finite(v$max_d))              # the full scan ran
  })
})

test_that("the criterion-only route scores without the sensitivity scan", {
  skip_if_not_installed("shiny"); skip_if_not_installed("DT")
  e <- app_env_c()
  shiny::testServer(e$server, {
    set_compound(session)
    session$setInputs(cmp_compute = 1)
    expect_null(cmp_computed()$error)
    session$setInputs(cmp_verify_manual = as_csv(cmp_design_df()))

    # the modal's "criterion only" button
    session$setInputs(cmp_verify_crit_only = 1)
    v <- rv$cmp_verify
    expect_false(inherits(v, "error"))
    expect_true(is.finite(v$criterion))
    expect_null(v$max_d)                          # no grid was built
    # and the panel says optimality was not assessed rather than erroring
    expect_match(output$cmp_verify_out$html, "Optimality was NOT assessed")

    # the modal's "proceed anyway" button reports the sensitivity, and the
    # SAME criterion -- skipping the scan changes nothing about the score
    session$setInputs(cmp_verify_full = 1)
    v2 <- rv$cmp_verify
    expect_true(is.finite(v2$max_d))
    expect_equal(v2$criterion, v$criterion, tolerance = 1e-10)
    expect_false(grepl("NOT assessed", output$cmp_verify_out$html, fixed = TRUE))
  })
})
