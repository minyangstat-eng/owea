# ===========================================================================
# owea web app -- point-and-click optimal designs.
#
# Launch locally with owea::run_owea_app(), or deploy this folder to
# shinyapps.io / Posit Connect (see deploy/deploy.R in the source repo).
#
# The app depends only on the exported owea API plus the internal helpers
# .ui_model_spec() / .ui_coef_names() / .plot_design() (see R/app.R, R/plot.R),
# which translate the friendly inputs into optimal_design()/exact_design() args.
# ===========================================================================

library(shiny)
library(owea)

MAX_COV   <- 6L
LINKS     <- c("Linear (normal)"        = "identity",
               "Logistic (binary)"      = "logit",
               "Poisson (counts)"       = "loglinear",
               "Multinomial (nominal)"  = "multinomial",
               "Ordinal (proportional odds)" = "cumulative")
MULTI_CAT <- c("multinomial", "cumulative")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("owea — Optimal Experimental Designs"),
  tags$p("Describe your model, click ", tags$b("Compute design"),
         ". No R code required."),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput("link", "Model family", choices = LINKS),
      conditionalPanel(
        "input.link == 'multinomial' || input.link == 'cumulative'",
        numericInput("ncat", "Number of response categories", value = 3,
                     min = 2, step = 1)),
      hr(),
      numericInput("ncov", "Number of covariates", value = 1, min = 1,
                   max = MAX_COV, step = 1),
      uiOutput("cov_ui"),
      hr(),
      uiOutput("interaction_ui"),
      uiOutput("quadratic_ui"),
      hr(),
      selectInput("crit", "Optimality criterion",
                  choices = c("D-optimal (overall precision)" = "0",
                              "A-optimal (average variance)"   = "1")),
      radioButtons("qoi", "Parameters of interest",
                   choices = c("All parameters" = "all",
                               "A subset"        = "subset")),
      conditionalPanel("input.qoi == 'subset'", uiOutput("subset_ui")),
      hr(),
      conditionalPanel("input.link != 'identity'",
                       tags$b("Assumed parameter values (theta)"),
                       helpText("Local designs depend on these guesses."),
                       uiOutput("theta_ui")),
      hr(),
      radioButtons("design_type", "Design type",
                   choices = c("Approximate (weights)" = "approx",
                               "Exact (integer runs)"   = "exact")),
      conditionalPanel(
        "input.design_type == 'exact'",
        numericInput("n_runs", "Total runs n", value = 20, min = 1, step = 1),
        numericInput("seed", "Random seed", value = 1, step = 1)),
      br(),
      actionButton("compute", "Compute design", class = "btn-primary",
                   width = "100%")
    ),
    mainPanel(
      width = 8,
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Design",
                 br(), DT::DTOutput("design_tbl"),
                 br(), downloadButton("dl_csv", "Download design (CSV)")),
        tabPanel("Plot",
                 br(), plotOutput("design_plot", height = "460px"),
                 br(), downloadButton("dl_png", "Download plot (PNG)")),
        tabPanel("Information matrix",
                 br(), helpText("The resulting Fisher information matrix ",
                                "(combined with any existing design)."),
                 tableOutput("info_tbl")),
        tabPanel("Model",
                 br(), verbatimTextOutput("model_txt"))
      ),
      uiOutput("post_result")
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- dynamic covariate blocks ------------------------------------------
  output$cov_ui <- renderUI({
    n <- max(1L, min(MAX_COV, as.integer(input$ncov %||% 1)))
    lapply(seq_len(n), function(i) {
      wellPanel(
        style = "padding:8px;",
        textInput(paste0("cov_name_", i), NULL, value = paste0("x", i),
                  placeholder = "covariate name"),
        selectInput(paste0("cov_type_", i), NULL,
                    choices = c("Continuous" = "continuous", "Factor" = "factor")),
        conditionalPanel(
          sprintf("input.cov_type_%d == 'continuous'", i),
          fluidRow(
            column(4, numericInput(paste0("cov_lo_", i), "low",  value = -1)),
            column(4, numericInput(paste0("cov_hi_", i), "high", value = 1)),
            column(4, numericInput(paste0("cov_step_", i), "grid step",
                                   value = 0.1, min = 1e-6)))),
        conditionalPanel(
          sprintf("input.cov_type_%d == 'factor'", i),
          numericInput(paste0("cov_nlev_", i), "number of levels", value = 2,
                       min = 2, step = 1))
      )
    })
  })

  # ---- gather covariates into the translator's input shape ----------------
  covariates <- reactive({
    n <- max(1L, min(MAX_COV, as.integer(input$ncov %||% 1)))
    out <- lapply(seq_len(n), function(i) {
      type <- input[[paste0("cov_type_", i)]] %||% "continuous"
      nm   <- input[[paste0("cov_name_", i)]] %||% paste0("x", i)
      if (identical(type, "factor")) {
        list(name = nm, type = "factor",
             nlevels = as.integer(input[[paste0("cov_nlev_", i)]] %||% 2))
      } else {
        list(name = nm, type = "continuous",
             lo = as.numeric(input[[paste0("cov_lo_", i)]]   %||% -1),
             hi = as.numeric(input[[paste0("cov_hi_", i)]]   %||%  1),
             step = as.numeric(input[[paste0("cov_step_", i)]] %||% 0.1))
      }
    })
    out
  })

  # ---- interaction & quadratic choices, built from current covariates -----
  output$interaction_ui <- renderUI({
    cv <- covariates(); n <- length(cv)
    if (n < 2) return(NULL)
    pairs <- combn(n, 2, simplify = FALSE)
    choices <- setNames(
      vapply(pairs, function(p) paste(p, collapse = "-"), character(1)),
      vapply(pairs, function(p) sprintf("%s × %s", cv[[p[1]]]$name, cv[[p[2]]]$name),
             character(1)))
    checkboxGroupInput("interactions", "Interactions", choices = choices)
  })
  output$quadratic_ui <- renderUI({
    cv <- covariates()
    cont <- which(vapply(cv, function(c) identical(c$type, "continuous"), logical(1)))
    if (!length(cont)) return(NULL)
    choices <- setNames(as.character(cont),
                        vapply(cont, function(i) cv[[i]]$name, character(1)))
    checkboxGroupInput("quadratics", "Quadratic terms", choices = choices)
  })

  # ---- the model spec (validated) ----------------------------------------
  spec <- reactive({
    inter <- lapply(input$interactions %||% character(0),
                    function(s) as.integer(strsplit(s, "-", fixed = TRUE)[[1]]))
    quad  <- as.integer(input$quadratics %||% character(0))
    ncat  <- if (input$link %in% MULTI_CAT) as.integer(input$ncat %||% 3) else NULL
    owea:::.ui_model_spec(covariates(), interactions = inter, quadratics = quad,
                          link = input$link, ncat = ncat)
  })

  coef_names <- reactive({
    tryCatch(owea:::.ui_coef_names(spec()), error = function(e) NULL)
  })

  # ---- dynamic theta and subset inputs -----------------------------------
  output$theta_ui <- renderUI({
    cn <- coef_names(); if (is.null(cn)) return(helpText("(define the model first)"))
    lapply(seq_along(cn), function(i)
      numericInput(paste0("theta_", i), cn[i], value = 0, step = 0.1))
  })
  output$subset_ui <- renderUI({
    cn <- coef_names(); if (is.null(cn)) return(NULL)
    checkboxGroupInput("subset_terms", "Include parameters", choices = cn)
  })

  # ---- compute on click, gated by a large-grid (>1e6) confirmation --------
  go    <- reactiveVal(0)
  gridN <- reactive(tryCatch(max(owea:::.ui_grid_sizes(covariates())),
                             error = function(e) 0))
  observeEvent(input$compute, {
    N <- gridN()
    if (is.finite(N) && N > 1e6)
      showModal(modalDialog(title = "Large candidate set",
        sprintf(paste0("The design space has %s design points at this step size; ",
          "building it may be slow."), format(round(N), big.mark = ",")),
        "Proceed anyway, or cancel and increase the step size?",
        footer = tagList(modalButton("Cancel — change step"),
                         actionButton("big_proceed", "Proceed anyway", class = "btn-warning"))))
    else go(go() + 1)
  })
  observeEvent(input$big_proceed, { removeModal(); go(go() + 1) })

  computed <- eventReactive(go(), {
    sp <- tryCatch(spec(), error = function(e)
      structure(list(msg = conditionMessage(e)), class = "owea_bad"))
    if (inherits(sp, "owea_bad")) return(list(error = sp$msg))

    cn <- owea:::.ui_coef_names(sp)
    theta <- NULL
    if (!identical(sp$link, "identity")) {
      theta <- vapply(seq_along(cn), function(i)
        as.numeric(input[[paste0("theta_", i)]] %||% 0), numeric(1))
    }
    p <- as.integer(input$crit)
    subset <- NULL
    if (identical(input$qoi, "subset")) {
      sel <- input$subset_terms %||% character(0)
      if (length(sel)) subset <- match(sel, cn)
    }
    args <- list(design_box = sp$design_box, step_sequence = sp$step_sequence,
                 link = sp$link, ncat = sp$ncat, f = sp$f, x = sp$x,
                 fx = sp$fx, xx = sp$xx, ff = sp$ff,
                 theta = theta, p = p, subset = subset)

    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch(
        if (identical(input$design_type, "exact"))
          do.call(exact_design,
                  c(list(n = as.integer(input$n_runs %||% 1),
                         seed = as.integer(input$seed %||% 1)), args))
        else
          do.call(optimal_design, args),
        error = function(e) structure(list(msg = conditionMessage(e)),
                                      class = "owea_bad")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w))
                              invokeRestart("muffleWarning") })
    if (inherits(res, "owea_bad")) return(list(error = res$msg))
    list(res = res, warns = warns, cov_names = names(sp$design_box),
         exact = identical(input$design_type, "exact"),
         sp = sp, theta = theta, p = p, subset = subset, finest = sp$finest,
         coef_names = cn, n = if (identical(input$design_type, "exact"))
                                as.integer(input$n_runs %||% 1) else NA_integer_)
  }, ignoreInit = TRUE)

  # ---- status banner ------------------------------------------------------
  output$status <- renderUI({
    c <- computed()
    if (!is.null(c$error))
      return(div(class = "alert alert-danger",
                 tags$b("Could not compute: "), c$error))
    res <- c$res
    conv <- isTRUE(res$converged)
    msg <- if (conv) "Design computed and verified optimal on the searched grid."
           else "Computed, but optimality was not fully certified — try a finer grid step."
    extra <- if (length(c$warns))
      tags$ul(lapply(unique(.friendly_warn(c$warns)), tags$li)) else NULL
    div(class = if (conv) "alert alert-success" else "alert alert-warning",
        tags$b(msg), extra)
  })

  # ---- design table -------------------------------------------------------
  # full precision: this frame feeds the download CSV and the Verify prefill.
  design_df <- reactive({
    c <- computed(); req(is.null(c$error))
    res <- c$res
    S <- as.matrix(res$support)
    colnames(S) <- c$cov_names
    df <- as.data.frame(S)
    if (c$exact) df[["count"]] <- res$counts
    else         df[["weight"]] <- res$weights
    df
  })
  output$design_tbl <- DT::renderDT({
    c <- computed(); req(is.null(c$error))
    res <- c$res
    cap <- if (c$exact)
      sprintf("n = %d | criterion (per-sample) = %.5f | criterion (total) = %.5f | efficiency >= %.1f%%",
              res$n, res$criterion, res$criterion_total, 100 * res$efficiency)
    else
      sprintf("criterion = %.5f | max sensitivity = %.2e (0 at optimum)",
              res$criterion, res$max_d)
    disp <- design_df()                              # round a display copy only
    for (nm in names(disp))
      if (is.numeric(disp[[nm]]) && !identical(nm, "count"))
        disp[[nm]] <- round(disp[[nm]], 4)
    DT::datatable(disp, rownames = FALSE, caption = cap,
                  options = list(dom = "t", paging = FALSE))
  })

  # ---- plot ---------------------------------------------------------------
  output$design_plot <- renderPlot({
    c <- computed(); req(is.null(c$error))
    plot_design(c$res, cov_names = c$cov_names,
                main = if (c$exact) "Exact design" else "Approximate design")
  })

  # ---- information matrix -------------------------------------------------
  output$info_tbl <- renderTable({
    c <- computed(); req(is.null(c$error))
    as.data.frame(round(c$res$information, 5))
  }, rownames = TRUE)

  # ---- model summary text -------------------------------------------------
  output$model_txt <- renderPrint({
    c <- computed(); req(is.null(c$error))
    model_summary(c$res)
  })

  # ---- downloads ----------------------------------------------------------
  output$dl_csv <- downloadHandler(
    filename = function() "owea_design.csv",
    content  = function(file) {
      old <- options(digits = 15); on.exit(options(old))   # full-precision round-trip
      utils::write.csv(design_df(), file, row.names = FALSE)
    })
  output$dl_png <- downloadHandler(
    filename = function() "owea_design.png",
    content  = function(file) {
      grDevices::png(file, width = 900, height = 650, res = 110)
      on.exit(grDevices::dev.off())
      c <- computed()
      if (is.null(c$error))
        plot_design(c$res, cov_names = c$cov_names)
    })

  # ---- post-result actions: Verify (approx) or Simulation (exact) + Reset --
  rv <- reactiveValues(verify = NULL, sim = list(), show_sim = FALSE)
  observeEvent(input$compute, { rv$verify <- NULL; rv$sim <- list(); rv$show_sim <- FALSE })
  observeEvent(input$reset_btn, session$reload())

  output$post_result <- renderUI({
    c <- computed(); req(is.null(c$error))
    if (c$exact)
      tagList(hr(),
        actionButton("sim_open", "Simulation study", class = "btn-info"),
        actionButton("reset_btn", "Reset"),
        uiOutput("sim_panel"))
    else
      tagList(hr(), uiOutput("verify_panel"))
  })

  # ---- Verify optimality of a GIVEN design (approximate design) -----------
  output$verify_panel <- renderUI({
    c <- computed(); req(is.null(c$error)); if (isTRUE(c$exact)) return(NULL)
    tagList(
      tags$h4("Verify optimality of a design"),
      helpText("Provide a design in the SAME format as the downloaded design CSV ",
               "(covariate columns + a 'weight' column). The box is prefilled with the ",
               "computed design — edit it, or upload a CSV, to verify a different design."),
      fileInput("verify_file", "Upload design CSV (optional)", accept = ".csv"),
      textAreaInput("verify_manual", "… or paste / edit the design",
                    value = design_csv_text(design_df()), rows = 6),
      actionButton("verify_btn", "Verify optimality", class = "btn-info"),
      actionButton("reset_btn", "Reset"),
      uiOutput("verify_out"))
  })

  observeEvent(input$verify_btn, {
    c <- computed(); if (is.null(c) || !is.null(c$error)) return()
    sp <- c$sp; warns <- character(0)
    v <- withCallingHandlers(
      tryCatch({
        dz <- read_design(input$verify_file, input$verify_manual, "weight", c$cov_names)
        w  <- as.numeric(dz$val); sw <- sum(w)
        if (!is.finite(sw) || sw <= 0)
          stop("weights must be positive and sum to a positive value.", call. = FALSE)
        if (abs(sw - 1) > 1e-6)
          warning(sprintf("Weights summed to %.6f; rescaled to sum to 1.", sw), call. = FALSE)
        verify_optimality(dz$support, w / sw,
          design_box = sp$design_box, step = c$finest,
          link = sp$link, f = sp$f, x = sp$x, fx = sp$fx, xx = sp$xx, ff = sp$ff,
          ncat = sp$ncat, theta = c$theta, p = c$p, subset = c$subset,
          max_points = Inf)
      },
      error = function(e) structure(list(msg = conditionMessage(e)), class = "owea_bad")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    rv$verify <- list(v = v, warns = warns, opt_crit = c$res$criterion)
  })

  output$verify_out <- renderUI({
    vr <- rv$verify; if (is.null(vr)) return(NULL)
    if (inherits(vr$v, "owea_bad"))
      return(div(class = "alert alert-danger", tags$b("Verify failed: "), vr$v$msg))
    v <- vr$v; opt <- isTRUE(v$is_optimal$value)
    tagList(br(),
      if (length(vr$warns))
        div(class = "alert alert-warning", tags$b("Warning: "),
            paste(unique(vr$warns), collapse = "  ")),
      div(class = if (opt) "alert alert-success" else "alert alert-warning",
          tags$p(tags$b("Max sensitivity: "),
                 sprintf("%.3e", v$max_sensitivity), " (0 at the optimum)"),
          tags$p(tags$b("Criterion (this design): "), sprintf("%.6f", v$criterion)),
          tags$p(tags$small(sprintf("Optimal design criterion (for reference): %.6f",
                                    vr$opt_crit))),
          tags$p(tags$b(if (opt) "The design IS optimal."
                        else "The design is NOT certified optimal.")),
          tags$p(tags$small(v$is_optimal$note))))
  })

  # ---- Simulation study (exact design) ------------------------------------
  observeEvent(input$sim_open, { rv$show_sim <- TRUE })

  output$sim_panel <- renderUI({
    if (!isTRUE(rv$show_sim)) return(NULL)
    c <- computed(); req(is.null(c$error)); cn <- c$coef_names
    tagList(hr(), tags$h4("Simulation study"),
      helpText(sprintf("Total sample size n = %d (matched across designs).", c$n)),
      tags$b("True parameter values (theta)"),
      do.call(tagList, lapply(seq_along(cn), function(i)
        numericInput(paste0("sim_theta_", i), cn[i],
                     value = if (!is.null(c$theta)) c$theta[i] else 0, step = 0.1))),
      if (identical(c$sp$link, "identity"))
        numericInput("sim_sigma", "residual SD (sigma)", value = 1, min = 1e-6),
      fluidRow(
        column(6, numericInput("sim_nsim", "simulations", value = 1000, min = 2, step = 1)),
        column(6, numericInput("sim_seed", "seed", value = 1, step = 1))),
      actionButton("run_sim", "Run simulation (exact design)", class = "btn-primary"),
      hr(),
      tags$b("Compare against another design of the same n (optional)"), br(), br(),
      actionButton("run_srs", "Run simple random sample (size n)"), br(), br(),
      helpText("Custom design: use the SAME format as the downloaded design CSV ",
               "(covariate columns + a 'count' column summing to n)."),
      fileInput("sim_custom_file", "Upload a custom design CSV (optional)", accept = ".csv"),
      textAreaInput("sim_custom", "… or paste a custom design",
        rows = 4, placeholder = "dose,conc,count\n-1,0.5,10\n1,-0.5,10"),
      actionButton("run_custom", "Run custom design"),
      br(), br(), tableOutput("sim_tbl"), uiOutput("sim_msg"))
  })

  # run a simulation for a given support/counts (or an exact_design object)
  sim_run <- function(support, counts) {
    c <- computed(); sp <- c$sp
    theta_true <- vapply(seq_along(c$coef_names),
      function(i) as.numeric(input[[paste0("sim_theta_", i)]] %||% 0), numeric(1))
    sigma <- if (identical(sp$link, "identity")) as.numeric(input$sim_sigma %||% 1) else 1
    simulate_design(support, counts = counts, design_box = sp$design_box,
      link = sp$link, f = sp$f, x = sp$x, fx = sp$fx, xx = sp$xx, ff = sp$ff,
      ncat = sp$ncat, theta = theta_true, sigma = sigma,
      nsim = as.integer(input$sim_nsim %||% 1000),
      seed = as.integer(input$sim_seed %||% 1))
  }

  observeEvent(input$run_sim, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim
    s$exact <- tryCatch(sim_run(c$res, NULL), error = function(e) e)
    rv$sim <- s
  })

  observeEvent(input$run_srs, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim; s$pool_warn <- NULL
    r <- tryCatch({
      pool <- candidate_grid(c$sp$design_box, c$finest)
      if (nrow(pool) > 1e6)
        s$pool_warn <- sprintf(paste0("The SRS candidate pool has %d points; consider ",
          "coarsening the step."), nrow(pool))
      idx  <- sample(nrow(pool), c$n, replace = TRUE)
      pts  <- pool[idx, , drop = FALSE]
      keys <- do.call(paste, c(as.data.frame(pts), sep = "\r"))
      uk   <- unique(keys)
      sup  <- pts[match(uk, keys), , drop = FALSE]
      cnt  <- as.integer(table(factor(keys, levels = uk)))
      sim_run(sup, cnt)
    }, error = function(e) e)
    s$srs <- r; rv$sim <- s
  })

  observeEvent(input$run_custom, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim
    r <- tryCatch({
      dz  <- read_design(input$sim_custom_file, input$sim_custom, "count", c$cov_names)
      cnt <- as.integer(dz$val)
      if (any(is.na(cnt))) stop("the 'count' column must be integers.", call. = FALSE)
      if (sum(cnt) != c$n)
        stop(sprintf("counts sum to %d but must sum to n = %d.", sum(cnt), c$n), call. = FALSE)
      sim_run(dz$support, cnt)
    }, error = function(e) e)
    s$custom <- r; rv$sim <- s
  })

  output$sim_tbl <- renderTable({
    c <- computed(); req(is.null(c$error)); s <- rv$sim
    getmse <- function(r) if (is.null(r) || inherits(r, "error") || is.null(r$mse)) NULL
                          else as.numeric(r$mse)
    e <- getmse(s$exact); sr <- getmse(s$srs); cu <- getmse(s$custom)
    if (is.null(e) && is.null(sr) && is.null(cu)) return(NULL)
    df <- data.frame(parameter = c$coef_names, stringsAsFactors = FALSE)
    if (!is.null(e))  df[["Exact design (MSE)"]] <- e
    if (!is.null(sr)) df[["SRS (MSE)"]] <- sr
    if (!is.null(cu)) df[["Custom (MSE)"]] <- cu
    df
  }, digits = 6)

  output$sim_msg <- renderUI({
    s <- rv$sim; msgs <- character(0)
    for (nm in c("exact", "srs", "custom"))
      if (inherits(s[[nm]], "error"))
        msgs <- c(msgs, sprintf("%s: %s", nm, conditionMessage(s[[nm]])))
    if (!is.null(s$pool_warn)) msgs <- c(msgs, s$pool_warn)
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-warning", lapply(msgs, tags$div))
  })
}

# small helpers ------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# a design data.frame -> CSV text (to prefill the verify box), full precision.
design_csv_text <- function(df) {
  old <- options(digits = 15); on.exit(options(old))
  paste(utils::capture.output(utils::write.csv(df, row.names = FALSE)), collapse = "\n")
}

# Read a design from an uploaded CSV (priority) or a pasted text box, in the
# downloaded-design-CSV format: covariate columns + one value column
# (`valcol` = "weight" or "count").  Returns list(support = matrix, val = vector).
read_design <- function(file, text, valcol, cov_names) {
  d <- if (!is.null(file) && !is.null(file$datapath) && nzchar(file$datapath))
         utils::read.csv(file$datapath, header = TRUE, check.names = FALSE)
       else {
         if (!nzchar(trimws(text %||% "")))
           stop("provide a design: upload a CSV or paste one.", call. = FALSE)
         utils::read.csv(text = text, header = TRUE, check.names = FALSE)
       }
  if (!valcol %in% names(d))
    stop(sprintf("the design needs a '%s' column (see the downloaded CSV format).", valcol),
         call. = FALSE)
  val <- as.numeric(d[[valcol]])
  if (all(cov_names %in% names(d))) {
    sup <- as.matrix(d[, cov_names, drop = FALSE])          # match by name
  } else {
    cols <- setdiff(names(d), valcol)                       # positional fallback
    if (length(cols) != length(cov_names))
      stop(sprintf("expected %d covariate column(s) plus '%s'.",
                   length(cov_names), valcol), call. = FALSE)
    sup <- as.matrix(d[, cols, drop = FALSE])
  }
  storage.mode(sup) <- "double"
  list(support = sup, val = val)
}

# turn owea's technical warnings into plain language for the banner.
.friendly_warn <- function(w) {
  vapply(w, function(m) {
    if (grepl("convergence is LOCAL", m, fixed = TRUE))
      "Optimality was certified only near the current support; enable a finer grid for a whole-region guarantee."
    else if (grepl("drawn from", m) || grepl("N\\(0", m))
      "No parameter values were given, so random ones were used — set theta for a meaningful local design."
    else if (grepl("did NOT converge", m, fixed = TRUE))
      "The solver did not fully converge; try a finer grid step or fewer parameters of interest."
    else m
  }, character(1), USE.NAMES = FALSE)
}

shinyApp(ui, server)
