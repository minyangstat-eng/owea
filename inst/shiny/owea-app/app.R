# ===========================================================================
# owea web app -- point-and-click optimal designs, as a step-by-step wizard.
#
#   model -> start -> [existing design | existing data] -> [theta] ->
#   criterion -> design type -> review -> results
#
# Steps that do not apply are skipped: no theta step for the identity link (its
# design is not local), and no existing-design step unless the user has one.
#
# Launch locally with owea::run_owea_app(), or deploy this folder to
# shinyapps.io / Posit Connect (see deploy/deploy.R in the source repo).
#
# The app depends only on the exported owea API plus the internal helpers in
# R/app.R (.ui_model_spec, .ui_wizard_steps, .ui_solver_args, ...), which
# translate the friendly inputs into optimal_design()/exact_design() args and
# are unit-tested independently of the app.  In particular EVERY solver call is
# assembled by .ui_solver_args(), so the model -- above all the factor `coding`
# -- can never drift between the design, the fit, the simulation and the
# verification.
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
CRITERIA  <- c("D-optimal (overall precision)" = "0",
               "A-optimal (average variance)"  = "1")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("owea — Optimal Experimental Designs"),
  tags$p("Build a D- or A-optimal design step by step: describe your model, ",
         "say whether you already have a design or a data set, set the assumed ",
         "parameter values and the criterion, then compute. No R code required."),
  uiOutput("progress"),
  hr(),

  tabsetPanel(
    id = "wiz", type = "hidden",

    # ---- 1. model ---------------------------------------------------------
    tabPanelBody(
      "model",
      h4("Step: the model"),
      fluidRow(
        column(
          5,
          selectInput("link", "Model family", choices = LINKS),
          conditionalPanel(
            "input.link == 'multinomial' || input.link == 'cumulative'",
            numericInput("ncat", "Number of response categories", value = 3,
                         min = 2, step = 1)),
          conditionalPanel(
            "input.link != 'identity'",
            div(class = "alert alert-info",
                tags$b("This will be a locally optimal design."), " For a ",
                "non-linear model the information matrix depends on the ",
                "parameters themselves, so the design is optimal only at an ",
                "assumed parameter value. You will be asked for those values ",
                "(or you can draw them at random) in a later step.")),
          uiOutput("coding_ui")
        ),
        column(
          7,
          numericInput("ncov", "Number of covariates", value = 1, min = 1,
                       max = MAX_COV, step = 1),
          uiOutput("cov_ui"),
          uiOutput("interaction_ui"),
          uiOutput("quadratic_ui")
        )
      )
    ),

    # ---- 2. starting point ------------------------------------------------
    tabPanelBody(
      "start",
      h4("Step: what do you already have?"),
      radioButtons(
        "start", NULL,
        choices = c("Nothing yet — design from scratch"          = "none",
                    "An existing design (support points + runs)" = "design",
                    "An existing data set (covariates + responses)" = "data"),
        selected = "none"),
      helpText("With an existing design or data set the new design is chosen to ",
               "complement what you already have (a second-stage design).")
    ),

    # ---- 3a. existing design ----------------------------------------------
    tabPanelBody(
      "design_in",
      h4("Step: your existing design"),
      helpText("Covariate columns plus a 'count' column (integer runs) or a ",
               "'weight' column (proportions) — the format the app downloads. ",
               "Factor covariates use integer levels 1..L."),
      fileInput("exist_file", "Upload a design CSV", accept = ".csv"),
      textAreaInput("exist_text", "… or paste / edit the design", rows = 6,
                    placeholder = "dose,count\n-1,10\n1,10"),
      numericInput("exist_n0", "n0 — sample size of the existing design",
                   value = NA, min = 1, step = 1),
      helpText("With a 'count' column, n0 is filled in from the counts. ",
               "With weights you must supply it: proportions carry no sample size."),
      uiOutput("exist_msg")
    ),

    # ---- 3b. existing data set --------------------------------------------
    tabPanelBody(
      "data_in",
      h4("Step: your existing data set"),
      helpText("Covariate columns plus a response column, one row per run. ",
               "Factor covariates use integer levels 1..L."),
      fileInput("data_file", "Upload a data CSV", accept = ".csv"),
      textAreaInput("data_text", "… or paste the data", rows = 6,
                    placeholder = "dose,y\n-1,0\n-1,1\n1,1"),
      uiOutput("data_response_ui"),
      hr(),
      radioButtons("use_cov",
                   "Use the data set's covariates as the existing design?",
                   choices = c("Yes" = "yes", "No" = "no"), selected = "yes",
                   inline = TRUE),
      radioButtons("use_theta",
                   "Estimate the assumed parameter values from the data set?",
                   choices = c("Yes" = "yes", "No" = "no"), selected = "yes",
                   inline = TRUE),
      helpText("Estimating the parameters fits the model you described in step 1 ",
               "to these data (maximum likelihood) and uses the estimates as the ",
               "assumed values — you can still edit them afterwards. The data set ",
               "is read as soon as you upload or paste it; the button below just ",
               "reads it again."),
      actionButton("data_load", "Re-read the data set"),
      uiOutput("data_msg"),
      uiOutput("data_existing_out"),
      uiOutput("data_fit_out")
    ),

    # ---- 4. assumed parameter values --------------------------------------
    tabPanelBody(
      "theta",
      h4("Step: assumed parameter values (theta)"),
      helpText("A locally optimal design depends on these values. Type them in, ",
               "or draw them at random and edit."),
      actionButton("theta_draw", "Draw from N(0,1)"),
      br(), br(),
      uiOutput("theta_ui"),
      uiOutput("theta_msg")
    ),

    # ---- 5. criterion -----------------------------------------------------
    tabPanelBody(
      "criterion",
      h4("Step: the criterion"),
      selectInput("crit", "Optimality criterion", choices = CRITERIA),
      radioButtons("qoi", "Parameters of interest",
                   choices = c("All parameters" = "all", "A subset" = "subset")),
      conditionalPanel("input.qoi == 'subset'", uiOutput("subset_ui"))
    ),

    # ---- 6. design type ---------------------------------------------------
    tabPanelBody(
      "design_type",
      h4("Step: the design"),
      radioButtons("design_type", "Design type",
                   choices = c("Approximate (weights)" = "approx",
                               "Exact (integer runs)"  = "exact")),
      uiOutput("n_ui"),
      conditionalPanel("input.design_type == 'exact'",
                       numericInput("seed", "Random seed", value = 1, step = 1))
    ),

    # ---- 7. review --------------------------------------------------------
    tabPanelBody(
      "review",
      h4("Step: review and compute"),
      uiOutput("review_out"),
      br(),
      actionButton("compute", "Compute design", class = "btn-primary",
                   width = "50%")
    ),

    # ---- 8. results -------------------------------------------------------
    tabPanelBody(
      "results",
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Design",
                 br(), DT::DTOutput("design_tbl"),
                 uiOutput("crit_note"),
                 br(), downloadButton("dl_csv", "Download design (CSV)")),
        tabPanel("Plot",
                 br(), plotOutput("design_plot", height = "460px"),
                 br(), downloadButton("dl_png", "Download plot (PNG)")),
        tabPanel("Information matrix",
                 br(), uiOutput("info_help"),
                 tableOutput("info_tbl")),
        tabPanel("Model",
                 br(), verbatimTextOutput("model_txt"))
      ),
      uiOutput("eff_panel"),
      uiOutput("post_result")
    )
  ),

  hr(),
  fluidRow(
    column(3, uiOutput("back_ui")),
    column(3, uiOutput("next_ui")),
    column(6, uiOutput("restart_ui"))
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    verify = NULL, sim = list(), show_sim = FALSE,   # post-result panels
    theta_prefill = NULL,                            # theta boxes' values
    data = NULL, data_existing = NULL, data_fit = NULL, data_msg = NULL,
    eff = NULL)

  # ---- model step ---------------------------------------------------------
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
            column(4, textInput(paste0("cov_step_", i), "grid step(s)",
                                value = "0.1", placeholder = "e.g. 0.5, 0.1"))),
          helpText("One number, or a comma-separated step sequence from coarse ",
                   "to fine (e.g. 0.5, 0.1, 0.02): the whole range is searched ",
                   "at the coarsest step, then refined locally at each finer one.")),
        conditionalPanel(
          sprintf("input.cov_type_%d == 'factor'", i),
          numericInput(paste0("cov_nlev_", i), "number of levels", value = 2,
                       min = 2, step = 1))
      )
    })
  })

  covariates <- reactive({
    n <- max(1L, min(MAX_COV, as.integer(input$ncov %||% 1)))
    lapply(seq_len(n), function(i) {
      type <- input[[paste0("cov_type_", i)]] %||% "continuous"
      nm   <- input[[paste0("cov_name_", i)]] %||% paste0("x", i)
      if (identical(type, "factor")) {
        list(name = nm, type = "factor",
             nlevels = as.integer(input[[paste0("cov_nlev_", i)]] %||% 2))
      } else {
        list(name = nm, type = "continuous",
             lo    = as.numeric(input[[paste0("cov_lo_", i)]] %||% -1),
             hi    = as.numeric(input[[paste0("cov_hi_", i)]] %||%  1),
             steps = owea:::.ui_parse_steps(input[[paste0("cov_step_", i)]] %||% "0.1"))
      }
    })
  })

  has_factor <- reactive(
    any(vapply(covariates(), function(cv) identical(cv$type, "factor"), logical(1))))

  # the factor coding only exists when there IS a factor (item 2)
  output$coding_ui <- renderUI({
    if (!has_factor()) return(NULL)
    tagList(
      radioButtons("coding", "Coding of the factor levels",
                   choices = c("Zero-sum (effect coding)" = "zero-sum",
                               "Baseline (reference level)" = "baseline"),
                   selected = "zero-sum"),
      helpText("The two codings describe the same model but parameterise it ",
               "differently, so the parameter values (and any estimates from a ",
               "data set) refer to whichever you pick here."))
  })

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

  spec <- reactive({
    inter <- lapply(input$interactions %||% character(0),
                    function(s) as.integer(strsplit(s, "-", fixed = TRUE)[[1]]))
    quad  <- as.integer(input$quadratics %||% character(0))
    ncat  <- if (input$link %in% MULTI_CAT) as.integer(input$ncat %||% 3) else NULL
    owea:::.ui_model_spec(covariates(), interactions = inter, quadratics = quad,
                          link = input$link, ncat = ncat,
                          coding = if (has_factor()) input$coding %||% "zero-sum"
                                   else "zero-sum")
  })
  spec_ok    <- reactive(tryCatch({ spec(); NULL }, error = conditionMessage))
  coef_names <- reactive(tryCatch(owea:::.ui_coef_names(spec()),
                                  error = function(e) NULL))
  cov_names  <- reactive(tryCatch(names(spec()$design_box),
                                  error = function(e) NULL))

  # ---- navigation ---------------------------------------------------------
  cur   <- reactiveVal("model")
  steps <- reactive(owea:::.ui_wizard_steps(input$link %||% "identity",
                                            input$start %||% "none"))
  step_i <- reactive({ i <- match(cur(), steps()); if (is.na(i)) 1L else i })

  output$progress <- renderUI({
    s <- steps(); i <- step_i()
    tagList(
      tags$b(sprintf("Step %d of %d", i, length(s))),
      div(class = "progress", style = "height:6px;",
          div(class = "progress-bar",
              style = sprintf("width:%.0f%%;", 100 * i / length(s)))))
  })
  output$back_ui <- renderUI(
    if (step_i() > 1L) actionButton("back_btn", "← Back", width = "100%"))
  output$next_ui <- renderUI(
    if (!identical(cur(), "review") && !identical(cur(), "results"))
      actionButton("next_btn", "Next →", class = "btn-primary", width = "100%"))
  output$restart_ui <- renderUI(
    if (identical(cur(), "results"))
      actionButton("restart_btn", "Start over", width = "40%"))

  goto <- function(id) { cur(id); updateTabsetPanel(session, "wiz", selected = id) }

  # per-step validation: a non-NULL message blocks Next
  step_error <- function(id) {
    if (identical(id, "model")) return(spec_ok())
    if (identical(id, "design_in"))
      return(tryCatch({ existing(); NULL }, error = conditionMessage))
    if (identical(id, "data_in")) {
      if (identical(input$use_cov, "no") && identical(input$use_theta, "no"))
        return(NULL)                                  # nothing needed from it
      if (!data_given())
        return("upload or paste your data set (or answer 'No' to both questions).")
      if (is.null(rv$data))
        return(rv$data_msg %||% "the data set could not be read.")
      if (identical(input$use_theta, "yes") && is.null(rv$data_fit))
        return("the model could not be fitted to these data — see the message above.")
      if (identical(input$use_cov, "yes") && is.null(rv$data_existing))
        return("the covariates could not be used as a design — see the message above.")
      return(NULL)
    }
    if (identical(id, "theta"))
      return(tryCatch(owea:::.ui_check_theta(theta(), spec()),
                      error = conditionMessage))
    if (identical(id, "design_type")) {
      if (need_n() && !is.finite(as.numeric(input$n_new %||% NA)))
        return("enter a sample size.")
      return(NULL)
    }
    NULL
  }

  # candidate-set size at each stage of the step sequence (coarsest first);
  # grid_ok remembers the counts the user already agreed to proceed with, so
  # the warning re-arms only when the grid inputs actually change
  grid_sizes <- reactive(tryCatch(owea:::.ui_grid_sizes(covariates()),
                                  error = function(e) NULL))
  grid_ok <- reactiveVal(NULL)

  advance <- function() {
    s <- steps(); i <- step_i()
    if (i < length(s)) goto(s[i + 1L])
  }

  observeEvent(input$next_btn, {
    err <- step_error(cur())
    if (!is.null(err)) {
      showNotification(err, type = "error", duration = 8); return()
    }
    if (identical(cur(), "model")) {
      gs <- grid_sizes(); N1 <- gs[1]              # first (coarsest) stage
      if (length(gs) && is.finite(N1) && N1 > 1e6 && !identical(grid_ok(), gs)) {
        n_txt <- format(round(N1), big.mark = ",")
        if (length(gs) == 1L)
          showModal(modalDialog(title = "Large candidate set",
            sprintf(paste0("The candidate set has %s design points at this grid ",
              "step; building it may be slow."), n_txt),
            footer = tagList(
              modalButton("Adjust the grid size"),
              actionButton("use_seq", "Use a step sequence"),
              actionButton("grid_proceed", "Proceed anyway",
                           class = "btn-warning"))))
        else
          showModal(modalDialog(title = "Large candidate set",
            sprintf(paste0("The candidate set has %s design points at the first ",
              "(coarsest) step of your step sequence; building it may be slow."),
              n_txt),
            footer = tagList(
              modalButton("Adjust the step sequence"),
              actionButton("grid_proceed", "Proceed anyway",
                           class = "btn-warning"))))
        return()
      }
    }
    advance()
  })
  observeEvent(input$use_seq, {
    removeModal()
    showNotification(paste0(
      "In the 'grid step(s)' box enter a comma-separated sequence of step sizes ",
      "from coarse to fine, e.g. 0.5, 0.1, 0.02. Make the first step coarse ",
      "enough that its grid stays below 1,000,000 points; the finest step sets ",
      "the final precision."), type = "message", duration = NULL)
  })
  observeEvent(input$grid_proceed, {
    removeModal(); grid_ok(grid_sizes()); advance()
  })
  observeEvent(input$back_btn, {
    s <- steps(); i <- step_i()
    if (i > 1L) goto(s[i - 1L])
  })
  observeEvent(input$restart_btn, session$reload())
  observeEvent(input$reset_btn,   session$reload())

  # changing the link or the starting point can remove the step we are standing
  # on (e.g. switching to the identity link while on the theta step)
  observeEvent(list(input$link, input$start), {
    if (!cur() %in% steps()) goto("model")
  }, ignoreInit = TRUE)

  # ---- existing design (pasted / uploaded CSV) ----------------------------
  # never throws: an empty or malformed box is an error VALUE, so that neither
  # the observer below nor an unrelated step blows up on it
  exist_csv <- reactive({
    tryCatch(read_design(input$exist_file, input$exist_text, cov_names(),
                         valcol = NULL),
             error = function(e) e)
  })
  # counts carry a sample size -- offer it as n0
  observeEvent(exist_csv(), {
    d <- exist_csv()
    if (inherits(d, "error")) return()
    if (identical(d$valcol, "count") && !is.finite(as.numeric(input$exist_n0 %||% NA)))
      updateNumericInput(session, "exist_n0", value = sum(round(d$val)))
  }, ignoreInit = TRUE)

  # the existing design in solver form: list(points, weights, n0), or NULL
  existing <- reactive({
    st <- input$start %||% "none"
    if (identical(st, "design")) {
      d  <- exist_csv()
      if (inherits(d, "error")) stop(conditionMessage(d), call. = FALSE)
      n0 <- as.numeric(input$exist_n0 %||% NA)
      e  <- owea:::.ui_existing_from_csv(d$support, d$val, d$valcol,
                                         n0 = if (is.finite(n0)) n0 else NULL)
      if (is.null(e$n0))
        stop("enter n0 — the sample size of the existing design. Weights are ",
             "proportions and carry no sample size.", call. = FALSE)
      e
    } else if (identical(st, "data") && identical(input$use_cov, "yes")) {
      rv$data_existing
    } else {
      NULL
    }
  })
  has_existing <- reactive(!is.null(tryCatch(existing(), error = function(e) NULL)))

  output$exist_msg <- renderUI({
    e <- tryCatch(existing(), error = function(e) e)
    if (inherits(e, "error"))
      return(div(class = "alert alert-warning", conditionMessage(e)))
    if (is.null(e)) return(NULL)
    div(class = "alert alert-success",
        sprintf("%d support point(s), n0 = %d.", nrow(e$points),
                as.integer(e$n0)),
        if (length(e$notes)) tags$ul(lapply(e$notes, tags$li)))
  })

  # ---- existing data set --------------------------------------------------
  # has the user actually given us data yet?
  data_given <- reactive({
    f <- input$data_file
    (!is.null(f) && !is.null(f$datapath) && nzchar(f$datapath)) ||
      nzchar(trimws(input$data_text %||% ""))
  })
  data_raw <- reactive({
    if (!data_given()) return(NULL)
    tryCatch(read_table(input$data_file, input$data_text), error = function(e) NULL)
  })
  output$data_response_ui <- renderUI({
    d <- data_raw(); if (is.null(d)) return(NULL)
    nm  <- names(d)
    sel <- if ("y" %in% nm) "y" else nm[length(nm)]
    selectInput("data_response", "Response column", choices = nm, selected = sel)
  })

  # Read the data set, use its covariates as the existing design and/or fit the
  # model to it -- whichever the two answers ask for.
  load_data <- function() {
    rv$data <- NULL; rv$data_existing <- NULL; rv$data_fit <- NULL
    rv$data_msg <- NULL
    if (!data_given()) return()             # nothing uploaded/pasted yet
    msgs <- character(0)
    d <- tryCatch(read_table(input$data_file, input$data_text),
                  error = function(e) e)
    if (inherits(d, "error")) { rv$data_msg <- conditionMessage(d); return() }
    rv$data <- d
    sp <- tryCatch(spec(), error = function(e) e)
    if (inherits(sp, "error")) { rv$data_msg <- conditionMessage(sp); return() }
    resp <- input$data_response %||% NULL

    if (identical(input$use_cov, "yes")) {
      ex <- tryCatch(owea:::.ui_existing_from_data(d, sp$design_box, resp),
                     error = function(e) e)
      if (inherits(ex, "error"))
        msgs <- c(msgs, paste("Covariates as a design:", conditionMessage(ex)))
      else rv$data_existing <- ex
    }
    if (identical(input$use_theta, "yes")) {
      fit <- tryCatch(
        suppressWarnings(suppressMessages(do.call(
          fit_design,
          c(list(data = d, response = resp),
            owea:::.ui_solver_args(sp, "fit"))))),
        error = function(e) e)
      if (inherits(fit, "error"))
        msgs <- c(msgs, paste("Fitting the model:", conditionMessage(fit)))
      else {
        rv$data_fit <- fit
        rv$theta_prefill <- as.numeric(fit$theta_hat)   # fills the theta boxes
      }
    }
    rv$data_msg <- if (length(msgs)) paste(msgs, collapse = "  ") else NULL
  }

  # The upload IS the load: process the data as soon as it arrives, and again
  # whenever anything it depends on changes (the response column, either answer,
  # or the model itself -- otherwise the fit would go stale).  Pasted text is
  # debounced so we do not refit on every keystroke.
  data_typed <- debounce(reactive(list(input$data_file$datapath, input$data_text)),
                         700)
  # NOT ignoreInit: the data can already be there the first time this runs (all
  # the inputs arriving in one flush, a restored session), and skipping that run
  # would leave it unread.  load_data() is a no-op until data is actually given.
  observeEvent(
    list(data_typed(), input$data_response, input$use_cov, input$use_theta,
         coef_names()),
    { if (identical(input$start, "data")) load_data() })
  observeEvent(input$data_load, load_data())   # explicit re-read, never required

  output$data_msg <- renderUI({
    if (!is.null(rv$data_msg))
      return(div(class = "alert alert-danger", rv$data_msg))
    if (is.null(rv$data)) return(NULL)
    div(class = "alert alert-success",
        sprintf("Loaded %d observation(s).", nrow(rv$data)))
  })

  output$data_existing_out <- renderUI({
    ex <- rv$data_existing; if (is.null(ex)) return(NULL)
    tagList(
      tags$b("The existing design taken from these covariates"),
      helpText(sprintf("n0 = %d (the number of observations).", ex$n0)),
      tags$pre(design_csv_text(ex$df)))
  })

  output$data_fit_tbl <- renderTable({
    fit <- rv$data_fit; req(!is.null(fit))
    data.frame(parameter    = fit$coef_names,
               estimate     = as.numeric(fit$theta_hat),
               `std. error` = as.numeric(fit$se),
               check.names  = FALSE)
  }, digits = 4)

  output$data_fit_out <- renderUI({
    fit <- rv$data_fit; if (is.null(fit)) return(NULL)
    tagList(
      tags$b("Estimated parameter values (used as the assumed values)"),
      if (!isTRUE(fit$converged))
        div(class = "alert alert-warning", "The fit did not converge."),
      tableOutput("data_fit_tbl"),
      helpText("You can edit these in the next step."))
  })

  # ---- assumed parameter values ------------------------------------------
  output$theta_ui <- renderUI({
    cn <- coef_names()
    if (is.null(cn)) return(helpText("(define the model first)"))
    pre <- rv$theta_prefill
    lapply(seq_along(cn), function(i)
      numericInput(paste0("theta_", i), cn[i],
                   value = if (!is.null(pre) && length(pre) >= i) pre[i] else 0,
                   step = 0.1))
  })
  # mirrors what the theta boxes SHOW: the box's value, else the prefill they
  # were rendered with (a fit or a draw), else 0 -- the box's own default
  theta <- reactive({
    cn <- coef_names(); if (is.null(cn)) return(NULL)
    if (identical(input$link, "identity")) return(NULL)
    pre <- rv$theta_prefill
    vapply(seq_along(cn), function(i) {
      v <- input[[paste0("theta_", i)]]
      if (!is.null(v) && !is.na(v)) return(as.numeric(v))
      if (!is.null(pre) && length(pre) >= i) as.numeric(pre[i]) else 0
    }, numeric(1))
  })
  observeEvent(input$theta_draw, {
    sp <- tryCatch(spec(), error = function(e) NULL); req(sp)
    th <- owea:::.ui_random_theta(sp); req(!is.null(th))
    rv$theta_prefill <- th                       # observable, and survives a re-render
    for (i in seq_along(th))
      updateNumericInput(session, paste0("theta_", i), value = round(th[i], 4))
  })
  output$theta_msg <- renderUI({
    msg <- tryCatch(owea:::.ui_check_theta(theta(), spec()),
                    error = conditionMessage)
    if (is.null(msg)) return(NULL)
    div(class = "alert alert-warning", msg)
  })

  # ---- criterion ----------------------------------------------------------
  output$subset_ui <- renderUI({
    cn <- coef_names(); if (is.null(cn)) return(NULL)
    checkboxGroupInput("subset_terms", "Include parameters", choices = cn)
  })
  subset_idx <- reactive({
    cn <- coef_names()
    if (!identical(input$qoi, "subset") || is.null(cn)) return(NULL)
    sel <- input$subset_terms %||% character(0)
    if (!length(sel)) return(NULL)
    match(sel, cn)
  })

  # ---- design type --------------------------------------------------------
  is_exact <- reactive(identical(input$design_type, "exact"))
  # an approximate design with no existing design needs no sample size at all
  need_n   <- reactive(is_exact() || has_existing())

  output$n_ui <- renderUI({
    if (!need_n()) return(helpText("An approximate design needs no sample size."))
    lab <- if (is_exact() && has_existing())
             "n — runs in the new stage (also used as n1)"
           else if (is_exact()) "n — total runs"
           else "n1 — sample size of the new stage"
    hlp <- if (is_exact() && has_existing())
             paste("The new stage's runs and the n1 used to weight the existing",
                   "design are tied together, so the reported criterion always",
                   "refers to the design you actually get.")
           else if (has_existing())
             "n0 : n1 sets how heavily the existing design counts."
           else NULL
    tagList(
      numericInput("n_new", lab, value = isolate(input$n_new) %||% 20,
                   min = 1, step = 1),
      if (!is.null(hlp)) helpText(hlp))
  })

  # ---- review -------------------------------------------------------------
  output$review_out <- renderUI({
    sp <- tryCatch(spec(), error = function(e) e)
    if (inherits(sp, "error"))
      return(div(class = "alert alert-danger", conditionMessage(sp)))
    ex <- tryCatch(existing(), error = function(e) NULL)
    gs <- grid_sizes()
    li <- list(
      sprintf("Model: %s%s", names(LINKS)[match(sp$link, LINKS)],
              if (is.null(sp$ncat)) "" else sprintf(" with %d categories", sp$ncat)),
      sprintf("Covariates: %s", paste(names(sp$design_box), collapse = ", ")),
      if (length(gs) && is.finite(gs[1]))
        sprintf("Candidate set: %s design points%s",
                format(round(gs[1]), big.mark = ","),
                if (length(gs) > 1L) " (first step of the step sequence)" else ""),
      if (has_factor()) sprintf("Factor coding: %s", owea:::.ui_coding(sp)),
      sprintf("Parameters (%d): %s", length(coef_names()),
              paste(coef_names(), collapse = ", ")),
      if (!identical(sp$link, "identity"))
        sprintf("Assumed theta: %s",
                paste(round(theta(), 4), collapse = ", ")),
      sprintf("Criterion: %s", names(CRITERIA)[match(input$crit, CRITERIA)]),
      sprintf("Parameters of interest: %s",
              if (is.null(subset_idx())) "all"
              else paste(coef_names()[subset_idx()], collapse = ", ")),
      if (!is.null(ex))
        sprintf("Existing design: %d support point(s), n0 = %d",
                nrow(ex$points), ex$n0),
      sprintf("Design type: %s", if (is_exact()) "exact (integer runs)"
                                 else "approximate (weights)"),
      if (need_n()) sprintf("Sample size: %s", input$n_new %||% "—"))
    tags$ul(lapply(Filter(Negate(is.null), li), tags$li))
  })

  # ---- compute (the large-grid warning now gates the model step's Next) ----
  go <- reactiveVal(0)
  observeEvent(input$compute, go(go() + 1))
  observeEvent(input$compute, {
    rv$verify <- NULL; rv$sim <- list(); rv$show_sim <- FALSE; rv$eff <- NULL
  })

  computed <- eventReactive(go(), {
    bad <- function(m) list(error = m)
    sp <- tryCatch(spec(), error = function(e) e)
    if (inherits(sp, "error")) return(bad(conditionMessage(sp)))
    ex <- tryCatch(existing(), error = function(e) e)
    if (inherits(ex, "error")) return(bad(conditionMessage(ex)))
    if (!is.null(ex)) ex$n1 <- as.numeric(input$n_new %||% NA)
    n <- if (is_exact()) as.numeric(input$n_new %||% NA) else NULL

    args <- tryCatch(
      owea:::.ui_solver_args(sp, if (is_exact()) "exact" else "optimal",
                             theta = theta(), p = as.integer(input$crit),
                             subset = subset_idx(), existing = ex, n = n),
      error = function(e) e)
    if (inherits(args, "error")) return(bad(conditionMessage(args)))
    if (is_exact()) args$seed <- as.integer(input$seed %||% 1)

    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch(
        withProgress(message = "Computing the design…", value = 0.5,
                     do.call(if (is_exact()) exact_design else optimal_design,
                             args)),
        error = function(e) structure(list(msg = conditionMessage(e)),
                                      class = "owea_bad")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w))
                              invokeRestart("muffleWarning") })
    if (inherits(res, "owea_bad")) return(bad(res$msg))
    list(res = res, warns = warns, cov_names = names(sp$design_box),
         exact = is_exact(), sp = sp, theta = theta(),
         p = as.integer(input$crit), subset = subset_idx(),
         existing = ex, args = args, coef_names = owea:::.ui_coef_names(sp),
         n = if (is_exact()) as.integer(n) else NA_integer_,
         n0 = if (is.null(ex)) 0L else as.integer(ex$n0))
  }, ignoreInit = TRUE)

  # a successful compute moves the wizard on; a failure keeps the user on Review
  observeEvent(computed(), {
    c <- computed()
    if (is.null(c$error)) goto("results")
    else showNotification(paste("Could not compute:", c$error), type = "error",
                          duration = 10)
  })

  # ---- status -------------------------------------------------------------
  output$status <- renderUI({
    c <- computed()
    if (!is.null(c$error))
      return(div(class = "alert alert-danger",
                 tags$b("Could not compute: "), c$error))
    conv <- isTRUE(c$res$converged)
    msg <- if (conv) "Design computed and verified optimal on the searched grid."
           else "Computed, but optimality was not fully certified — try a finer grid step."
    extra <- if (length(c$warns))
      tags$ul(lapply(unique(.friendly_warn(c$warns)), tags$li)) else NULL
    div(class = if (conv) "alert alert-success" else "alert alert-warning",
        tags$b(msg), extra)
  })

  # ---- design table -------------------------------------------------------
  design_df <- reactive({
    c <- computed(); req(is.null(c$error))
    S <- as.matrix(c$res$support); colnames(S) <- c$cov_names
    df <- as.data.frame(S)
    if (c$exact) df[["count"]] <- c$res$counts else df[["weight"]] <- c$res$weights
    df
  })
  output$design_tbl <- DT::renderDT({
    c <- computed(); req(is.null(c$error)); res <- c$res
    cap <- if (c$exact)
      sprintf("n = %d | criterion (per-sample) = %.5f | criterion (total) = %.5f | efficiency >= %.1f%%",
              res$n, res$criterion, res$criterion_total, 100 * res$efficiency)
    else
      sprintf("criterion = %.5f | max sensitivity = %.2e (0 at optimum)",
              res$criterion, res$max_d)
    disp <- design_df()
    for (nm in names(disp))
      if (is.numeric(disp[[nm]]) && !identical(nm, "count"))
        disp[[nm]] <- round(disp[[nm]], 4)
    DT::datatable(disp, rownames = FALSE, caption = cap,
                  options = list(dom = "t", paging = FALSE))
  })

  # the per-sample vs total criterion note (exact designs only)
  output$crit_note <- renderUI({
    c <- computed(); req(is.null(c$error)); if (!c$exact) return(NULL)
    res <- c$res
    N   <- c$n + c$n0                      # runs behind the TOTAL information
    arith <- if (c$p == 0L)
      sprintf("Here: %.5f − log(%d) = %.5f.", res$criterion, N, res$criterion_total)
    else
      sprintf("Here: %.5f / %d = %.5f.", res$criterion, N, res$criterion_total)
    div(class = "alert alert-info",
        tags$b("Why the two criteria do not simply differ by a factor of n."),
        tags$p(sprintf(paste0("The per-sample criterion is evaluated at the ",
                              "per-observation information matrix; the total one ",
                              "at the information from all N = %d run(s)%s. ",
                              "Multiplying the information by N divides the ",
                              "variances by N — but the two criteria absorb that ",
                              "differently:"),
                       N,
                       if (c$n0 > 0) sprintf(" (n0 = %d existing + n = %d new)",
                                             c$n0, c$n) else "")),
        tags$ul(
          tags$li(tags$b("A-optimality"), " averages variances, so the factor ",
                  "carries through: total = per-sample / N."),
          tags$li(tags$b("D-optimality"), " is a ", tags$i("log"),
                  "-determinant, so the factor becomes a shift: ",
                  "total = per-sample − log(N). Doubling N subtracts log 2 ≈ ",
                  "0.693; it does not halve the value.")),
        tags$p(arith),
        tags$small("Either number ranks designs of the same size identically — ",
                   "only the scale differs."))
  })

  output$info_help <- renderUI({
    c <- computed(); req(is.null(c$error))
    helpText(if (c$exact)
      paste("The TOTAL Fisher information from all runs (including any existing",
            "design) — the matrix behind the 'criterion (total)' value.")
      else
      paste("The per-observation Fisher information matrix (combined with any",
            "existing design)."))
  })

  output$design_plot <- renderPlot({
    c <- computed(); req(is.null(c$error))
    plot_design(c$res, cov_names = c$cov_names,
                main = if (c$exact) "Exact design" else "Approximate design")
  })
  output$info_tbl <- renderTable({
    c <- computed(); req(is.null(c$error))
    as.data.frame(round(c$res$information, 5))
  }, rownames = TRUE)
  output$model_txt <- renderPrint({
    c <- computed(); req(is.null(c$error))
    model_summary(c$res)
  })

  output$dl_csv <- downloadHandler(
    filename = function() "owea_design.csv",
    content  = function(file) {
      old <- options(digits = 15); on.exit(options(old))
      utils::write.csv(design_df(), file, row.names = FALSE)
    })
  output$dl_png <- downloadHandler(
    filename = function() "owea_design.png",
    content  = function(file) {
      grDevices::png(file, width = 900, height = 650, res = 110)
      on.exit(grDevices::dev.off())
      c <- computed()
      if (is.null(c$error)) plot_design(c$res, cov_names = c$cov_names)
    })

  # ---- efficiency under a DIFFERENT criterion -----------------------------
  output$eff_panel <- renderUI({
    c <- computed(); req(is.null(c$error))
    other <- if (identical(input$crit, "0")) "1" else "0"
    tagList(
      hr(),
      tags$h4("Efficiency under a different criterion"),
      helpText("How well does the design you just computed do under another ",
               "criterion, or for another set of parameters? The efficiency is ",
               "measured against the design that is optimal for THAT criterion, ",
               "so this runs a second optimisation."),
      selectInput("eff_crit", "Optimality criterion", choices = CRITERIA,
                  selected = other),
      radioButtons("eff_qoi", "Parameters of interest",
                   choices = c("All parameters" = "all", "A subset" = "subset"),
                   selected = if (is.null(c$subset)) "all" else "subset"),
      conditionalPanel("input.eff_qoi == 'subset'",
                       checkboxGroupInput("eff_subset", "Include parameters",
                                          choices = c$coef_names,
                                          selected = c$coef_names[c$subset])),
      actionButton("eff_btn", "Compute efficiency", class = "btn-info"),
      uiOutput("eff_out"))
  })

  observeEvent(input$eff_btn, {
    c <- computed(); req(is.null(c$error))
    sub <- if (identical(input$eff_qoi, "subset")) {
      sel <- input$eff_subset %||% character(0)
      if (length(sel)) match(sel, c$coef_names) else NULL
    } else NULL
    rv$eff <- withProgress(
      message = "Optimising under the other criterion…", value = 0.5,
      tryCatch(suppressWarnings(
        owea:::.ui_efficiency(c$res, c$sp, c$theta, as.integer(input$eff_crit),
                              sub, c$existing)),
        error = function(e) e))
  })

  output$eff_out <- renderUI({
    e <- rv$eff; if (is.null(e)) return(NULL)
    if (inherits(e, "error"))
      return(div(class = "alert alert-danger", tags$b("Could not compute: "),
                 conditionMessage(e)))
    lab <- names(CRITERIA)[match(as.character(e$p), CRITERIA)]
    if (!is.finite(e$efficiency))
      return(div(class = "alert alert-warning",
                 sprintf(paste0("The design is not estimable under %s (the ",
                                "criterion is infinite) — usually because those ",
                                "parameters are not identified by its support."),
                         lab)))
    tagList(br(), div(
      class = "alert alert-success",
      tags$p(tags$b(sprintf("Efficiency under %s: %.1f%%", lab,
                            100 * min(e$efficiency, 1)))),
      tags$p(sprintf("criterion of this design = %.6f; best achievable = %.6f",
                     e$crit_design, e$crit_ref)),
      if (!e$converged)
        tags$p(tags$small("The reference optimisation did not fully converge, ",
                          "so this is approximate.")),
      tags$small("100% would mean the design is also optimal for this criterion.")))
  })

  # ---- post-result: verify (approx) or simulate (exact) -------------------
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

  output$verify_panel <- renderUI({
    c <- computed(); req(is.null(c$error)); if (isTRUE(c$exact)) return(NULL)
    tagList(
      tags$h4("Verify optimality of a design"),
      helpText("Provide a design in the SAME format as the downloaded design CSV ",
               "(covariate columns + a 'weight' column). The box is prefilled with ",
               "the computed design — edit it, or upload a CSV, to verify a ",
               "different design. Optimality is checked under the ORIGINAL ",
               "criterion and parameters of interest, over the design box at ",
               "the finest grid step you specified."),
      fileInput("verify_file", "Upload design CSV (optional)", accept = ".csv"),
      textAreaInput("verify_manual", "… or paste / edit the design",
                    value = design_csv_text(design_df()), rows = 6),
      actionButton("verify_btn", "Verify optimality", class = "btn-info"),
      actionButton("reset_btn", "Reset"),
      uiOutput("verify_out"))
  })

  # the verify itself: original criterion (c$p) and parameters of interest
  # (c$subset), on the finest step of the original step sequence (via
  # .ui_solver_args' "verify" target)
  run_verify <- function(criterion_only = FALSE) {
    c <- computed(); if (is.null(c) || !is.null(c$error)) return()
    warns <- character(0)
    v <- withCallingHandlers(
      tryCatch({
        dz <- read_design(input$verify_file, input$verify_manual, c$cov_names,
                          valcol = "weight")
        w  <- as.numeric(dz$val); sw <- sum(w)
        if (!is.finite(sw) || sw <= 0)
          stop("weights must be positive and sum to a positive value.", call. = FALSE)
        if (abs(sw - 1) > 1e-6)
          warning(sprintf("Weights summed to %.6f; rescaled to sum to 1.", sw),
                  call. = FALSE)
        do.call(verify_optimality,
                c(list(support = dz$support, weights = w / sw,
                       criterion_only = criterion_only),
                  owea:::.ui_solver_args(c$sp, "verify", c$theta, c$p, c$subset,
                                         c$existing)))
      },
      error = function(e) structure(list(msg = conditionMessage(e)), class = "owea_bad")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w))
                              invokeRestart("muffleWarning") })
    rv$verify <- list(v = v, warns = warns, opt_crit = c$res$criterion)
  }

  # gate on the size of the finest-step grid: checking optimality scans the
  # WHOLE design space at the finest step, so past 1e6 points offer the same
  # three choices as verify_optimality() itself
  observeEvent(input$verify_btn, {
    c <- computed(); if (is.null(c) || !is.null(c$error)) return()
    N <- tryCatch(owea:::.ui_verify_points(c$sp), error = function(e) 0)
    if (is.finite(N) && N > 1e6)
      showModal(modalDialog(title = "Large design space",
        sprintf(paste0("Checking optimality evaluates the sensitivity at all ",
          "%s design points of the finest-step grid; building it may be slow."),
          format(round(N), big.mark = ",")),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("verify_full", "Proceed anyway (full check)"),
          actionButton("verify_crit_only",
                       "Criterion only (skip the optimality check)",
                       class = "btn-warning"))))
    else run_verify()
  })
  observeEvent(input$verify_full,      { removeModal(); run_verify(FALSE) })
  observeEvent(input$verify_crit_only, { removeModal(); run_verify(TRUE) })

  output$verify_out <- renderUI({
    vr <- rv$verify; if (is.null(vr)) return(NULL)
    if (inherits(vr$v, "owea_bad"))
      return(div(class = "alert alert-danger", tags$b("Verify failed: "), vr$v$msg))
    v <- vr$v
    warn_ui <- if (length(vr$warns))
      div(class = "alert alert-warning", tags$b("Warning: "),
          paste(unique(vr$warns), collapse = "  "))
    crit_ui <- tagList(
      tags$p(tags$b("Criterion (this design): "), sprintf("%.6f", v$criterion)),
      tags$p(tags$small(sprintf("Optimal design criterion (for reference): %.6f",
                                vr$opt_crit))))
    if (is.na(v$is_optimal$value))                  # criterion-only run
      return(tagList(br(), warn_ui,
        div(class = "alert alert-info", crit_ui,
            tags$p(tags$b("Optimality was NOT assessed"),
                   " (the max-sensitivity check was skipped)."),
            tags$p(tags$small(v$is_optimal$note)))))
    opt <- isTRUE(v$is_optimal$value)
    tagList(br(), warn_ui,
      div(class = if (opt) "alert alert-success" else "alert alert-warning",
          tags$p(tags$b("Max sensitivity: "),
                 sprintf("%.3e", v$max_sensitivity), " (0 at the optimum)"),
          crit_ui,
          tags$p(tags$b(if (opt) "The design IS optimal."
                        else "The design is NOT certified optimal.")),
          tags$p(tags$small(v$is_optimal$note))))
  })

  # ---- simulation study (exact design) ------------------------------------
  observeEvent(input$sim_open, { rv$show_sim <- TRUE })

  output$sim_panel <- renderUI({
    if (!isTRUE(rv$show_sim)) return(NULL)
    c <- computed(); req(is.null(c$error)); cn <- c$coef_names
    tagList(hr(), tags$h4("Simulation study"),
      helpText(sprintf("Total sample size n = %d (matched across designs).", c$n)),
      if (!is.null(c$existing))
        div(class = "alert alert-info",
            tags$b("Both stages are simulated and analysed together."),
            sprintf(paste0(" Each replicate estimates theta from all %d run(s) ",
                           "— the existing %d plus the new %d — because that is ",
                           "the analysis the design was optimised for."),
                    c$n0 + c$n, c$n0, c$n),
            if (!is.null(sim_obs()))
              tags$div(tags$small(paste(
                "The existing runs keep their OBSERVED responses (they are real",
                "data); only the new runs are simulated, from the true values",
                "below. Every design being compared shares that same observed",
                "half.")))
            else
              tags$div(tags$small(paste(
                "The existing design has no responses, so they are simulated",
                "from the true values below along with the new runs.")))),
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

  # Simulate a design and estimate theta the way the experiment actually will:
  # POOLED with the first stage, if there is one.  The design was optimised for
  # the combined information, so simulating the new runs alone would understate
  # every design's precision (and unequally, since the designs differ only in
  # their new half).
  sim_run <- function(support, counts) {
    c <- computed()
    theta_true <- vapply(seq_along(c$coef_names),
      function(i) as.numeric(input[[paste0("sim_theta_", i)]] %||% 0), numeric(1))
    sigma <- if (identical(c$sp$link, "identity"))
               as.numeric(input$sim_sigma %||% 1) else 1
    owea:::.ui_simulate(
      c$sp, theta = theta_true, sigma = sigma,
      support = support, counts = counts,
      # an OBSERVED first stage keeps its real responses; a design-only first
      # stage has none, so they are simulated too
      obs      = sim_obs(),
      existing = if (is.null(sim_obs())) c$existing else NULL,
      nsim = as.integer(input$sim_nsim %||% 1000),
      seed = as.integer(input$sim_seed %||% 1))
  }

  # the observed first stage, when the data set's covariates ARE the existing
  # design (otherwise the data set plays no part in the second stage)
  sim_obs <- reactive({
    c <- computed()
    if (is.null(c$existing) || !identical(input$start, "data") ||
        !identical(input$use_cov, "yes") || is.null(rv$data)) return(NULL)
    list(data = rv$data, response = input$data_response %||% NULL)
  })

  observeEvent(input$run_sim, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim
    s$exact <- tryCatch(sim_run(c$res$support, c$res$counts),
                        error = function(e) e)
    rv$sim <- s
  })
  observeEvent(input$run_srs, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim; s$pool_warn <- NULL
    r <- tryCatch({
      pool <- candidate_grid(c$sp$design_box, c$sp$finest)
      if (nrow(pool) > 1e6)
        s$pool_warn <- sprintf(paste0("The SRS candidate pool has %d points; ",
          "consider coarsening the step."), nrow(pool))
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
      dz  <- read_design(input$sim_custom_file, input$sim_custom, c$cov_names,
                         valcol = "count")
      cnt <- as.integer(dz$val)
      if (any(is.na(cnt))) stop("the 'count' column must be integers.", call. = FALSE)
      if (sum(cnt) != c$n)
        stop(sprintf("counts sum to %d but must sum to n = %d.", sum(cnt), c$n),
             call. = FALSE)
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

# a design data.frame -> CSV text (to prefill / display), full precision.
design_csv_text <- function(df) {
  old <- options(digits = 15); on.exit(options(old))
  paste(utils::capture.output(utils::write.csv(df, row.names = FALSE)), collapse = "\n")
}

# Read a table from an uploaded CSV (priority) or a pasted text box.
read_table <- function(file, text) {
  if (!is.null(file) && !is.null(file$datapath) && nzchar(file$datapath))
    return(utils::read.csv(file$datapath, header = TRUE, check.names = FALSE))
  if (!nzchar(trimws(text %||% "")))
    stop("provide the data: upload a CSV or paste one.", call. = FALSE)
  utils::read.csv(text = text, header = TRUE, check.names = FALSE)
}

# Read a design in the downloaded-design-CSV format: covariate columns plus one
# value column.  `valcol` names it; with valcol = NULL either a 'count' or a
# 'weight' column is accepted (counts win, as they also carry the sample size).
# Returns list(support = matrix, val = numeric, valcol = character).
read_design <- function(file, text, cov_names, valcol = NULL) {
  d <- read_table(file, text)
  if (is.null(valcol)) {
    valcol <- if ("count" %in% names(d)) "count"
              else if ("weight" %in% names(d)) "weight"
              else stop("the design needs a 'count' or a 'weight' column ",
                        "(see the downloaded CSV format).", call. = FALSE)
  } else if (!valcol %in% names(d)) {
    stop(sprintf("the design needs a '%s' column (see the downloaded CSV format).",
                 valcol), call. = FALSE)
  }
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
  list(support = sup, val = val, valcol = valcol)
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
