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
# lpSolve is an optional dependency of owea (the exact linear program behind
# Elfving's certificate of c-optimality, used when a single parameter is of
# interest).  Naming it here makes a deployment (shinyapps.io) install it; the
# package falls back to a derivative-free bound without it.
requireNamespace("lpSolve", quietly = TRUE)

MAX_COV   <- 6L
LINKS     <- c("Linear (normal)"        = "identity",
               "Logistic (binary)"      = "logit",
               "Poisson (counts)"       = "loglinear",
               "Multinomial (nominal)"  = "multinomial",
               "Ordinal (proportional odds)" = "cumulative")
MULTI_CAT <- c("multinomial", "cumulative")
CRITERIA  <- c("D-optimal (overall precision)" = "0",
               "A-optimal (average variance)"  = "1")
# the sample-size box means two different things, so it has two defaults: the
# TOTAL runs of an exact design, or the size of the new stage when augmenting
# an existing design with an approximate one.
N_DEFAULT_EXACT <- 200
N_DEFAULT_STAGE <- 20

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("owea — Optimal Experimental Designs"),
  tags$p("Build an optimal design step by step: describe your model (or models), ",
         "set the assumed parameter values and the criterion, then compute. ",
         "No R code required."),
  uiOutput("progress"),
  hr(),

  tabsetPanel(
    id = "wiz", type = "hidden",

    # ---- 0. mode ----------------------------------------------------------
    tabPanelBody(
      "mode",
      h4("Step: what kind of design problem is this?"),
      radioButtons(
        "mode_choice", NULL, selected = "classical",
        choiceValues = c("classical", "compound"),
        choiceNames = list(
          tags$span(tags$b("One criterion (classical)."),
                    " You have a single model and a single objective. The design ",
                    "is D-optimal (all parameters) or A-optimal (average ",
                    "variance), possibly for a subset of the parameters."),
          tags$span(tags$b("Several criteria at once (compound)."),
                    " You have more than one objective, and they can come from ",
                    "different models. One design is found that serves them all, ",
                    "weighted as you choose."))),
      div(class = "alert alert-info",
          tags$b("Which one do I want?"),
          tags$ul(
            tags$li(tags$b("Classical"), " — \"I know the model, and I want the ",
                    "most precise estimates of its parameters.\""),
            tags$li(tags$b("Compound"), " — \"I am not sure whether the ",
                    "interaction is real, and I want a design that is good ",
                    "either way\"; or \"I care about the slopes of one model ",
                    "and the overall fit of another\"; or \"my collaborators ",
                    "disagree about the model.\"")),
          tags$small("The compound branch reports how efficient the design is ",
                     "for each objective, and what a design optimised for any ",
                     "one objective alone would have cost under the others.")),
      div(class = "alert alert-warning",
          tags$b("Note."), " The compound branch builds approximate ",
          "(continuous-weight) and exact (integer-run) designs, and runs the ",
          "same simulation study. Fitting the parameters from a data set is ",
          "available in the classical branch only — fit them there, or ",
          "elsewhere, and type the estimates in as each objective's assumed ",
          "values.")
    ),

    # ---- 1. model ---------------------------------------------------------
    tabPanelBody(
      "model",
      h4("Step: the model"),
      div(class = "alert alert-info",
          tags$b("Not seeing your model? "),
          "The app covers the built-in model families and terms below, with ",
          "either all parameters or a subset of them as the parameters of interest. ",
          "Any other model, or another quantity of interest, can be handled in the ",
          "R package directly: pass your own per-point information as ",
          tags$code("info_vector"), " or ", tags$code("info_matrix"),
          ", and the parameters of interest as a matrix ", tags$code("wb"),
          " or a gradient function ", tags$code("grad_g"), ", to ",
          tags$code("optimal_design()"), " -- see the package README."),
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
          uiOutput("search_ui"),
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
                       numericInput("seed", "Random seed", value = 1, step = 1),
                       div(class = "alert alert-info",
                           tags$b("Tip: an exact design can be checked by simulation. "),
                           "Once it is computed, the results page offers a ",
                           tags$b("Simulation study"), " button. It simulates responses ",
                           "from your model at the design, fits the model to each ",
                           "simulated data set and reports the mean squared errors of the ",
                           "estimates, side by side with a simple random sample of the ",
                           "same size and, if you like, a design of your own -- a direct ",
                           "check of what the design delivers."))
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
      uiOutput("post_result"),
      # last, so it can include every step above (the design and its check)
      uiOutput("code_panel")
    ),

    # =====================================================================
    #  COMPOUND BRANCH
    #  One design, several criteria, possibly from different models.  The
    #  covariates (and hence the design region) are shared; each component
    #  carries its own model, criterion and parameters of interest.
    # =====================================================================

    # ---- C1. shared covariates -------------------------------------------
    tabPanelBody(
      "cmp_cov",
      h4("Step: the covariates"),
      helpText("One experiment, so one set of covariates and one design ",
               "region — shared by every objective. The models come next."),
      fluidRow(
        column(5,
               numericInput("cmp_ncov", "How many covariates?", value = 2,
                            min = 1, max = MAX_COV, step = 1)),
        column(7, uiOutput("cmp_cov_ui"))),
      uiOutput("cmp_cov_msg")
    ),

    # ---- C2. the components ----------------------------------------------
    tabPanelBody(
      "cmp_models",
      h4("Step: the objectives"),
      helpText("Each objective is a model plus a criterion. They may be ",
               "entirely different models — different families, different ",
               "terms, different numbers of parameters."),
      div(class = "alert alert-info",
          tags$b("Not seeing your model? "),
          "The app offers the built-in model families, with all parameters or a ",
          "subset as the parameters of interest. Any other model or quantity of ",
          "interest can be used in the R package directly: each component of ",
          tags$code("compound_design()"), " accepts its own ",
          tags$code("info_vector"), " or ", tags$code("info_matrix"),
          " and a matrix ", tags$code("wb"), " or function ", tags$code("grad_g"),
          " -- see the package README."),
      numericInput("ncomp", "How many objectives?", value = 2, min = 1,
                   max = 4, step = 1),
      uiOutput("cmp_models_ui"),
      uiOutput("cmp_models_msg")
    ),

    # ---- C3. options ------------------------------------------------------
    tabPanelBody(
      "cmp_options",
      h4("Step: weighting and any existing design"),
      radioButtons(
        "cmp_efficiency", "How should the weights be applied?",
        choiceValues = c("eff", "raw"), selected = "eff",
        choiceNames = list(
          tags$span(tags$b("To efficiencies (recommended)."),
                    " Each objective is divided by the best it could achieve on ",
                    "its own, so the weights compare like with like and the ",
                    "reported value is a weighted average of the objectives' ",
                    "efficiencies. It can never exceed 1, and the compound-optimal ",
                    "design reaches the highest value possible, which is below 1 ",
                    "unless one design is optimal for every objective at once."),
          tags$span(tags$b("To the raw criterion values."),
                    " Criteria of different models are on unrelated scales, so ",
                    "the weights will absorb that difference rather than express ",
                    "a preference."))),
      hr(),
      radioButtons(
        "cmp_design_type", "What kind of design?",
        choiceValues = c("approx", "exact"), selected = "approx",
        choiceNames = list(
          tags$span(tags$b("Approximate."),
                    " Continuous weights summing to 1 — the proportion of the ",
                    "effort each point should receive."),
          tags$span(tags$b("Exact."),
                    " An integer number of runs at each point, for a sample ",
                    "size you fix."))),
      conditionalPanel(
        "input.cmp_design_type == 'exact'",
        fluidRow(
          column(6, numericInput("cmp_n", "Number of runs (n)",
                                 value = N_DEFAULT_EXACT, min = 1, step = 1)),
          column(6, numericInput("cmp_seed", "Random seed", value = 1,
                                 step = 1))),
        helpText("The approximate optimum is computed first, rounded to whole ",
                 "runs, then improved by random exchanges."),
        div(class = "alert alert-info",
            tags$b("Tip: an exact design can be checked by simulation. "),
            "Once it is computed, the results page offers a ",
            tags$b("Simulation study"), " button. For every objective it simulates ",
            "responses from that objective's model at the design, fits the model to ",
            "each simulated data set and reports the mean squared errors of the ",
            "estimates -- so you can see, model by model, what one shared design ",
            "delivers.")),
      hr(),
      radioButtons("cmp_start", "Do you already have a design to build on?",
                   choices = c("No, start from scratch" = "none",
                               "Yes, augment an existing design" = "design"),
                   selected = "none"),
      conditionalPanel(
        "input.cmp_start == 'design'",
        helpText("Same format as the download: one column per covariate plus a ",
                 "'count' or 'weight' column."),
        fileInput("cmp_exist_file", "Upload a CSV", accept = ".csv"),
        textAreaInput("cmp_exist_text", "…or paste it here", rows = 5,
                      placeholder = "x1,x2,count\n-1,-1,10\n1,1,10"),
        fluidRow(
          column(6, numericInput("cmp_exist_n0", "Runs already made (n0)",
                                 value = NA, min = 1, step = 1)),
          column(6, numericInput("cmp_exist_n1", "New runs to add (n1)",
                                 value = NA, min = 1, step = 1))),
        uiOutput("cmp_exist_msg"))
    ),

    # ---- C4. review -------------------------------------------------------
    tabPanelBody(
      "cmp_review",
      h4("Step: review and compute"),
      uiOutput("cmp_review_out"),
      br(),
      actionButton("cmp_compute", "Compute design", class = "btn-primary",
                   width = "50%")
    ),

    # ---- C5. results ------------------------------------------------------
    tabPanelBody(
      "cmp_results",
      uiOutput("cmp_status"),
      tabsetPanel(
        tabPanel("Design",
                 br(), DT::DTOutput("cmp_design_tbl"),
                 br(), downloadButton("cmp_dl_csv", "Download design (CSV)")),
        tabPanel("Efficiency",
                 br(), uiOutput("cmp_eff_help"),
                 tableOutput("cmp_summary_tbl"),
                 br(), uiOutput("cmp_cross_help"),
                 tableOutput("cmp_cross_tbl")),
        tabPanel("Plot",
                 br(), plotOutput("cmp_plot", height = "460px"),
                 br(), downloadButton("cmp_dl_png", "Download plot (PNG)")),
        tabPanel("Information matrices",
                 br(), uiOutput("cmp_info_pick"),
                 tableOutput("cmp_info_tbl"))
      ),
      uiOutput("cmp_post_result")
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
    eff = NULL,
    cmp_verify = NULL, cmp_sim = list(),             # compound branch
    cmp_show_sim = FALSE)

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
            # the grid step(s) only matter when the continuous covariates are
            # searched on a grid; the continuous search needs none
            column(4, conditionalPanel(
              "input.search_mode != 'continuous'",
              textInput(paste0("cov_step_", i), "grid step(s)",
                        value = "0.1", placeholder = "e.g. 0.5, 0.1")))),
          conditionalPanel(
            "input.search_mode != 'continuous'",
            helpText("One number, or a comma-separated step sequence from coarse ",
                     "to fine (e.g. 0.5, 0.1, 0.02): the whole range is searched ",
                     "at the coarsest step, then refined locally at each finer one."))),
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
  has_continuous <- reactive(
    any(vapply(covariates(), function(cv) identical(cv$type, "continuous"), logical(1))))

  # how the continuous covariates are searched: on a grid of steps (the classic
  # path) or directly in the continuous region -- optimal_design(continuous =
  # TRUE), which needs no grid steps and does not grow with the number of
  # continuous covariates
  is_continuous <- reactive(has_continuous() &&
                            identical(input$search_mode %||% "grid", "continuous"))
  output$search_ui <- renderUI({
    if (!has_continuous()) return(NULL)
    tagList(
      radioButtons("search_mode", "How to search the continuous covariates",
                   choices = c("On a grid — enter the grid step(s) below" = "grid",
                               "Directly in the continuous region — no grid steps needed" = "continuous"),
                   selected = isolate(input$search_mode) %||% "grid"),
      conditionalPanel(
        "input.search_mode == 'continuous'",
        helpText("The continuous search starts from a coarse grid with three ",
                 "levels per continuous covariate, moves the support points to ",
                 "their best locations, and adds points where the sensitivity is ",
                 "largest, until a multi-start search of the region finds no ",
                 "violation of the equivalence theorem. Recommended when several ",
                 "covariates are continuous, where a fine grid grows too large. ",
                 "Once the design is computed, the results page lets you verify it ",
                 "on a grid or by a random audit of the region.")))
  })

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
                                   else "zero-sum",
                          continuous = is_continuous(),
                          # the search itself runs without a random audit; the
                          # audit is offered afterwards in the verify panel
                          n_audit = 0)
  })
  spec_ok    <- reactive(tryCatch({ spec(); NULL }, error = conditionMessage))
  coef_names <- reactive(tryCatch(owea:::.ui_coef_names(spec()),
                                  error = function(e) NULL))
  cov_names  <- reactive(tryCatch(names(spec()$design_box),
                                  error = function(e) NULL))

  # ---- navigation ---------------------------------------------------------
  # The wizard starts on the mode screen and then walks one of two chains.  The
  # classical chain is exactly .ui_wizard_steps() as before -- "mode" is
  # prepended HERE rather than inside that helper, so the helper (and the tests
  # that assert it verbatim) are untouched.
  cur  <- reactiveVal("mode")
  mode <- reactive(input$mode_choice %||% "classical")
  is_compound <- reactive(identical(mode(), "compound"))

  steps <- reactive({
    if (is_compound()) c("mode", owea:::.uic_wizard_steps())
    else c("mode", owea:::.ui_wizard_steps(input$link %||% "identity",
                                           input$start %||% "none"))
  })
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
    if (!cur() %in% c("review", "results", "cmp_review", "cmp_results"))
      actionButton("next_btn", "Next →", class = "btn-primary", width = "100%"))
  output$restart_ui <- renderUI(
    if (cur() %in% c("results", "cmp_results"))
      actionButton("restart_btn", "Start over", width = "40%"))

  goto <- function(id) { cur(id); updateTabsetPanel(session, "wiz", selected = id) }

  # per-step validation: a non-NULL message blocks Next
  step_error <- function(id) {
    if (identical(id, "mode")) return(NULL)
    # ---- compound branch ----
    if (identical(id, "cmp_cov"))
      return(tryCatch({ owea:::.ui_design_box(cmp_covariates()); NULL },
                      error = conditionMessage))
    if (identical(id, "cmp_models"))
      return(tryCatch(owea:::.uic_check(cmp_specs(), cmp_comps(), cmp_alpha()),
                      error = conditionMessage))
    if (identical(id, "cmp_options")) {
      if (!identical(input$cmp_start, "design")) return(NULL)
      return(tryCatch({ cmp_existing(); NULL }, error = conditionMessage))
    }
    # ---- classical branch ----
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
  grid_sizes <- reactive(tryCatch(owea:::.ui_grid_sizes(covariates(), is_continuous()),
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
  # on (e.g. switching to the identity link while on the theta step); switching
  # branch replaces the chain wholesale
  observeEvent(list(input$link, input$start, input$mode_choice, input$cmp_start), {
    if (!cur() %in% steps())
      goto(if (is_compound()) "cmp_cov" else "model")
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
    # this box is re-rendered whenever the design type flips, so a plain
    # isolate(input$n_new) would carry the approximate default into the exact
    # branch: swap between the two defaults, but never discard a typed n.
    def  <- if (is_exact()) N_DEFAULT_EXACT else N_DEFAULT_STAGE
    prev <- isolate(input$n_new)
    val  <- if (is.null(prev) || !is.finite(prev) ||
                prev %in% c(N_DEFAULT_STAGE, N_DEFAULT_EXACT)) def else prev
    tagList(
      numericInput("n_new", lab, value = val, min = 1, step = 1),
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
      if (isTRUE(sp$continuous))
        sprintf(paste0("Design space: the continuous region is searched directly ",
                       "(no grid; the search starts from a coarse %s-point grid)"),
                format(round(gs[1]), big.mark = ","))
      else if (length(gs) && is.finite(gs[1]))
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
    # the continuous search certifies over the points it examined (multi-start
    # search + random audit), the grid path over the searched grid
    cont <- isTRUE(c$sp$continuous)
    n_au <- if (cont && !is.null(c$res$n_audit)) c$res$n_audit else 0L
    msg <- if (conv && cont)
      sprintf(paste0("Design computed by the continuous search: a multi-start search of ",
                     "the region%s found no violation of the equivalence theorem. ",
                     "Verify it below on a grid or by a random audit."),
              if (n_au > 0) sprintf(" plus a random audit of %s points",
                                    format(n_au, big.mark = ","))
              else "")
    else if (conv) "Design computed and verified optimal on the searched grid."
    else if (cont)
      paste0("Computed, but optimality was not fully certified — the continuous ",
             "search stopped short. In the R package, raise n_starts or max_iter, ",
             "or supply init_points.")
    else "Computed, but optimality was not fully certified — try a finer grid step."
    extra <- if (length(c$warns))
      tags$ul(lapply(unique(.friendly_warn(c$warns)), tags$li)) else NULL
    # certified efficiency bound (Becker & Yang, Theorems 4.5 / 4.6): valid
    # even when the design did not fully converge
    eb <- if (c$exact) NULL else c$res$efficiency_lower_bound
    cert <- if (!is.null(eb) && is.finite(eb))
      tags$p(tags$small(sprintf(paste0(
        "Efficiency: at least %.4f%% of the optimum over %s ",
        "(a guaranteed lower bound computed from the max sensitivity, valid even ",
        "when convergence was not reached)."), 100 * eb,
        if (cont) "the points examined by the search and the audit"
        else "the searched grid")))
    div(class = if (conv) "alert alert-success" else "alert alert-warning",
        tags$b(msg), cert, extra)
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
      sprintf("n = %d | criterion (per-sample) = %.5f | criterion (total) = %.5f | efficiency >= %.4f%%",
              res$n, res$criterion, res$criterion_total, 100 * res$efficiency_lower_bound)
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
  # ---- the R code behind the result (rendered LAST on the results page) ----
  # the computation the app made and, for an approximate design, the
  # equivalence-theorem check at the grid step entered in the verify panel
  # for an exact design whose simulation study has been run: the study's
  # settings and the designs it compared, so the script can repeat it with
  # simulate_design()
  sim_code_info <- reactive({
    c <- computed(); s <- rv$sim
    if (!isTRUE(c$exact) || is.null(s$exact) || inherits(s$exact, "error")) return(NULL)
    theta_true <- vapply(seq_along(c$coef_names),
      function(i) as.numeric(input[[paste0("sim_theta_", i)]] %||% 0), numeric(1))
    ex <- if (is.null(sim_obs())) c$existing else NULL
    keep <- function(res, dsg) if (!is.null(res) && !inherits(res, "error")) dsg else NULL
    list(theta = theta_true,
         sigma = if (identical(c$sp$link, "identity")) as.numeric(input$sim_sigma %||% 1) else NULL,
         nsim  = as.integer(input$sim_nsim %||% 1000),
         seed  = as.integer(input$sim_seed %||% 1),
         existing = if (!is.null(ex))
           list(points = ex$points,
                counts = ex$counts %||% owea:::.apportion(ex$weights, as.integer(ex$n0)))
           else NULL,
         obs = sim_obs(),
         designs = Filter(Negate(is.null),
                          list(SRS = keep(s$srs, s$srs_design),
                               Custom = keep(s$custom, s$custom_design))))
  })

  r_code <- reactive({
    c <- computed(); req(is.null(c$error))
    vargs <- if (c$exact) NULL else {
      has_cont <- any(!owea:::.parse_design_box(c$sp$design_box)$is_factor)
      if (has_cont && verify_audit_mode()) {
        na <- verify_audit_eff()
        tryCatch(owea:::.ui_solver_args(c$sp, "verify", c$theta, c$p, c$subset,
                                        c$existing,
                                        verify_audit = if (is.na(na)) 20000L else na),
                 error = function(e) NULL)
      } else {
        se <- verify_step_eff(c)
        tryCatch(owea:::.ui_solver_args(c$sp, "verify", c$theta, c$p, c$subset,
                                        c$existing,
                                        verify_step = if (se$ok) se$step else NULL),
                 error = function(e) NULL)
      }
    }
    owea:::.ui_r_code(c$args, if (c$exact) "exact" else "optimal", verify_args = vargs,
                      sim = sim_code_info())
  })
  output$code_txt <- renderText(r_code())
  output$dl_code <- downloadHandler(
    filename = function() "owea_design.R",
    content  = function(file) writeLines(r_code(), file))
  output$code_panel <- renderUI({
    c <- computed(); req(is.null(c$error))
    has_sim <- !is.null(sim_code_info())
    tagList(
      hr(),
      tags$h4("R code for this analysis"),
      helpText("The exact call the app made",
               if (!c$exact) ", followed by the optimality check as set up in the verify section above (grid or random audit)",
               if (has_sim) ", followed by the simulation study as run above (with the designs it compared)",
               ". Copy it, or download it, to rerun and adapt the computation in R."),
      verbatimTextOutput("code_txt"),
      downloadButton("dl_code", "Download R code (.R)"))
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
    if (!is.finite(e$efficiency_lower_bound))
      return(div(class = "alert alert-warning",
                 sprintf(paste0("The design is not estimable under %s (the ",
                                "criterion is infinite) — usually because those ",
                                "parameters are not identified by its support."),
                         lab)))
    tagList(br(), div(
      class = "alert alert-success",
      tags$p(tags$b(sprintf("Efficiency under %s: at least %.4f%%", lab,
                            100 * e$efficiency_lower_bound))),
      tags$p(sprintf("criterion of this design = %.6f; optimal criterion = %.6f",
                     e$crit_design, e$crit_ref)),
      tags$p(tags$small("The efficiency is a guaranteed lower bound relative to the ",
                        "true optimum for this criterion.")),
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
    has_cont <- any(!owea:::.parse_design_box(c$sp$design_box)$is_factor)
    cont <- isTRUE(c$sp$continuous)
    tagList(
      tags$h4("Verify optimality of a design"),
      helpText("Provide a design in the SAME format as the downloaded design CSV ",
               "(covariate columns + a 'weight' column). The box is prefilled with ",
               "the computed design — edit it, or upload a CSV, to verify a ",
               "different design. Optimality is checked under the ORIGINAL ",
               "criterion and parameters of interest (the general equivalence ",
               "theorem)",
               if (has_cont) ", over a grid or by a random audit of the region."
               else " at every level combination of the factors."),
      fileInput("verify_file", "Upload design CSV (optional)", accept = ".csv"),
      textAreaInput("verify_manual", "… or paste / edit the design",
                    value = design_csv_text(design_df()), rows = 6),
      # how to check: a grid at a step of the user's choice, or a random audit
      # of the continuous region with a chosen number of points
      if (has_cont) radioButtons(
        "verify_mode", "Check the equivalence theorem",
        choices = c("On a grid over the design box (step(s) below)" = "grid",
                    "By a random audit of the design region (number of points below)" = "audit"),
        selected = isolate(input$verify_mode) %||% "grid"),
      if (has_cont) conditionalPanel(
        "input.verify_mode != 'audit'",
        textInput("verify_step", "grid step(s) for the check",
                  value = isolate(input$verify_step) %||% owea:::.ui_verify_step_default(c$sp),
                  placeholder = "one number, or one per continuous covariate"),
        if (!cont)
          div(class = "alert alert-info",
              sprintf(paste0("Optimality has already been verified on the grid of the ",
                             "computation (finest step %s). Choose a different step ",
                             "here to avoid repeating that computation."),
                      owea:::.ui_verify_step_default(c$sp)))
        else
          helpText("The design was found by a continuous search, so its support ",
                   "points are off the grid: the sensitivity is evaluated at every ",
                   "grid point with the design taken exactly as given."),
        helpText("One number, or one per continuous covariate. A finer step is a ",
                 "stronger check but a larger grid; beyond 1,000,000 points you ",
                 "are asked before it is built.")),
      if (has_cont) conditionalPanel(
        "input.verify_mode == 'audit'",
        numericInput("verify_audit", "Random audit points",
                     value = isolate(input$verify_audit) %||% 20000, min = 1, step = 1000),
        helpText("The sensitivity is evaluated at this many random points of the ",
                 "region, the best few are polished locally, and a multi-start ",
                 "search from the support points and the box corners is added. Each ",
                 "point costs one model evaluation (about a second per million ",
                 "points for the built-in models).")),
      actionButton("verify_btn", "Verify optimality", class = "btn-info"),
      actionButton("reset_btn", "Reset"),
      uiOutput("verify_out"))
  })

  # the audit mode of the verify panel, and its number of points (NA when
  # something unusable was typed)
  verify_audit_mode <- reactive(identical(input$verify_mode %||% "grid", "audit"))
  verify_audit_eff  <- function() {
    n <- suppressWarnings(as.integer(round(as.numeric(input$verify_audit %||% 20000))))
    if (length(n) != 1L || is.na(n) || n < 1L) NA_integer_ else n
  }

  # the grid step(s) of the check: what was typed, or the suggested default when
  # the box is empty (NULL for an all-factor box, which is enumerated); ok is
  # FALSE only when something unusable was typed
  verify_step_eff <- function(c) {
    if (!any(!owea:::.parse_design_box(c$sp$design_box)$is_factor))
      return(list(step = NULL, ok = TRUE))
    raw <- trimws(as.character(input$verify_step %||% ""))
    if (!nzchar(raw)) raw <- owea:::.ui_verify_step_default(c$sp)
    st <- owea:::.ui_parse_steps(raw)
    list(step = st, ok = length(st) > 0 && all(st > 0))
  }

  # a one-line description of the design space a check runs over
  verify_space_label <- function(st, N, audit = NA) {
    if (!is.na(audit))
      return(sprintf(paste0("a random audit of the design region with %s points (plus ",
                            "a multi-start search from the support points and the box ",
                            "corners)"), format(audit, big.mark = ",")))
    npts <- if (is.finite(N)) format(round(N), big.mark = ",") else "?"
    if (is.null(st))
      sprintf("all level combinations of the factors (%s design points)", npts)
    else
      sprintf("a grid over the design box with step(s) %s (%s design points)",
              paste(format(st), collapse = ", "), npts)
  }

  # the verify itself: original criterion (c$p) and parameters of interest
  # (c$subset), over the grid at the entered step or by the random audit (via
  # .ui_solver_args' "verify" target with verify_step / verify_audit)
  run_verify <- function(criterion_only = FALSE) {
    c <- computed(); if (is.null(c) || !is.null(c$error)) return()
    has_cont <- any(!owea:::.parse_design_box(c$sp$design_box)$is_factor)
    bad <- function(msg) {
      rv$verify <- list(v = structure(list(msg = msg), class = "owea_bad"),
                        warns = character(0), opt_crit = c$res$criterion)
    }
    audit <- NA_integer_; st <- NULL; N <- NA_real_
    if (has_cont && verify_audit_mode()) {
      audit <- verify_audit_eff()
      if (is.na(audit)) return(bad("enter a positive number of random audit points."))
    } else {
      se <- verify_step_eff(c)
      if (!se$ok)
        return(bad("enter positive grid step(s) for the check: one number, or one per continuous covariate."))
      st <- se$step
      N <- tryCatch(owea:::.ui_verify_points(c$sp, step = if (is.null(st)) 1 else st),
                    error = function(e) NA_real_)
    }
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
                                         c$existing, verify_step = st,
                                         verify_audit = if (is.na(audit)) NULL else audit)))
      },
      error = function(e) structure(list(msg = conditionMessage(e)), class = "owea_bad")),
      warning = function(w) { warns <<- c(warns, conditionMessage(w))
                              invokeRestart("muffleWarning") })
    rv$verify <- list(v = v, warns = warns, opt_crit = c$res$criterion,
                      space = verify_space_label(st, N, audit))
  }

  # gate on the size of the grid to scan (grid checks only): past 1e6 points
  # offer the same three choices as verify_optimality() itself
  observeEvent(input$verify_btn, {
    c <- computed(); if (is.null(c) || !is.null(c$error)) return()
    has_cont <- any(!owea:::.parse_design_box(c$sp$design_box)$is_factor)
    if (has_cont && verify_audit_mode()) { run_verify(); return() }
    se <- verify_step_eff(c)
    N <- if (se$ok)
      tryCatch(owea:::.ui_verify_points(c$sp, step = if (is.null(se$step)) 1 else se$step),
               error = function(e) 0)
    else 0
    if (is.finite(N) && N > 1e6)
      showModal(modalDialog(title = "Large design space",
        sprintf(paste0("Checking optimality evaluates the sensitivity at all ",
          "%s design points of this grid; building it may be slow."),
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
      div(class = "alert alert-warning", tags$b("Note: "),
          paste(unique(.friendly_warn(vr$warns)), collapse = "  "))
    space_ui <- if (!is.null(vr$space))
      tags$p(tags$b("Checked over: "), vr$space)
    crit_ui <- tagList(
      space_ui,
      tags$p(tags$b("Criterion (this design): "), sprintf("%.6f", v$criterion)),
      tags$p(tags$small(sprintf("Optimal design criterion (for reference): %.6f",
                                vr$opt_crit))),
      if (!is.null(v$efficiency_lower_bound) && is.finite(v$efficiency_lower_bound))
        tags$p(tags$b("Efficiency: "),
               sprintf(paste0("at least %.4f%% of the true optimum over this design ",
                              "space (a guaranteed lower bound from the max sensitivity)"),
                       100 * v$efficiency_lower_bound)))
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
      br(), br(), tableOutput("sim_tbl"),
      uiOutput("sim_note"), uiOutput("sim_msg"))
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
  # the simple random sample: n runs from the finest-step grid after a grid
  # computation, or uniformly from the continuous region after a continuous
  # search (which has no grid) -- see .ui_srs_design()
  observeEvent(input$run_srs, {
    c <- computed(); req(is.null(c$error))
    s <- rv$sim; s$pool_warn <- NULL
    r <- tryCatch({
      d <- owea:::.ui_srs_design(c$sp, c$n)
      if (is.finite(d$pool_size) && d$pool_size > 1e6)
        s$pool_warn <- sprintf(paste0("The SRS candidate pool has %d points; ",
          "consider coarsening the step."), d$pool_size)
      s$srs_design <- list(support = d$support, counts = d$counts)
      sim_run(d$support, d$counts)
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
      s$custom_design <- list(support = dz$support, counts = cnt)
      sim_run(dz$support, cnt)
    }, error = function(e) e)
    s$custom <- r; rv$sim <- s
  })

  # The comparison, assembled ONCE, so the error columns and the ratio columns
  # derived from them can never disagree about which designs were compared.
  # Mean AND median squared error: the mean is the usual summary, the median
  # says what a typical replicate achieved when a tail is distorting the mean.
  sim_mats <- reactive({
    c <- computed(); req(is.null(c$error)); s <- rv$sim
    slots <- list(`Exact design` = s$exact, SRS = s$srs, Custom = s$custom)
    grab <- function(what) owea:::.ui_sim_matrix(
      lapply(slots, function(r)
        if (is.null(r) || inherits(r, "error")) NULL else r[[what]]),
      param_names = c$coef_names)
    list(`Mean squared error` = grab("mse"),
         `Median squared error` = grab("medse"))
  })

  # every design compared here is measured against the one the app computed
  SIM_REF <- "Exact design"

  output$sim_tbl <- renderTable({
    m <- sim_mats(); M <- m[["Mean squared error"]]
    if (is.null(M)) return(NULL)
    df <- data.frame(parameter = colnames(M), stringsAsFactors = FALSE)
    # each measure: one column per design, then that design's ratio to the
    # computed design -- above 1, the computed design has the smaller error
    add <- function(df, A, unit) {
      if (is.null(A)) return(df)
      for (nm in rownames(A)) df[[sprintf("%s (%s)", nm, unit)]] <- A[nm, ]
      R <- owea:::.ui_sim_ratio(A, SIM_REF)
      if (!is.null(R))
        for (nm in rownames(R))
          df[[sprintf("%s / %s (%s)", nm, SIM_REF, unit)]] <- R[nm, ]
      df
    }
    df <- add(df, M, "MSE")
    add(df, m[["Median squared error"]], "median SE")
  }, digits = 6)

  # How to read the ratio columns, and a warning when a single replicate is
  # carrying an MSE -- see .ui_sim_dominance().
  output$sim_note <- renderUI({
    s <- rv$sim
    if (!length(s)) return(NULL)
    ratio_ui <- if (!is.null(owea:::.ui_sim_ratio(sim_mats()[[1]], SIM_REF)))
      helpText(tags$b(sprintf("The \"/ %s\" columns are ratios.", SIM_REF)),
               paste("Above 1 the computed design has the smaller error for",
                     "that parameter (2.0 = twice the error); below 1 the other",
                     "design does. The errors themselves are on each",
                     "parameter's own scale, so only the ratio compares",
                     "designs across a row."))
    lab <- c(exact = "the exact design", srs = "SRS", custom = "custom")
    conv <- character(0); wild <- character(0)
    for (nm in c("exact", "srs", "custom")) {
      r <- s[[nm]]
      if (is.null(r) || inherits(r, "error")) next
      conv <- c(conv, sprintf("%s %d/%d", lab[[nm]], r$n_converged, r$nsim))
      d <- owea:::.ui_sim_dominance(r)
      if (is.finite(d) && d > 0.2)
        wild <- c(wild, sprintf("%s (one replicate contributes %.0f%% of an MSE)",
                                lab[[nm]], 100 * d))
    }
    tagList(
      ratio_ui,
      if (length(conv))
        helpText(sprintf("Replications that converged: %s.",
                         paste(conv, collapse = ", "))),
      if (length(wild))
        div(class = "alert alert-warning",
            tags$b("Read the MSE column with care."),
            sprintf(" Some fits returned extreme estimates — %s. ",
                    paste(wild, collapse = ", ")),
            "That is separation: a design concentrated on few distinct points ",
            "can, in a given replicate, leave a binary or count response ",
            "perfectly predicted, and the fitted coefficients then run away. ",
            "One such replicate can dominate a mean squared error, so a design ",
            "that is better on average can look worse here. ",
            tags$small("The median-squared-error columns are not affected by a ",
                       "handful of wild replicates; compare those, and the ",
                       "criterion on the Design tab.")))
  })

  output$sim_msg <- renderUI({
    s <- rv$sim; msgs <- character(0)
    for (nm in c("exact", "srs", "custom"))
      if (inherits(s[[nm]], "error"))
        msgs <- c(msgs, sprintf("%s: %s", nm, conditionMessage(s[[nm]])))
    if (!is.null(s$pool_warn)) msgs <- c(msgs, s$pool_warn)
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-warning", lapply(msgs, tags$div))
  })

  # =======================================================================
  #  COMPOUND BRANCH
  #  Everything below drives compound_design()/compound_criterion() through
  #  owea:::.uic_solver_args(), the compound counterpart of the invariant the
  #  classical branch keeps with .ui_solver_args(): one assembler, so the
  #  models cannot drift between the design and the verification.
  # =======================================================================

  # memoised reference values: list(key = .uic_psi_key(args), psi_star = ...)
  # from the last successful solve.  Deliberately a plain local, not a
  # reactiveValue -- nothing renders it, and no output should invalidate when
  # it changes.
  cmp_psi_cache <- NULL

  # ---- C1. the shared covariates ------------------------------------------
  output$cmp_cov_ui <- renderUI({
    n <- max(1L, min(MAX_COV, as.integer(input$cmp_ncov %||% 2)))
    lapply(seq_len(n), function(i) {
      wellPanel(
        style = "padding:8px;",
        textInput(paste0("ccov_name_", i), NULL, value = paste0("x", i),
                  placeholder = "covariate name"),
        selectInput(paste0("ccov_type_", i), NULL,
                    choices = c("Continuous" = "continuous", "Factor" = "factor")),
        conditionalPanel(
          sprintf("input.ccov_type_%d == 'continuous'", i),
          fluidRow(
            column(4, numericInput(paste0("ccov_lo_", i), "low",  value = -1)),
            column(4, numericInput(paste0("ccov_hi_", i), "high", value = 1)),
            column(4, textInput(paste0("ccov_step_", i), "grid step(s)",
                                value = "0.1", placeholder = "e.g. 0.5, 0.1")))),
        conditionalPanel(
          sprintf("input.ccov_type_%d == 'factor'", i),
          numericInput(paste0("ccov_nlev_", i), "number of levels", value = 2,
                       min = 2, step = 1)))
    })
  })

  cmp_covariates <- reactive({
    n <- max(1L, min(MAX_COV, as.integer(input$cmp_ncov %||% 2)))
    lapply(seq_len(n), function(i) {
      type <- input[[paste0("ccov_type_", i)]] %||% "continuous"
      nm   <- input[[paste0("ccov_name_", i)]] %||% paste0("x", i)
      if (identical(type, "factor")) {
        list(name = nm, type = "factor",
             nlevels = as.integer(input[[paste0("ccov_nlev_", i)]] %||% 2))
      } else {
        list(name = nm, type = "continuous",
             lo    = as.numeric(input[[paste0("ccov_lo_", i)]] %||% -1),
             hi    = as.numeric(input[[paste0("ccov_hi_", i)]] %||%  1),
             steps = owea:::.ui_parse_steps(input[[paste0("ccov_step_", i)]] %||% "0.1"))
      }
    })
  })

  cmp_has_factor <- reactive(
    any(vapply(cmp_covariates(), function(cv) identical(cv$type, "factor"),
               logical(1))))
  cmp_grid_sizes <- reactive(tryCatch(owea:::.ui_grid_sizes(cmp_covariates()),
                                      error = function(e) NULL))

  output$cmp_cov_msg <- renderUI({
    err <- tryCatch({ owea:::.ui_design_box(cmp_covariates()); NULL },
                    error = conditionMessage)
    if (!is.null(err))
      return(div(class = "alert alert-danger", err))
    gs <- cmp_grid_sizes()
    if (is.null(gs) || !length(gs)) return(NULL)
    div(class = if (gs[1] > 1e6) "alert alert-warning" else "alert alert-info",
        sprintf("Candidate points at the coarsest stage: %s.",
                format(round(gs[1]), big.mark = ",")),
        if (gs[1] > 1e6)
          tags$span(" That is a very large grid — consider a coarser first step ",
                    "or a step sequence such as \"0.5, 0.1\"."))
  })

  # ---- C2. the components --------------------------------------------------
  output$cmp_models_ui <- renderUI({
    J  <- max(1L, min(4L, as.integer(input$ncomp %||% 2)))
    cv <- cmp_covariates(); n <- length(cv)
    pairs <- if (n >= 2) combn(n, 2, simplify = FALSE) else list()
    cont  <- which(vapply(cv, function(z) identical(z$type, "continuous"),
                          logical(1)))
    lapply(seq_len(J), function(j) {
      wellPanel(
        tags$b(sprintf("Objective %d", j)),
        fluidRow(
          column(4, textInput(paste0("cmp_name_", j), "Name",
                              value = paste0("objective ", j))),
          column(4, selectInput(paste0("cmp_link_", j), "Model family",
                                choices = LINKS)),
          column(4, numericInput(paste0("cmp_alpha_", j), "Weight", value = 1,
                                 min = 0, step = 0.05))),
        conditionalPanel(
          sprintf("input.cmp_link_%d == 'multinomial' || input.cmp_link_%d == 'cumulative'",
                  j, j),
          numericInput(paste0("cmp_ncat_", j), "Number of response categories",
                       value = 3, min = 2, step = 1)),
        fluidRow(
          column(6,
                 if (length(pairs))
                   checkboxGroupInput(
                     paste0("cmp_int_", j), "Interactions",
                     choices = setNames(
                       vapply(pairs, function(p) paste(p, collapse = "-"),
                              character(1)),
                       vapply(pairs, function(p)
                         sprintf("%s × %s", cv[[p[1]]]$name, cv[[p[2]]]$name),
                         character(1))))),
          column(6,
                 if (length(cont))
                   checkboxGroupInput(
                     paste0("cmp_quad_", j), "Quadratic terms",
                     choices = setNames(as.character(cont),
                                        vapply(cont, function(i) cv[[i]]$name,
                                               character(1)))))),
        fluidRow(
          column(6, selectInput(paste0("cmp_crit_", j), "Criterion",
                                choices = CRITERIA)),
          column(6, radioButtons(paste0("cmp_qoi_", j), "Parameters of interest",
                                 choices = c("All" = "all", "A subset" = "subset"),
                                 inline = TRUE))),
        conditionalPanel(sprintf("input.cmp_qoi_%d == 'subset'", j),
                         uiOutput(paste0("cmp_subset_ui_", j))),
        uiOutput(paste0("cmp_theta_ui_", j)))
    })
  })

  cmp_J <- reactive(max(1L, min(4L, as.integer(input$ncomp %||% 2))))

  # The MODEL-defining half of a component: link, family size, terms, coding.
  # Deliberately does NOT read the theta boxes or the subset checkboxes -- the
  # spec does not depend on them, and reading them here would make the widgets
  # that SET them depend on themselves: ticking a subset box would invalidate
  # the spec, re-render the checkbox group, and wipe the tick.
  cmp_model_comp <- function(j) {
    lk <- input[[paste0("cmp_link_", j)]] %||% "identity"
    list(
      link  = lk,
      ncat  = if (lk %in% MULTI_CAT)
                as.integer(input[[paste0("cmp_ncat_", j)]] %||% 3) else NULL,
      coding = if (cmp_has_factor()) input$cmp_coding %||% "zero-sum" else "zero-sum",
      interactions = lapply(input[[paste0("cmp_int_", j)]] %||% character(0),
                            function(s) as.integer(strsplit(s, "-", fixed = TRUE)[[1]])),
      quadratics = as.integer(input[[paste0("cmp_quad_", j)]] %||% character(0)))
  }
  # the specs depend on the models only, so the theta/subset widgets can safely
  # be built from them
  cmp_specs <- reactive(
    owea:::.uic_specs(cmp_covariates(), lapply(seq_len(cmp_J()), cmp_model_comp)))

  cmp_coef_names <- function(j)
    tryCatch(owea:::.ui_coef_names(cmp_specs()[[j]]), error = function(e) NULL)

  # the FULL component: the model plus the criterion, the parameters of
  # interest and the assumed values.  Used for validation and for solving --
  # never for rendering the inputs it reads.
  cmp_comp <- function(j) {
    m  <- cmp_model_comp(j)
    cn <- cmp_coef_names(j)
    k  <- length(cn %||% character(0))
    c(m, list(
      name = input[[paste0("cmp_name_", j)]] %||% sprintf("objective %d", j),
      p = as.integer(input[[paste0("cmp_crit_", j)]] %||% "0"),
      subset = if (identical(input[[paste0("cmp_qoi_", j)]] %||% "all", "subset"))
                 as.integer(input[[paste0("cmp_subset_", j)]] %||% integer(0))
               else NULL,
      theta = if (identical(m$link, "identity") || !k) NULL
              else vapply(seq_len(k), function(i)
                as.numeric(input[[sprintf("cmp_theta_%d_%d", j, i)]] %||% 0),
                numeric(1))))
  }
  cmp_comps <- reactive(lapply(seq_len(cmp_J()), cmp_comp))
  cmp_alpha <- reactive(vapply(seq_len(cmp_J()), function(j)
    as.numeric(input[[paste0("cmp_alpha_", j)]] %||% 1), numeric(1)))

  # per-component theta boxes and subset pickers.  Both are built from
  # cmp_specs() only; the current selection is carried over with isolate() so a
  # re-render caused by an unrelated change (adding an interaction, say) does
  # not silently discard what the user typed or ticked.
  observe({
    lapply(seq_len(cmp_J()), function(j) {
      local({
        jj <- j
        output[[paste0("cmp_theta_ui_", jj)]] <- renderUI({
          sp <- tryCatch(cmp_specs()[[jj]], error = function(e) NULL)
          if (is.null(sp) || identical(sp$link, "identity")) return(NULL)
          cn <- cmp_coef_names(jj)
          if (is.null(cn)) return(NULL)
          tagList(
            tags$b("Assumed parameter values"),
            helpText("This family gives a locally optimal design, so it needs ",
                     "assumed values."),
            actionButton(paste0("cmp_draw_", jj), "Draw from N(0,1)"),
            br(), br(),
            fluidRow(lapply(seq_along(cn), function(i)
              column(3, numericInput(
                sprintf("cmp_theta_%d_%d", jj, i), cn[i], step = 0.1,
                value = isolate(input[[sprintf("cmp_theta_%d_%d", jj, i)]]) %||% 0)))))
        })
        output[[paste0("cmp_subset_ui_", jj)]] <- renderUI({
          cn <- cmp_coef_names(jj)
          if (is.null(cn)) return(NULL)
          sel <- isolate(input[[paste0("cmp_subset_", jj)]])
          checkboxGroupInput(
            paste0("cmp_subset_", jj), NULL, inline = TRUE,
            choices  = setNames(as.character(seq_along(cn)), cn),
            selected = sel[sel %in% as.character(seq_along(cn))])
        })
      })
    })
  })

  # "draw from N(0,1)" for each component
  observe({
    J <- max(1L, min(4L, as.integer(input$ncomp %||% 2)))
    lapply(seq_len(J), function(j) {
      local({
        jj <- j
        observeEvent(input[[paste0("cmp_draw_", jj)]], {
          sp <- tryCatch(cmp_specs()[[jj]], error = function(e) NULL)
          if (is.null(sp)) return()
          th <- tryCatch(owea:::.ui_random_theta(sp), error = function(e) NULL)
          if (is.null(th)) return()
          for (i in seq_along(th))
            updateNumericInput(session, sprintf("cmp_theta_%d_%d", jj, i),
                               value = round(th[i], 4))
        }, ignoreInit = TRUE)
      })
    })
  })

  output$cmp_models_msg <- renderUI({
    msg <- tryCatch(owea:::.uic_check(cmp_specs(), cmp_comps(), cmp_alpha()),
                    error = conditionMessage)
    if (is.null(msg)) return(NULL)
    div(class = "alert alert-warning", msg)
  })

  # ---- C3. the existing design --------------------------------------------
  cmp_exist_csv <- reactive({
    tryCatch(read_design(input$cmp_exist_file, input$cmp_exist_text,
                         names(owea:::.ui_design_box(cmp_covariates())$design_box),
                         valcol = NULL),
             error = function(e) e)
  })

  cmp_existing <- reactive({
    if (!identical(input$cmp_start, "design")) return(NULL)
    d <- cmp_exist_csv()
    if (inherits(d, "error")) stop(conditionMessage(d), call. = FALSE)
    if (is.null(d)) stop("upload or paste your existing design.", call. = FALSE)
    n0 <- input$cmp_exist_n0
    e  <- owea:::.ui_existing_from_csv(d$support, d$val, d$valcol,
                                       n0 = if (is.finite(n0 %||% NA)) n0 else NULL)
    if (is.null(e$n0))
      stop("this design has weights, not counts, so it carries no sample size: ",
           "enter the number of runs already made (n0).", call. = FALSE)
    n1 <- input$cmp_exist_n1
    e$n1 <- if (is.finite(n1 %||% NA)) as.integer(n1) else e$n0
    e
  })

  observeEvent(cmp_exist_csv(), {
    d <- cmp_exist_csv()
    if (inherits(d, "error") || is.null(d)) return()
    if (identical(d$valcol, "count") && !is.finite(input$cmp_exist_n0 %||% NA))
      updateNumericInput(session, "cmp_exist_n0", value = sum(round(d$val)))
  })

  output$cmp_exist_msg <- renderUI({
    e <- tryCatch(cmp_existing(), error = function(e) e)
    if (inherits(e, "error"))
      return(div(class = "alert alert-danger", conditionMessage(e)))
    if (is.null(e)) return(NULL)
    div(class = "alert alert-success",
        sprintf("%d support point(s); n0 = %d already made, n1 = %d to add.",
                nrow(e$points), e$n0, e$n1))
  })

  # ---- C4. review ----------------------------------------------------------
  output$cmp_review_out <- renderUI({
    sp <- tryCatch(cmp_specs(), error = function(e) e)
    if (inherits(sp, "error"))
      return(div(class = "alert alert-danger", conditionMessage(sp)))
    cm <- cmp_comps(); al <- cmp_alpha(); al <- al / sum(al)
    gs <- cmp_grid_sizes()
    ex <- tryCatch(cmp_existing(), error = function(e) NULL)
    li <- list(
      sprintf("Covariates: %s.",
              paste(names(sp[[1]]$design_box), collapse = ", ")),
      sprintf("Candidate points at the coarsest stage: %s.",
              if (is.null(gs)) "?" else format(round(gs[1]), big.mark = ",")),
      sprintf("Objectives: %d.", length(cm)),
      sprintf("Weighting: %s.",
              if (identical(input$cmp_efficiency %||% "eff", "eff"))
                "applied to efficiencies (each objective divided by its own optimum)"
              else "applied to the raw criterion values"),
      if (!is.null(ex))
        sprintf("Augmenting an existing design: n0 = %d, n1 = %d.", ex$n0, ex$n1))
    comp_li <- lapply(seq_along(cm), function(j) {
      k <- length(tryCatch(owea:::.ui_coef_names(sp[[j]]), error = function(e) ""))
      tags$li(sprintf("%s — %s, %d parameter(s), %s-optimal%s, weight %.3f",
                      cm[[j]]$name, cm[[j]]$link, k,
                      if (cm[[j]]$p == 0) "D" else "A",
                      if (length(cm[[j]]$subset))
                        sprintf(" for parameter(s) %s",
                                paste(cm[[j]]$subset, collapse = ", ")) else "",
                      al[j]))
    })
    tagList(tags$ul(lapply(Filter(Negate(is.null), li), tags$li)),
            tags$b("The objectives:"), tags$ul(comp_li),
            if (identical(input$cmp_efficiency %||% "eff", "eff"))
              helpText(sprintf("This runs %d single-objective designs first (to ",
                               length(cm)),
                       "get each objective's own best value), then the compound ",
                       "design itself."))
  })

  # ---- compute -------------------------------------------------------------
  cmp_is_exact <- reactive(identical(input$cmp_design_type %||% "approx", "exact"))

  cmp_go <- reactiveVal(0)
  observeEvent(input$cmp_compute, cmp_go(cmp_go() + 1))
  observeEvent(input$cmp_compute, {
    rv$cmp_verify <- NULL; rv$cmp_sim <- list(); rv$cmp_show_sim <- FALSE
  })

  cmp_computed <- eventReactive(cmp_go(), {
    bad <- function(m) list(error = m)
    sp <- tryCatch(cmp_specs(), error = function(e) e)
    if (inherits(sp, "error")) return(bad(conditionMessage(sp)))
    cm <- cmp_comps()
    msg <- tryCatch(owea:::.uic_check(sp, cm, cmp_alpha()),
                    error = conditionMessage)
    if (!is.null(msg)) return(bad(msg))
    ex <- tryCatch(cmp_existing(), error = function(e) e)
    if (inherits(ex, "error")) return(bad(conditionMessage(ex)))
    exact <- cmp_is_exact()
    nn <- if (exact) as.numeric(input$cmp_n %||% NA) else NULL
    if (exact && !is.finite(nn))
      return(bad("enter the number of runs n for an exact design."))

    args <- tryCatch(
      owea:::.uic_solver_args(sp, cm, if (exact) "exact" else "design",
                              alpha = cmp_alpha(),
                              efficiency = identical(input$cmp_efficiency %||% "eff",
                                                     "eff"),
                              existing = ex, n = nn),
      error = function(e) e)
    if (inherits(args, "error")) return(bad(conditionMessage(args)))
    if (exact) args$seed <- as.integer(input$cmp_seed %||% 1)

    # The reference values Psi*_j cost one full multistage solve per component
    # and dominate a compound run, but they depend only on the objectives, the
    # region and any existing design -- not on alpha, n or the seed.  So if
    # none of that has changed since the last successful solve, hand the stored
    # psi_star back to the solver and skip the reference solves entirely.
    # Same numbers either way: this only reuses values already computed here.
    key <- owea:::.uic_psi_key(args)
    reused <- FALSE
    if (isTRUE(args$efficiency) && !is.null(cmp_psi_cache) &&
        identical(cmp_psi_cache$key, key)) {
      args$psi_star <- cmp_psi_cache$psi_star
      args$reference_bound <- cmp_psi_cache$reference_bound
      reused <- TRUE
    }

    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch(withProgress(message = "Computing the compound design…",
                            value = 0.5,
                            do.call(if (exact) compound_exact_design
                                    else compound_design, args)),
               error = function(e) structure(list(msg = conditionMessage(e)),
                                             class = "owea_bad")),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    if (inherits(res, "owea_bad")) return(bad(res$msg))
    ps <- suppressWarnings(as.numeric(res$psi_star))
    if (isTRUE(args$efficiency) && length(ps) && all(is.finite(ps)) && all(ps > 0))
      cmp_psi_cache <<- list(key = key, psi_star = ps,
                             reference_bound = suppressWarnings(as.numeric(res$reference_bound)))
    list(res = res, warns = warns, specs = sp, comps = cm, exact = exact,
         n = if (exact) as.integer(nn) else NA_integer_,
         cov_names = names(sp[[1]]$design_box), existing = ex, args = args,
         psi_reused = reused)
  }, ignoreInit = TRUE)

  observeEvent(cmp_computed(), {
    c2 <- cmp_computed()
    if (is.null(c2$error)) goto("cmp_results")
    else showNotification(c2$error, type = "error", duration = 10)
  })

  # ---- C5. results ---------------------------------------------------------
  output$cmp_status <- renderUI({
    c2 <- cmp_computed()
    if (!is.null(c2$error))
      return(div(class = "alert alert-danger", tags$b("Could not compute. "),
                 c2$error))
    r <- c2$res
    # why a re-run can be much faster than the first one
    reuse <- if (isTRUE(c2$psi_reused))
      tags$div(tags$small(paste(
        "The objectives, the region and the existing design are unchanged since",
        "the last run, so the reference values (Psi*) were reused and the",
        "per-objective reference solves were skipped. Same numbers, less work.")))
    if (c2$exact)
      return(div(class = "alert alert-success",
        tags$b("Exact compound design found. "),
        sprintf(paste0("%s = %.6f (approximate compound design: %.6f); efficiency ",
                       "at least %.4f%% of the compound optimum. %d exchange(s) accepted."),
                if (isTRUE(r$efficiency_weighted))
                  "Weighted average efficiency" else "Weighted criterion",
                r$criterion, r$criterion_approx, 100 * r$efficiency_exact_lower_bound,
                r$exchanges),
        if (length(c2$warns))
          tags$ul(lapply(unique(c2$warns),
                         function(w) tags$li(.friendly_warn(w)))),
        reuse))
    ok <- isTRUE(r$converged)
    # Psi_alpha is a weighted average of the component efficiencies, so it is
    # bounded by 1, but 1 is reached only if one design were optimal for every
    # objective at once.  Say what the number is, not just its scale.
    val_ui <- if (isTRUE(r$efficiency_weighted))
      sprintf(paste0("Weighted average of the objectives' efficiencies = %.6f. ",
                     if (ok) paste0("This is the highest value any single design can ",
                                    "reach for these objectives with these weights ",
                                    "(certified: max sensitivity %.2e, 0 at the optimum). ")
                     else "(max sensitivity %.2e; 0 at the optimum). ",
                     "The average cannot exceed 1 and would equal 1 only if one design ",
                     "were optimal for every objective simultaneously; the shortfall ",
                     "from 1 is the price of serving all objectives with one design."),
              r$criterion, r$max_d)
    else
      sprintf("Weighted criterion = %.6f  (max sensitivity %.2e; 0 at the optimum).",
              r$criterion, r$max_d)
    div(class = if (ok) "alert alert-success" else "alert alert-warning",
        tags$b(if (ok) "Compound design found. " else "Did not converge. "),
        val_ui,
        if (length(c2$warns))
          tags$ul(lapply(unique(c2$warns), function(w) tags$li(.friendly_warn(w)))),
        reuse)
  })

  cmp_design_df <- reactive({
    c2 <- cmp_computed(); req(is.null(c2$error))
    d <- as.data.frame(c2$res$support)
    names(d) <- c2$cov_names
    if (c2$exact) d$count  <- as.integer(c2$res$counts)
    else          d$weight <- as.numeric(c2$res$weights)
    d
  })

  output$cmp_design_tbl <- DT::renderDT({
    c2 <- cmp_computed(); req(is.null(c2$error))
    d <- cmp_design_df()
    num <- vapply(d, is.numeric, logical(1)); num["count"] <- FALSE
    d[num] <- lapply(d[num], function(z) round(z, 4))
    cap <- if (c2$exact)
             sprintf("Exact compound design: %d run(s) over %d support point(s).",
                     c2$n, nrow(d))
           else
             sprintf(paste0("Approximate compound design: %d support point(s); ",
                            "weights sum to 1."), nrow(d))
    DT::datatable(d, rownames = FALSE, caption = cap,
                  options = list(dom = "t", paging = FALSE))
  })

  output$cmp_dl_csv <- downloadHandler(
    filename = function() "owea_compound_design.csv",
    content  = function(file) {
      old <- options(digits = 15); on.exit(options(old))
      utils::write.csv(cmp_design_df(), file, row.names = FALSE)
    })

  output$cmp_eff_help <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    div(class = "alert alert-info",
        tags$b("What each objective got."),
        " 'criterion' is the objective's own criterion value at this design, ",
        "computed exactly as the single-criterion Design tab reports it ",
        "(D: log det of the covariance / number of parameters; A: trace / number ",
        "of parameters; smaller is better)",
        if (isTRUE(c2$res$efficiency_weighted))
          tags$span("; 'optimal' is the best value that objective could achieve ",
                    "on its own, and 'efficiency' is the fraction of that best this ",
                    "one design delivers, reported as a guaranteed lower bound ",
                    "relative to the true optimum (it accounts for the reference ",
                    "solve possibly stopping short of the optimum).")
        else tags$span("."))
  })
  output$cmp_summary_tbl <- renderTable({
    c2 <- cmp_computed(); req(is.null(c2$error))
    owea:::.uic_summary_table(c2$res)
  }, rownames = TRUE, digits = 6)

  output$cmp_cross_help <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    if (is.null(c2$res$cross_efficiency)) return(NULL)
    div(class = "alert alert-info",
        tags$b("What a single-objective design would have cost."),
        " Each row is a design, each column an objective; the entries are ",
        "efficiencies. Row j is the design optimised for objective j alone — ",
        "1.000 on its own diagonal, and whatever it happens to give under the ",
        "others. The last row is the compound design you just computed.")
  })
  output$cmp_cross_tbl <- renderTable({
    c2 <- cmp_computed(); req(is.null(c2$error))
    ct <- owea:::.uic_cross_table(c2$res)
    if (is.null(ct)) return(NULL)
    ct
  }, rownames = TRUE, digits = 3)

  output$cmp_plot <- renderPlot({
    c2 <- cmp_computed(); req(is.null(c2$error))
    plot_design(c2$res, cov_names = c2$cov_names, main = "Compound design")
  })
  output$cmp_dl_png <- downloadHandler(
    filename = function() "owea_compound_design.png",
    content  = function(file) {
      grDevices::png(file, width = 900, height = 650, res = 110)
      on.exit(grDevices::dev.off())
      c2 <- cmp_computed()
      plot_design(c2$res, cov_names = c2$cov_names, main = "Compound design")
    })

  output$cmp_info_pick <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    nm <- names(c2$res$information)
    tagList(
      helpText("Each objective has its own model, so its own information ",
               "matrix — and they need not be the same size."),
      selectInput("cmp_info_which", "Objective", choices = nm))
  })
  output$cmp_info_tbl <- renderTable({
    c2 <- cmp_computed(); req(is.null(c2$error))
    w <- input$cmp_info_which %||% names(c2$res$information)[1]
    M <- c2$res$information[[w]]
    if (is.null(M)) return(NULL)
    as.data.frame(round(M, 5))
  }, rownames = TRUE, digits = 5)

  # ---- post-result: verify (approximate) or simulate (exact) --------------
  output$cmp_post_result <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    if (c2$exact)
      tagList(hr(),
              actionButton("cmp_sim_open", "Simulation study", class = "btn-info"),
              actionButton("reset_btn", "Reset"),
              uiOutput("cmp_sim_panel"))
    else
      uiOutput("cmp_verify_panel")
  })

  # ---- verify a design under the compound criterion ------------------------
  output$cmp_verify_panel <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    tagList(
      hr(),
      h4("Score another design"),
      helpText("Paste or upload any design — an edited version of this one, a ",
               "factorial, a design from the literature — and see what it scores ",
               "under the same compound criterion."),
      fileInput("cmp_verify_file", "Upload a CSV", accept = ".csv"),
      textAreaInput("cmp_verify_manual", "…or edit the design here", rows = 6,
                    value = design_csv_text(cmp_design_df())),
      actionButton("cmp_verify_btn", "Score this design", class = "btn-info"),
      actionButton("reset_btn", "Reset"),
      uiOutput("cmp_verify_out"))
  })

  # scoring a design scans the WHOLE region at the finest step, exactly as the
  # classical branch's verify does, so it is gated the same way
  run_cmp_verify <- function(criterion_only = FALSE) {
    c2 <- cmp_computed(); if (is.null(c2) || !is.null(c2$error)) return()
    out <- tryCatch({
      d <- read_design(input$cmp_verify_file, input$cmp_verify_manual,
                       c2$cov_names, valcol = "weight")
      if (is.null(d)) stop("paste or upload a design first.", call. = FALSE)
      w <- as.numeric(d$val)
      if (any(w < 0)) stop("weights must be nonnegative.", call. = FALSE)
      w <- w / sum(w)
      args <- owea:::.uic_solver_args(
        c2$specs, c2$comps, "verify", alpha = cmp_alpha(),
        efficiency = isTRUE(c2$res$efficiency_weighted),
        existing = c2$existing,
        psi_star = if (isTRUE(c2$res$efficiency_weighted)) c2$res$psi_star else NULL,
        reference_bound = if (isTRUE(c2$res$efficiency_weighted)) c2$res$reference_bound else NULL)
      withProgress(message = "Scoring the design…", value = 0.5,
        do.call(compound_criterion,
                c(list(support = d$support, weights = w,
                       criterion_only = isTRUE(criterion_only)), args)))
    }, error = function(e) e)
    rv$cmp_verify <- out
  }

  # gate on the size of the finest-step grid, offering the same three choices
  # compound_criterion() itself offers at the console
  observeEvent(input$cmp_verify_btn, {
    c2 <- cmp_computed(); req(is.null(c2$error))
    N <- tryCatch(owea:::.uic_verify_points(c2$specs), error = function(e) 0)
    if (is.finite(N) && N > 1e6)
      showModal(modalDialog(title = "Large design space",
        sprintf(paste0("Scoring a design evaluates the compound sensitivity at ",
          "all %s design points of the finest-step grid; building it may be ",
          "slow. The criterion itself needs no grid."),
          format(round(N), big.mark = ",")),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("cmp_verify_full", "Proceed anyway (full check)"),
          actionButton("cmp_verify_crit_only",
                       "Criterion only (skip the optimality check)",
                       class = "btn-warning"))))
    else run_cmp_verify()
  })
  observeEvent(input$cmp_verify_full,      { removeModal(); run_cmp_verify(FALSE) })
  observeEvent(input$cmp_verify_crit_only, { removeModal(); run_cmp_verify(TRUE) })

  output$cmp_verify_out <- renderUI({
    v <- rv$cmp_verify
    if (is.null(v)) return(NULL)
    if (inherits(v, "error"))
      return(div(class = "alert alert-danger", conditionMessage(v)))
    c2 <- cmp_computed()
    eff <- v$efficiency_lower_bound
    # the headline value, shared by both kinds of run
    head_ui <- tagList(
      tags$b(sprintf("%s = %.6f. ",
                     if (isTRUE(c2$res$efficiency_weighted))
                       "Weighted average efficiency" else "Weighted criterion",
                     v$criterion)),
      if (all(is.finite(eff)))
        tags$span(sprintf("Per objective: %s.",
                          paste(sprintf("%s %.3f", names(eff), eff),
                                collapse = ", "))))
    share_ui <- if (is.finite(c2$res$criterion))
      tags$span(sprintf(" It reaches %.1f%% of the computed design's value.",
                        100 * v$criterion / c2$res$criterion))
    # criterion-only run: no grid was built, so optimality was not assessed
    if (is.null(v$max_d))
      return(div(class = "alert alert-info", head_ui, tags$br(),
                 tags$b("Optimality was NOT assessed"),
                 " (the sensitivity scan was skipped, so no grid was built).",
                 share_ui))
    div(class = if (isTRUE(v$is_optimal)) "alert alert-success"
                else "alert alert-warning",
        head_ui,
        tags$br(),
        sprintf("Maximum sensitivity %.2e — %s", v$max_d,
                if (isTRUE(v$is_optimal))
                  "this design is compound-optimal over the grid."
                else "this design is NOT compound-optimal; a better one exists."),
        if (!isTRUE(v$is_optimal)) share_ui)
  })

  # ---- simulation study (exact compound designs) ---------------------------
  # A compound design serves several models, so the study is run once PER
  # component: generate under that component's model at the design, fit the
  # same model, and report the per-parameter MSE.  The same three designs the
  # classical branch compares are offered -- the computed design, a random
  # (SRS) design, and one the user supplies.
  observeEvent(input$cmp_sim_open, { rv$cmp_show_sim <- TRUE })

  output$cmp_sim_panel <- renderUI({
    if (!isTRUE(rv$cmp_show_sim)) return(NULL)
    c2 <- cmp_computed(); req(is.null(c2$error))
    sp <- c2$specs; cm <- c2$comps
    any_identity <- any(vapply(sp, function(s) identical(s$link, "identity"),
                               logical(1)))
    tagList(
      hr(),
      h4("Simulation study"),
      helpText("Each objective is simulated under its OWN model: data are ",
               "generated from that model at the design, the same model is ",
               "fitted back, and the mean squared error of each parameter is ",
               "averaged over the replications. Objectives are judged ",
               "separately, so the numbers are comparable across designs."),
      div(class = "alert alert-info",
          tags$b("True parameter values."),
          " These need not equal the assumed values used to build the design — ",
          "setting them apart is how you check how much the design suffers when ",
          "the guess was wrong."),
      lapply(seq_along(sp), function(j) {
        cn <- tryCatch(owea:::.ui_coef_names(sp[[j]]), error = function(e) NULL)
        if (is.null(cn)) return(NULL)
        th <- cm[[j]]$theta
        wellPanel(
          tags$b(sprintf("%s (%s)", cm[[j]]$name, sp[[j]]$link)),
          fluidRow(lapply(seq_along(cn), function(i)
            column(3, numericInput(
              sprintf("cmp_sim_theta_%d_%d", j, i), cn[i], step = 0.1,
              value = isolate(input[[sprintf("cmp_sim_theta_%d_%d", j, i)]])
                      %||% (if (!is.null(th) && length(th) >= i) th[i] else 0))))))
      }),
      if (any_identity)
        numericInput("cmp_sim_sigma", "Residual SD (linear objectives)",
                     value = 1, min = 1e-8, step = 0.1),
      fluidRow(
        column(6, numericInput("cmp_sim_nsim", "Replications", value = 500,
                               min = 10, step = 50)),
        column(6, numericInput("cmp_sim_seed", "Random seed", value = 1,
                               step = 1))),
      actionButton("cmp_run_sim", "Run for this design", class = "btn-primary"),
      actionButton("cmp_run_srs", "Compare with a random design"),
      br(), br(),
      helpText("…or supply your own design with the same number of runs:"),
      fileInput("cmp_sim_custom_file", "Upload a CSV", accept = ".csv"),
      textAreaInput("cmp_sim_custom", NULL, rows = 5,
                    placeholder = "x1,x2,count\n-2,-2,10\n2,2,20"),
      actionButton("cmp_run_custom", "Compare with my design"),
      br(), br(),
      uiOutput("cmp_sim_pick"),
      tableOutput("cmp_sim_tbl"),
      uiOutput("cmp_sim_note"),
      uiOutput("cmp_sim_msg"))
  })

  # the true values, one vector per component
  cmp_sim_theta <- reactive({
    c2 <- cmp_computed(); req(is.null(c2$error))
    lapply(seq_along(c2$specs), function(j) {
      cn <- tryCatch(owea:::.ui_coef_names(c2$specs[[j]]), error = function(e) NULL)
      if (is.null(cn)) return(NULL)
      vapply(seq_along(cn), function(i)
        as.numeric(input[[sprintf("cmp_sim_theta_%d_%d", j, i)]] %||% 0),
        numeric(1))
    })
  })

  # one simulation run over all components, for a given (support, counts)
  cmp_sim_run <- function(support, counts) {
    c2 <- cmp_computed()
    withProgress(message = "Simulating…", value = 0.5,
      owea:::.uic_simulate(c2$specs, c2$comps, support = support,
                           counts = counts, true_theta = cmp_sim_theta(),
                           sigma = as.numeric(input$cmp_sim_sigma %||% 1),
                           existing = c2$existing,
                           nsim = as.integer(input$cmp_sim_nsim %||% 500),
                           seed = as.integer(input$cmp_sim_seed %||% 1)))
  }

  observeEvent(input$cmp_run_sim, {
    c2 <- cmp_computed(); req(is.null(c2$error))
    rv$cmp_sim$design <- tryCatch(
      cmp_sim_run(c2$res$support, as.integer(c2$res$counts)),
      error = function(e) e)
  })

  observeEvent(input$cmp_run_srs, {
    c2 <- cmp_computed(); req(is.null(c2$error))
    rv$cmp_sim$srs <- tryCatch({
      pool <- candidate_grid(c2$specs[[1]]$design_box, c2$specs[[1]]$finest)
      set.seed(as.integer(input$cmp_sim_seed %||% 1))
      pick <- pool[sample.int(nrow(pool), c2$n, replace = TRUE), , drop = FALSE]
      key <- apply(pick, 1, paste, collapse = "\r")
      u   <- !duplicated(key)
      sup <- pick[u, , drop = FALSE]
      cnt <- as.integer(table(factor(key, levels = key[u])))
      cmp_sim_run(sup, cnt)
    }, error = function(e) e)
  })

  observeEvent(input$cmp_run_custom, {
    c2 <- cmp_computed(); req(is.null(c2$error))
    rv$cmp_sim$custom <- tryCatch({
      d <- read_design(input$cmp_sim_custom_file, input$cmp_sim_custom,
                       c2$cov_names, valcol = "count")
      if (is.null(d)) stop("paste or upload a design first.", call. = FALSE)
      cnt <- as.numeric(d$val)
      if (any(cnt < 0) || any(abs(cnt - round(cnt)) > 1e-8))
        stop("the counts must be nonnegative whole numbers.", call. = FALSE)
      if (sum(cnt) != c2$n)
        stop(sprintf("the counts must sum to %d, the same number of runs as the ",
                     c2$n), "computed design.", call. = FALSE)
      cmp_sim_run(d$support, as.integer(round(cnt)))
    }, error = function(e) e)
  })

  output$cmp_sim_pick <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    if (!length(rv$cmp_sim)) return(NULL)
    selectInput("cmp_sim_which", "Show the MSEs for",
                choices = vapply(c2$comps, function(z) z$name, character(1)))
  })

  # which objective the table and the chart are showing
  cmp_sim_which <- reactive({
    c2 <- cmp_computed(); req(is.null(c2$error))
    nms <- vapply(c2$comps, function(z) z$name, character(1))
    j <- match(input$cmp_sim_which %||% nms[1], nms)
    list(j = if (is.na(j)) 1L else j, names = nms)
  })

  # the comparison for that ONE objective, assembled once for both renderings
  # (mean and median squared error), exactly as in the classical branch
  cmp_sim_mats <- reactive({
    s <- rv$cmp_sim
    if (!length(s)) return(list(`Mean squared error` = NULL,
                                `Median squared error` = NULL))
    j <- cmp_sim_which()$j
    slots <- list(`This design` = s$design, Random = s$srs, Custom = s$custom)
    grab <- function(what) owea:::.ui_sim_matrix(
      lapply(slots, function(v)
        if (!is.list(v) || inherits(v, "error")) NULL
        else owea:::.uic_sim_mse(v[[j]], what)))
    list(`Mean squared error` = grab("mse"),
         `Median squared error` = grab("medse"))
  })

  # the compound counterpart of SIM_REF: the design the app computed
  CMP_SIM_REF <- "This design"

  output$cmp_sim_tbl <- renderTable({
    m <- cmp_sim_mats(); M <- m[["Mean squared error"]]
    if (is.null(M)) return(NULL)
    df <- data.frame(Parameter = colnames(M), stringsAsFactors = FALSE)
    # one column per design, then its ratio to this design, exactly as in the
    # classical branch -- above 1, this design has the smaller error
    add <- function(df, A, unit) {
      if (is.null(A)) return(df)
      for (nm in rownames(A)) df[[sprintf("%s (%s)", nm, unit)]] <- A[nm, ]
      R <- owea:::.ui_sim_ratio(A, CMP_SIM_REF)
      if (!is.null(R))
        for (nm in rownames(R))
          df[[sprintf("%s / %s (%s)", nm, CMP_SIM_REF, unit)]] <- R[nm, ]
      df
    }
    df <- add(df, M, "MSE")
    add(df, m[["Median squared error"]], "median SE")
  }, digits = 6)

  # How many replicates converged, and a warning when a few wild fits are
  # driving the MSE.  A D-optimal design concentrates on few distinct points,
  # so for a binary or count response some replicates can separate; the fitter
  # still "converges" but returns enormous coefficients, and a handful of those
  # dominate a mean squared error.  Without this note the table can look as
  # though the optimal design were the worse one, when the asymptotic criterion
  # (and the median estimate) say otherwise.
  output$cmp_sim_note <- renderUI({
    c2 <- cmp_computed(); req(is.null(c2$error))
    s <- rv$cmp_sim
    if (!length(s)) return(NULL)
    j <- cmp_sim_which()$j
    ratio_ui <- if (!is.null(owea:::.ui_sim_ratio(cmp_sim_mats()[[1]],
                                                  CMP_SIM_REF)))
      helpText(tags$b(sprintf("The \"/ %s\" columns are ratios.", CMP_SIM_REF)),
               paste("Above 1 the computed design has the smaller error for",
                     "that parameter (2.0 = twice the error); below 1 the other",
                     "design does. The errors themselves are on each",
                     "parameter's own scale, so only the ratio compares",
                     "designs across a row."))
    lab <- c(design = "this design", srs = "random", custom = "custom")
    conv <- character(0); wild <- character(0)
    for (slot in names(s)) {
      v <- s[[slot]]
      if (!is.list(v) || inherits(v, "error")) next
      sj <- v[[j]]
      if (is.null(sj) || inherits(sj, "error")) next
      conv <- c(conv, sprintf("%s %d/%d", lab[[slot]], sj$n_converged, sj$nsim))
      dom <- owea:::.ui_sim_dominance(sj)     # shared with the classical branch
      if (is.finite(dom) && dom > 0.2)
        wild <- c(wild, sprintf("%s (one replicate contributes %.0f%% of an MSE)",
                                lab[[slot]], 100 * dom))
    }
    tagList(
      ratio_ui,
      if (length(conv))
        helpText(sprintf("Replications that converged: %s.",
                         paste(conv, collapse = ", "))),
      if (length(wild))
        div(class = "alert alert-warning",
            tags$b("Read these MSEs with care."),
            sprintf(" Some fits returned extreme estimates — %s. ",
                    paste(wild, collapse = ", ")),
            "That is separation: a design concentrated on few distinct points ",
            "can, in a given replicate, leave a binary or count response ",
            "perfectly predicted, and the fitted coefficients then run away. ",
            "A handful of such replicates dominates a mean squared error, so ",
            "the design that is better asymptotically can look worse here. ",
            tags$small("The median-squared-error columns are not affected by a ",
                       "handful of wild replicates; compare those, and the ",
                       "Efficiency tab, which reports the large-sample ",
                       "criterion.")))
  })

  output$cmp_sim_msg <- renderUI({
    s <- rv$cmp_sim
    if (!length(s)) return(NULL)
    c2 <- cmp_computed()
    nms <- vapply(c2$comps, function(z) z$name, character(1))
    msgs <- character(0)
    for (slot in names(s)) {
      v <- s[[slot]]
      if (inherits(v, "error")) {
        msgs <- c(msgs, sprintf("%s: %s", slot, conditionMessage(v)))
      } else if (is.list(v)) {
        for (j in seq_along(v))
          if (inherits(v[[j]], "error"))
            msgs <- c(msgs, sprintf("%s, %s: %s", slot, nms[j],
                                    conditionMessage(v[[j]])))
      }
    }
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
    else if (grepl("not among the design points", m, fixed = TRUE))
      "The design's support points are not grid points (for example, a design from the continuous search). The sensitivity was evaluated at every grid point with the design taken exactly as given, so the check is valid."
    else if (grepl("continuous search did NOT converge", m, fixed = TRUE))
      "The continuous search stopped before certifying optimality; the design is close to optimal but not certified. In the R package, raise n_starts or max_iter, or supply init_points."
    else if (grepl("did NOT converge", m, fixed = TRUE))
      "The solver did not fully converge; try a finer grid step or fewer parameters of interest."
    else m
  }, character(1), USE.NAMES = FALSE)
}

shinyApp(ui, server)
