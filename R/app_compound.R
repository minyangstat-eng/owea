# ===========================================================================
# app_compound.R -- the translator layer for the web app's COMPOUND branch.
#
# ADDITIVE FILE.  Nothing here touches the classical `.ui_*` helpers in
# R/app.R; it reuses them.  All helpers are internal (`.uic_*`, reached from
# the app via owea:::) and pure, so they are unit-tested without launching
# Shiny -- exactly the arrangement R/app.R uses for the classical branch.
#
# The compound problem the app poses:
#
#   * ONE experiment, so ONE set of covariates and one design region, shared
#     by every component;
#   * J components, each with its OWN model (family, terms, assumed theta),
#     its OWN criterion (D or A) and its OWN parameters of interest;
#   * weights alpha_j saying how much each component matters.
#
# The key reuse: a component's model is nothing more than .ui_model_spec()
# called with the SHARED covariates and that component's own link and terms.
# So the whole classical translation layer -- factor coding, within-kind
# indices, interaction pairs, coefficient names -- carries over per component
# and cannot drift between the two branches.
# ===========================================================================

# Cap on the number of components the app offers (see MAX_COV for covariates).
.UIC_MAX_COMP <- 4L

# ---- wizard step graph ----------------------------------------------------
# The compound branch's ordered step ids.  Unlike the classical graph this one
# is fixed: an existing design is an option INSIDE the options step rather than
# a step of its own, so nothing reshapes the chain underneath the user.
.uic_wizard_steps <- function() {
  c("cmp_cov", "cmp_models", "cmp_options", "cmp_review", "cmp_results")
}

# ---- one component --------------------------------------------------------
# `covariates` is the SHARED covariate list (as .ui_design_box() takes it).
# `comp` is one component's UI state:
#   link, ncat, coding      -- the model family
#   interactions, quadratics-- terms, in the same encoding .ui_model_spec() takes
#   p                       -- 0 = D, 1 = A
#   subset                  -- integer parameter indices, or NULL/empty = all
#   theta                   -- assumed values (ignored for the identity link)
#   alpha                   -- the component's weight
#   name                    -- label used in the printout and the results tables
#
# Returns the .ui_model_spec() spec for that component.
.uic_component_spec <- function(covariates, comp) {
  .ui_model_spec(covariates,
                 interactions = if (is.null(comp$interactions)) list()
                                else comp$interactions,
                 quadratics   = if (is.null(comp$quadratics)) integer(0)
                                else comp$quadratics,
                 link   = if (is.null(comp$link)) "identity" else comp$link,
                 ncat   = comp$ncat,
                 coding = if (is.null(comp$coding)) "zero-sum" else comp$coding)
}

# All component specs, in order.  Errors if there are none, or if the
# components somehow disagree about the design region (they cannot, since they
# are built from one covariate list, but the check documents the invariant).
.uic_specs <- function(covariates, comps) {
  if (!length(comps))
    stop("add at least one component.", call. = FALSE)
  specs <- lapply(comps, function(cm) .uic_component_spec(covariates, cm))
  db1 <- specs[[1]]$design_box
  for (s in specs)
    if (!identical(names(s$design_box), names(db1)))
      stop("every component must use the same covariates.", call. = FALSE)
  specs
}

# ---- validation -----------------------------------------------------------
# Pre-flight checks so the wizard shows an inline message instead of a hard
# solver error.  Returns NULL when everything is fine, else one message.
.uic_check <- function(specs, comps, alpha) {
  J <- length(specs)
  if (J < 1L) return("add at least one component.")
  a <- suppressWarnings(as.numeric(alpha))
  if (length(a) != J)
    return(sprintf("give one weight per component (%d).", J))
  if (any(!is.finite(a)) || any(a < 0))
    return("every weight must be a finite, nonnegative number.")
  if (sum(a) <= 0) return("the weights cannot all be zero.")
  for (j in seq_len(J)) {
    sp <- specs[[j]]; cm <- comps[[j]]
    nm <- if (is.null(cm$name)) sprintf("component %d", j) else cm$name
    k  <- length(.ui_coef_names(sp))
    msg <- .ui_check_theta(cm$theta, sp)
    if (!is.null(msg)) return(sprintf("%s: %s", nm, msg))
    if (length(cm$subset)) {
      s <- suppressWarnings(as.integer(cm$subset))
      if (any(is.na(s)) || any(s < 1L) || any(s > k))
        return(sprintf("%s: the parameters of interest must be among 1..%d.",
                       nm, k))
    }
    if (!is.null(cm$p) && !(as.integer(cm$p) %in% c(0L, 1L)))
      return(sprintf("%s: the criterion must be D (0) or A (1).", nm))
  }
  NULL
}

# ---- one component's argument list ----------------------------------------
# The entry compound_design()/compound_criterion() expects in `components`.
# NOTE it deliberately carries NO design_box: the design region is a top-level
# argument of compound_design(), shared by every component.
.uic_component_args <- function(spec, comp) {
  k <- length(.ui_coef_names(spec))
  a <- list(link = spec$link, ncat = spec$ncat, coding = .ui_coding(spec),
            f = spec$f, x = spec$x, fx = spec$fx, xx = spec$xx, ff = spec$ff,
            p = .check_criterion(if (is.null(comp$p)) 0L else comp$p))

  # theta: NULL for the identity link, required and length-checked otherwise
  if (identical(spec$link, "identity")) {
    a$theta <- NULL
  } else {
    th <- comp$theta
    if (is.null(th) || length(th) != k || any(!is.finite(th)))
      stop(sprintf("the %s model needs %d finite assumed parameter value(s).",
                   spec$link, k), call. = FALSE)
    a$theta <- as.numeric(th)
  }

  # parameters of interest; an EMPTY subset stays absent (= all parameters)
  if (length(comp$subset)) {
    s <- sort(unique(as.integer(comp$subset)))
    if (any(is.na(s)) || any(s < 1L) || any(s > k))
      stop(sprintf("the parameters of interest must be among 1..%d.", k),
           call. = FALSE)
    a$subset <- s
  }
  if (!is.null(comp$name)) a$name <- as.character(comp$name)
  a
}

# ---- the single argument assembler ----------------------------------------
# Every compound solver call the app makes goes through here, so the models --
# above all each component's factor `coding` -- can never drift between the
# design and the verification.  Mirrors .ui_solver_args() for the classical
# branch, including its rule that with no existing design the xi0_* / n0 / n1
# arguments are OMITTED entirely.
#
#   specs      : .uic_specs() output
#   comps      : the parallel component UI list
#   target     : "design" -> compound_design(), "verify" -> compound_criterion()
#   alpha      : component weights (normalised by the solver)
#   efficiency : TRUE = weight efficiencies Psi_j/Psi_j*, FALSE = raw Psi_j
#   existing   : NULL, or .ui_existing_from_csv() output
#   psi_star   : reference values to REUSE for a verify call, so a pasted design
#                is scored on exactly the scale the computed design was
#
# Returns a named list to do.call() into.
#   n          : exact designs only, the runs to allocate
.uic_solver_args <- function(specs, comps,
                             target = c("design", "exact", "verify"),
                             alpha = NULL, efficiency = TRUE,
                             existing = NULL, psi_star = NULL,
                             reference_bound = NULL, n = NULL) {
  target <- match.arg(target)
  J <- length(specs)
  if (J != length(comps))
    stop("every component needs its own specification.", call. = FALSE)

  args <- list(
    components = lapply(seq_len(J),
                        function(j) .uic_component_args(specs[[j]], comps[[j]])),
    alpha      = if (is.null(alpha)) rep(1 / J, J) else as.numeric(alpha),
    efficiency = isTRUE(efficiency))

  # the design region -- shared, so taken from the first spec
  sp1 <- specs[[1]]
  args$design_box <- sp1$design_box
  if (target %in% c("design", "exact")) {
    args$step_sequence <- sp1$step_sequence
  } else {
    args$step <- sp1$finest            # a single step, as verify wants
    # the app gates the grid size itself (with a modal, via .uic_verify_points),
    # so the package cap must not fire a SECOND time -- inside a running Shiny
    # session interactive() is TRUE and compound_criterion()'s menu() would
    # block on the console with nobody to answer it
    args$max_points <- Inf
    if (!is.null(psi_star)) args$psi_star <- as.numeric(psi_star)
    # the reference designs' certified bounds travel with psi_star, so the
    # verify scores are the same lower bounds the design step reported
    if (!is.null(reference_bound)) args$reference_bound <- as.numeric(reference_bound)
  }

  # exact designs allocate n runs
  if (identical(target, "exact")) {
    if (is.null(n) || !is.finite(n) || n < 1)
      stop("an exact design needs the number of runs n (>= 1).", call. = FALSE)
    args$n <- as.integer(round(n))
  }

  # the existing design
  if (!is.null(existing)) {
    e <- .ui_check_existing(existing, sp1)
    args$xi0_points  <- e$points
    args$xi0_weights <- e$weights
    args$n0          <- e$n0
    args$n1          <- e$n1
  }
  args
}

# ---- R code that reproduces a compound computation -------------------------
# The compound counterpart of .ui_r_code() (app.R): the objectives as a
# `components` variable, the compound_design() / compound_exact_design() call
# the app made (its cached psi_star / reference_bound dropped, so the script
# recomputes the reference designs), the compound_criterion() scoring of the
# result for an approximate design, and, for an exact design whose simulation
# study has been run, that study -- each objective under its own model with
# simulate_design(), for the design and every design it was compared with.
#   args        : .uic_solver_args() output (+ seed for an exact design)
#   verify_args : .uic_solver_args(..., "verify") output, or NULL
#   sim         : NULL, or list(theta = list of true values per objective,
#                 sigma, nsim, seed, existing = list(points, counts) or NULL,
#                 designs = named list of list(support, counts))
.uic_r_code <- function(args, exact = FALSE, verify_args = NULL, sim = NULL) {
  ver <- tryCatch(as.character(utils::packageVersion("owea")), error = function(e) "")
  fn  <- if (exact) "compound_exact_design" else "compound_design"
  comps <- args$components
  comp_lines <- unlist(lapply(seq_along(comps), function(j) {
    cj <- comps[[j]]; nm <- cj$name; cj$name <- NULL
    cl <- .ui_fmt_call("list", cj,
                       first = if (!is.null(nm)) sprintf("name = \"%s\"", nm)
                               else character(0))
    cl <- paste0("  ", cl)
    cl[length(cl)] <- paste0(cl[length(cl)], if (j < length(comps)) "," else "")
    cl
  }))
  lines <- c(
    sprintf("## R code generated by the owea app (owea %s): it reproduces the computed compound design.", ver),
    sprintf("## Edit any argument and rerun; see ?%s for all options.", fn),
    "library(owea)", "",
    "## the objectives: one model, criterion and set of parameters of interest each",
    "components <- list(", comp_lines, ")", "")
  dargs <- args[setdiff(names(args), c("components", "psi_star", "reference_bound"))]
  dc <- .ui_fmt_call(fn, dargs, first = "components = components")
  lines <- c(lines, paste0("res <- ", dc[1]), dc[-1],
             "print(res)   # the design, and each objective's criterion and efficiency")
  if (!exact && !is.null(verify_args)) {
    vargs <- verify_args[setdiff(names(verify_args),
                                 c("components", "psi_star", "reference_bound"))]
    raw <- if (isTRUE(args$efficiency))
             c("psi_star = res$psi_star", "reference_bound = res$reference_bound")
           else character(0)
    vc <- .ui_fmt_call("compound_criterion", vargs,
                       first = c("support = res$support", "weights = res$weights",
                                 "components = components"), raw = raw)
    lines <- c(lines, "",
               "## Score the design (or any other) under the same compound criterion; the",
               "## reference values of the objectives are reused from the design above",
               paste0("v <- ", vc[1]), vc[-1],
               "v$criterion                # the weighted average of the objectives' efficiencies",
               "v$efficiency_lower_bound   # one certified value per objective",
               "v$max_d                    # max sensitivity over the grid (0 at the optimum)")
  }
  if (exact && !is.null(sim)) lines <- c(lines, "", .uic_r_code_sim(args, sim))
  paste(lines, collapse = "\n")
}

# The simulation block of .uic_r_code(): for every objective, simulate_design()
# at the exact design and at every compared design, then the mean squared
# errors side by side.
.uic_r_code_sim <- function(args, sim) {
  comps <- args$components; J <- length(comps)
  out <- c("## Simulation study: each objective under its OWN model -- responses are",
           "## generated from that model at the design, the same model is fitted back,",
           "## and the mean squared error of each parameter is averaged over replications")
  common <- list(design_box = args$design_box, nsim = as.integer(sim$nsim),
                 seed = as.integer(sim$seed))
  if (!is.null(sim$existing)) {
    common$existing <- list(points = unname(as.matrix(sim$existing$points)),
                            counts = as.integer(sim$existing$counts))
    out <- c(out, "## the existing first stage is simulated and analysed together with the new runs")
  }
  for (j in seq_len(J))
    out <- c(out, sprintf("theta_true_%d <- %s   # true values, objective %d", j,
                          .ui_fmt_arg(as.numeric(
                            if (is.null(sim$theta[[j]])) comps[[j]]$theta
                            else sim$theta[[j]])), j))
  vars <- character(0)
  for (nm in names(sim$designs)) {
    dsg <- sim$designs[[nm]]
    var <- tolower(gsub("[^A-Za-z0-9]+", "_", nm)); vars[nm] <- var
    out <- c(out, "",
             sprintf("## the %s design the app compared with (%d runs)", nm, sum(dsg$counts)),
             paste0(var, "_support <- ", .ui_fmt_arg(unname(as.matrix(dsg$support)))),
             paste0(var, "_counts  <- ", .ui_fmt_arg(as.integer(dsg$counts))))
  }
  for (j in seq_len(J)) {
    cj  <- comps[[j]]
    mod <- cj[intersect(names(cj), c("link", "ncat", "coding", "f", "x", "fx", "xx", "ff"))]
    a   <- c(mod, common)
    if (identical(cj$link, "identity") && !is.null(sim$sigma)) a$sigma <- as.numeric(sim$sigma)
    raw <- sprintf("theta = theta_true_%d", j)
    out <- c(out, "", sprintf("## objective %d%s", j,
                              if (is.null(cj$name)) "" else paste0(": ", cj$name)))
    sc <- .ui_fmt_call("simulate_design", a,
                       first = c("support = res$support", "counts = res$counts"), raw = raw)
    out <- c(out, paste0(sprintf("sim_%d <- ", j), sc[1]), sc[-1])
    rows <- sprintf("`This design` = sim_%d$mse", j)
    for (nm in names(sim$designs)) {
      var <- vars[[nm]]
      sc2 <- .ui_fmt_call("simulate_design", a,
                          first = c(sprintf("support = %s_support", var),
                                    sprintf("counts = %s_counts", var)), raw = raw)
      out  <- c(out, paste0(sprintf("sim_%s_%d <- ", var, j), sc2[1]), sc2[-1])
      rows <- c(rows, sprintf("%s = sim_%s_%d$mse", .ui_bt(nm), var, j))
    }
    out <- c(out, sprintf("mse_%d <- rbind(%s)", j, paste(rows, collapse = ", ")),
             sprintf("mse_%d   # per-parameter mean squared errors, objective %d", j, j))
  }
  out
}

# ---- reference values: what they depend on --------------------------------
# Psi*_j is the value the design optimal for component j ALONE attains, so
# compound_design() pays one full multistage solve per component before it can
# even start -- the dominant cost of a compound run.  It depends on the
# components, the design region and stage structure, and any existing design;
# it does NOT depend on alpha (each reference solve runs at weight 1), on the
# efficiency switch, on n, or on the exact-design seed.
#
# This returns the part of a .uic_solver_args() list that Psi* depends on, so
# the app can hold on to a computed psi_star and reuse it -- identical() on
# this key is exactly the condition under which the reference solves would
# reproduce the values it already has.
.uic_psi_key <- function(args) {
  volatile <- c("alpha", "efficiency", "n", "seed", "psi_star", "reference_bound")
  args[setdiff(names(args), volatile)]
}

# ---- results helpers ------------------------------------------------------
# The per-component summary table shown on the results step (and the same
# information print.compound_design() puts on the console).
.uic_summary_table <- function(res) {
  # criterion values on the scale the single-criterion tab reports
  # (D: log det Sigma / v, A: tr(Sigma) / v; smaller is better), so the two
  # tabs agree number for number.  Psi itself stays in the result object.
  d <- data.frame(alpha = round(as.numeric(res$alpha), 4),
                  type = ifelse(res$p == 0, "D", "A"),
                  criterion = signif(as.numeric(res$component_criterion), 7),
                  stringsAsFactors = FALSE)
  if (isTRUE(res$efficiency_weighted)) {
    d$optimal    <- signif(as.numeric(res$component_criterion_star), 7)
    # a guaranteed lower bound vs the TRUE optimum (Becker & Yang, Thm 4.5 / 4.6)
    d$efficiency_lower_bound <- round(as.numeric(res$efficiency_lower_bound), 6)
  }
  rownames(d) <- names(res$psi)
  d
}

# The cross-efficiency table as a data.frame ready for rendering, or NULL when
# no reference solves were run (raw weighting, or supplied psi_star).
.uic_cross_table <- function(res) {
  ce <- res$cross_efficiency
  if (is.null(ce)) return(NULL)
  as.data.frame(round(ce, 4))
}

# Number of design points a compound verify will scan -- the shared grid at the
# finest step.  Closed form; no grid is built.
.uic_verify_points <- function(specs) .ui_verify_points(specs[[1]])

# ---- simulation study ------------------------------------------------------
# A compound design serves several models, so "how well does it do?" has one
# answer PER MODEL.  For each component the study is exactly the classical one
# -- generate the response under that component's model at the design, fit the
# SAME model, repeat -- so this is .ui_simulate() run once per component and
# the results collected.  Nothing about the compound criterion enters: each
# component is judged on its own terms, which is what makes the per-component
# MSEs comparable across competing designs.
#
#   specs, comps : as elsewhere
#   support,counts: the design being studied (integer runs)
#   true_theta   : list of the TRUE parameter values, one vector per component
#                  (the assumed values are only a default; the study may use
#                  different ones)
#   sigma        : residual SD for identity-link components
#   existing     : an existing design to pool in, as in .ui_simulate()
#
# Returns a list of length J; element j is .ui_simulate()'s output for
# component j, or the error condition if that component could not be simulated
# (one failing component must not lose the others).
.uic_simulate <- function(specs, comps, support, counts, true_theta = NULL,
                          sigma = 1, existing = NULL, nsim = 1000L,
                          seed = NULL) {
  J <- length(specs)
  lapply(seq_len(J), function(j) {
    th <- if (!is.null(true_theta) && !is.null(true_theta[[j]])) true_theta[[j]]
          else comps[[j]]$theta
    tryCatch(
      .ui_simulate(specs[[j]], theta = th, sigma = sigma,
                   support = support, counts = counts,
                   existing = existing, obs = NULL,
                   nsim = nsim, seed = seed),
      error = function(e) e)
  })
}

# Per-parameter error summary of one component's simulation, as a named vector;
# NULL when that component failed.  `what` selects the mean squared error
# ("mse", the usual summary) or the median squared error ("medse", which a
# handful of separated replicates cannot distort).
.uic_sim_mse <- function(s, what = c("mse", "medse")) {
  what <- match.arg(what)
  if (is.null(s) || inherits(s, "error") || is.null(s[[what]])) return(NULL)
  stats::setNames(as.numeric(s[[what]]), s$coef_names)
}
