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
                             existing = NULL, psi_star = NULL, n = NULL) {
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
    d$efficiency <- round(as.numeric(res$efficiency), 6)
    # guaranteed lower bound vs the TRUE optimum (Becker & Yang, Thm 4.5 / 4.6)
    d$certified  <- round(as.numeric(res$efficiency_certified), 6)
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
