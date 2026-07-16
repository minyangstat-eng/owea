# ===========================================================================
# app.R -- point-and-click Shiny front end for owea.
#
# Practitioners describe the model with plain inputs (covariate names, types
# and ranges; interactions picked by name; a criterion; parameter values) and
# never touch the package's internal conventions (design_box factor coding,
# within-kind indices, two-digit interaction codes).  The pure helpers below
# translate that friendly description into the arguments optimal_design() /
# exact_design() expect, and are unit-tested independently of the app.
# ===========================================================================

# ---- UI covariate list -> design_box + step_sequence ----------------------
# `covariates` is a list, one per covariate in display order, each a list with
#   name  : character label
#   type  : "continuous" or "factor"
#   lo,hi : numeric range        (continuous)
#   step  : numeric finest grid step (continuous)
#   nlevels: integer >= 2        (factor)
#   steps  : numeric vector, the multistage step SEQUENCE for a continuous
#            covariate (any order; a scalar 'step' is accepted for back-compat).
#            The app's "grid step(s)" text box feeds this via .ui_parse_steps().
# Returns a NAMED design_box, the coarse-to-fine step_sequence (list form, one
# per-continuous-covariate vector per stage; numeric(0) if all-factor), and the
# per-continuous finest steps ('finest').  Each continuous covariate's sequence
# is sorted coarsest-first; unequal-length sequences are padded at the coarse end
# with that covariate's coarsest step so every stage covers all covariates.
.ui_design_box <- function(covariates) {
  db <- vector("list", length(covariates))
  nm <- character(length(covariates))
  seqs <- list()                                   # per-continuous step sequences
  for (i in seq_along(covariates)) {
    cv <- covariates[[i]]
    nm[i] <- cv$name
    if (identical(cv$type, "factor")) {
      L <- as.integer(cv$nlevels)
      if (is.na(L) || L < 2L)
        stop(sprintf("factor '%s' must have at least 2 levels.", cv$name),
             call. = FALSE)
      db[[i]] <- L
    } else {
      if (!is.finite(cv$lo) || !is.finite(cv$hi) || cv$lo >= cv$hi)
        stop(sprintf("continuous '%s' needs low < high.", cv$name), call. = FALSE)
      s <- as.numeric(if (!is.null(cv$steps)) cv$steps else cv$step)
      s <- s[is.finite(s)]
      if (!length(s) || any(s <= 0))
        stop(sprintf("continuous '%s' needs a step sequence of positive numbers.",
                     cv$name), call. = FALSE)
      db[[i]] <- c(cv$lo, cv$hi)
      seqs[[length(seqs) + 1L]] <- sort(unique(s), decreasing = TRUE)  # coarsest first
    }
  }
  names(db) <- nm
  if (length(seqs) == 0L) {
    step_sequence <- numeric(0); finest <- numeric(0)
  } else {
    K <- max(vapply(seqs, length, integer(1)))
    seqs <- lapply(seqs, function(s)              # pad at the coarse end
                   c(rep(s[1], K - length(s)), s))
    step_sequence <- lapply(seq_len(K), function(j) vapply(seqs, `[`, numeric(1), j))
    finest <- vapply(seqs, function(s) s[length(s)], numeric(1))
  }
  list(design_box = db, step_sequence = step_sequence, finest = finest)
}

# Parse the "grid step(s)" text box: one number or a comma-separated step
# sequence ("0.5, 0.1, 0.02").  Numeric input passes through unchanged (for
# programmatic use).  An empty field or ANY unparseable token gives numeric(0),
# so .ui_design_box() raises its "needs a step sequence of positive numbers."
# error instead of silently dropping a typo.
.ui_parse_steps <- function(txt) {
  if (is.numeric(txt)) return(as.numeric(txt))
  tok <- trimws(strsplit(paste(as.character(txt), collapse = ","), ",")[[1]])
  tok <- tok[nzchar(tok)]
  if (!length(tok)) return(numeric(0))
  s <- suppressWarnings(as.numeric(tok))
  if (any(is.na(s))) return(numeric(0))
  s
}

# Full-box grid size at each stage of the step sequence, for the large-grid
# warning: prod(continuous: round((hi-lo)/step)+1) * prod(factor nlevels).
# Returns one count per stage (integer-ish numeric), or 1 for an all-factor box.
.ui_grid_sizes <- function(covariates) {
  db <- .ui_design_box(covariates)
  fac <- vapply(covariates, function(cv) identical(cv$type, "factor"), logical(1))
  nlev <- prod(vapply(covariates[fac],
                      function(cv) as.numeric(cv$nlevels), numeric(1)))
  if (length(db$step_sequence) == 0L) return(nlev)               # all-factor
  rng <- vapply(covariates[!fac],
                function(cv) as.numeric(cv$hi) - as.numeric(cv$lo), numeric(1))
  vapply(db$step_sequence, function(by) nlev * prod(round(rng / by) + 1), numeric(1))
}

# Number of design points verify_optimality() will scan for a model spec:
# the full design_box grid at the FINEST step of the step sequence (the step
# .ui_solver_args() passes for the "verify" target), times all factor levels.
# Closed form -- no grid is built.
.ui_verify_points <- function(spec) {
  db  <- spec$design_box
  fac <- vapply(db, function(b) length(b) == 1L, logical(1))
  nlev <- if (any(fac)) prod(vapply(db[fac], as.numeric, numeric(1))) else 1
  if (!length(spec$finest)) return(nlev)                         # all-factor
  rng <- vapply(db[!fac], function(b) b[2] - b[1], numeric(1))
  nlev * prod(round(rng / spec$finest) + 1)
}

# within-kind index + kind ("f"/"x") for each covariate, in display order.
.ui_kind_index <- function(covariates) {
  fi <- 0L; xi <- 0L
  lapply(covariates, function(cv) {
    if (identical(cv$type, "factor")) { fi <<- fi + 1L; list(kind = "f", idx = fi) }
    else                              { xi <<- xi + 1L; list(kind = "x", idx = xi) }
  })
}

# Translate the friendly UI state into the full formula-style model spec.
#
#   covariates   : as in .ui_design_box()
#   interactions : list of integer pairs c(a, b) -- covariate positions (1-based)
#   quadratics   : integer vector of continuous-covariate positions to square
#   link, ncat   : model family
#
# Returns a list ready to splice into optimal_design()/exact_design():
#   design_box, step_sequence, link, ncat, f, x, fx, xx, ff
# fx / xx / ff are emitted in the LIST-of-pairs form (robust to > 9 covariates).
.ui_model_spec <- function(covariates, interactions = list(), quadratics = integer(0),
                           link = "identity", ncat = NULL,
                           coding = "zero-sum") {
  if (length(covariates) == 0L)
    stop("add at least one covariate.", call. = FALSE)
  coding <- match.arg(coding, c("zero-sum", "baseline"))
  db <- .ui_design_box(covariates)
  ki <- .ui_kind_index(covariates)

  is_factor <- vapply(covariates, function(cv) identical(cv$type, "factor"), logical(1))
  f <- unlist(lapply(which(is_factor),  function(i) ki[[i]]$idx))
  x <- unlist(lapply(which(!is_factor), function(i) ki[[i]]$idx))

  fx <- list(); xx <- list(); ff <- list()
  for (pr in interactions) {
    a <- pr[1]; b <- pr[2]
    ka <- ki[[a]]; kb <- ki[[b]]
    if (ka$kind == "f" && kb$kind == "f") {
      if (ka$idx == kb$idx)
        stop("a factor cannot interact with itself.", call. = FALSE)
      ff[[length(ff) + 1L]] <- c(ka$idx, kb$idx)
    } else if (ka$kind == "x" && kb$kind == "x") {
      xx[[length(xx) + 1L]] <- c(ka$idx, kb$idx)
    } else {
      fac <- if (ka$kind == "f") ka else kb
      con <- if (ka$kind == "x") ka else kb
      fx[[length(fx) + 1L]] <- c(fac$idx, con$idx)
    }
  }
  for (q in quadratics) {
    if (ki[[q]]$kind != "x")
      stop("only continuous covariates can have a quadratic term.", call. = FALSE)
    xx[[length(xx) + 1L]] <- c(ki[[q]]$idx, ki[[q]]$idx)
  }

  nz <- function(v) if (length(v) == 0L) NULL else v
  list(design_box = db$design_box, step_sequence = db$step_sequence,
       finest = db$finest, link = link, ncat = ncat, coding = coding,
       f = nz(f), x = nz(x), fx = nz(fx), xx = nz(xx), ff = nz(ff))
}

# The spec's factor coding, defaulting for specs built before it existed.
.ui_coding <- function(spec)
  if (is.null(spec$coding)) "zero-sum" else spec$coding

# Coefficient labels for the spec (order matches the theta the solvers expect).
# Uses model_info_matrix(), which carries "coef_names" for every link and needs
# no theta.  (The labels do not depend on the coding -- only their meaning does
# -- but we pass it so the spec drives every model call.)
.ui_coef_names <- function(spec) {
  mf <- model_info_matrix(design_box = spec$design_box, link = spec$link,
                          f = spec$f, x = spec$x, fx = spec$fx, xx = spec$xx,
                          ff = spec$ff, ncat = spec$ncat,
                          coding = .ui_coding(spec))
  attr(mf, "coef_names")
}

# ---- wizard step graph ----------------------------------------------------
# The ordered step ids the app walks through.  Steps that do not apply are
# simply absent: the theta step for the identity link (no local design), and the
# existing-design / existing-dataset steps unless the user said they have one.
.ui_wizard_steps <- function(link = "identity", start = "none") {
  start <- match.arg(start, c("none", "design", "data"))
  c("model", "start",
    if (identical(start, "design")) "design_in",
    if (identical(start, "data"))   "data_in",
    if (!identical(link, "identity")) "theta",
    "criterion", "design_type", "review", "results")
}

# ---- the single argument assembler ----------------------------------------
# Every solver call the app makes goes through here, so the model (and above all
# the factor `coding`) can never drift between the design, the fit, the
# simulation and the verification.
#
#   spec     : .ui_model_spec() output
#   target   : which function the args are for
#   theta    : assumed parameters (NULL/ignored for the identity link)
#   p        : 0 = D, 1 = A
#   subset   : integer indices of the parameters of interest (NULL = all)
#   existing : NULL, or .ui_existing_from_csv()/.ui_existing_from_data() output
#              -- list(points, weights, n0, n1)
#   n        : exact designs only, the runs to allocate
#
# Returns a named list to do.call() into.  With no existing design the xi0_* /
# n0 / n1 arguments are OMITTED entirely (passing xi0_points with n0 = 0 would
# make the solver warn that the existing design is ignored).
.ui_solver_args <- function(spec, target = c("optimal", "exact", "verify",
                                             "fit", "simulate"),
                            theta = NULL, p = 0L, subset = NULL,
                            existing = NULL, n = NULL) {
  target <- match.arg(target)
  cn <- .ui_coef_names(spec)
  k  <- length(cn)

  # model block -- identical for every target
  args <- list(design_box = spec$design_box, link = spec$link, ncat = spec$ncat,
               f = spec$f, x = spec$x, fx = spec$fx, xx = spec$xx, ff = spec$ff,
               coding = .ui_coding(spec))

  # theta: NULL for the identity link, required and length-checked otherwise
  if (!identical(target, "fit")) {
    if (identical(spec$link, "identity")) {
      args$theta <- NULL
    } else {
      if (is.null(theta) || length(theta) != k || any(!is.finite(theta)))
        stop(sprintf("the %s model needs %d finite assumed parameter value(s).",
                     spec$link, k), call. = FALSE)
      args$theta <- as.numeric(theta)
    }
  }

  # criterion + parameters of interest
  if (target %in% c("optimal", "exact", "verify")) {
    args$p <- .check_criterion(p)
    if (length(subset)) {
      s <- sort(unique(as.integer(subset)))
      if (any(is.na(s)) || any(s < 1L) || any(s > k))
        stop(sprintf("the parameters of interest must be among 1..%d.", k),
             call. = FALSE)
      args$subset <- s                      # an EMPTY subset stays NULL (= all)
    }
  }

  # the design space
  if (target %in% c("optimal", "exact")) args$step_sequence <- spec$step_sequence
  if (identical(target, "verify")) {
    args$step       <- spec$finest
    args$max_points <- Inf
  }

  # the existing design
  if (!is.null(existing) && target %in% c("optimal", "exact", "verify")) {
    e <- .ui_check_existing(existing, spec)
    args$xi0_points  <- e$points
    args$xi0_weights <- e$weights
    args$n0          <- e$n0
    args$n1          <- e$n1
  }

  # exact designs allocate n runs.  With an existing design we TIE n = n1: the
  # criterion is then reported for exactly the balance that is realised, and the
  # total-sample information really is (n0 + n) times the per-sample one (see
  # the note in the app's Design tab).
  if (identical(target, "exact")) {
    if (is.null(n) || !is.finite(n) || n < 1)
      stop("an exact design needs the number of runs n (>= 1).", call. = FALSE)
    args$n <- as.integer(round(n))
    if (!is.null(args$n1)) args$n1 <- args$n
  }
  args
}

# Validate an existing design against the model spec.
.ui_check_existing <- function(existing, spec) {
  pts <- as.matrix(existing$points); storage.mode(pts) <- "double"
  w   <- as.numeric(existing$weights)
  n0  <- existing$n0
  n1  <- if (is.null(existing$n1)) 1 else existing$n1
  d   <- length(spec$design_box)
  if (nrow(pts) != length(w))
    stop("the existing design has ", nrow(pts), " point(s) but ", length(w),
         " weight(s).", call. = FALSE)
  if (ncol(pts) != d)
    stop(sprintf(paste0("the existing design has %d column(s) but the model ",
                        "has %d covariate(s)."), ncol(pts), d), call. = FALSE)
  if (is.null(n0) || !is.finite(n0) || n0 < 1)
    stop("the existing design needs a sample size n0 >= 1.", call. = FALSE)
  if (!is.finite(n1) || n1 < 1)
    stop("the new stage needs a sample size n1 >= 1.", call. = FALSE)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0 || any(w < 0))
    stop("the existing design's weights must be nonnegative and sum to a ",
         "positive value.", call. = FALSE)
  meta <- .parse_design_box(spec$design_box)
  .validate_factor_columns(pts, meta$is_factor, meta$nlevels,
                           "the existing design")
  list(points = pts, weights = w / sw, n0 = as.integer(round(n0)),
       n1 = as.integer(round(n1)))
}

# ---- assumed parameter values ---------------------------------------------
# The model "plan" behind a spec (parameter count, per-link parameterisation).
.ui_plan <- function(spec) {
  meta <- .parse_design_box(spec$design_box)
  .build_model_terms(meta, spec$f, spec$x, spec$fx, spec$xx, spec$ff,
                     intercept = TRUE, coding = .ui_coding(spec),
                     link = spec$link, cov_names = names(spec$design_box),
                     ncat = spec$ncat)
}

# A draw from the standard normal, in the spec's parameterisation.  NULL for the
# identity link (which needs no assumed values).  Reuses .random_theta(), which
# keeps the cumulative model's thresholds strictly increasing.
.ui_random_theta <- function(spec) {
  if (identical(spec$link, "identity")) return(NULL)
  .random_theta(.ui_plan(spec))
}

# Check hand-edited assumed values BEFORE they reach the solver, so the wizard
# can show an inline message instead of a hard error.  NULL when they are fine.
.ui_check_theta <- function(theta, spec) {
  if (identical(spec$link, "identity")) return(NULL)
  k <- length(.ui_coef_names(spec))
  if (is.null(theta) || length(theta) != k)
    return(sprintf("this model needs %d parameter value(s); got %d.", k,
                   length(theta)))
  if (any(!is.finite(theta)))
    return("every parameter value must be a finite number.")
  if (identical(spec$link, "cumulative")) {
    a <- theta[seq_len(as.integer(spec$ncat) - 1L)]
    if (length(a) > 1L && any(diff(a) <= 0))
      return(sprintf(paste0("the first %d value(s) are the thresholds of the ",
                            "ordinal model and must be strictly increasing."),
                     length(a)))
  }
  NULL
}

# ---- existing designs ------------------------------------------------------
# A design pasted or uploaded as CSV (the format the app downloads): covariate
# columns plus a 'count' or a 'weight' column.  Counts carry a sample size, so
# they give n0 for free; weights do not, so n0 must be supplied.
.ui_existing_from_csv <- function(support, val, valcol, n0 = NULL) {
  valcol <- match.arg(valcol, c("count", "weight"))
  pts <- as.matrix(support); storage.mode(pts) <- "double"
  val <- as.numeric(val)
  notes <- character(0)
  if (!nrow(pts) || any(!is.finite(val)) || any(val < 0))
    stop("the existing design needs nonnegative ", valcol, "s.", call. = FALSE)
  if (identical(valcol, "count")) {
    if (any(abs(val - round(val)) > 1e-8))
      stop("the 'count' column must contain whole numbers.", call. = FALSE)
    if (is.null(n0)) {
      n0 <- sum(round(val))
      notes <- c(notes, sprintf("n0 = %d taken from the counts.", as.integer(n0)))
    }
  } else if (abs(sum(val) - 1) > 1e-6) {
    notes <- c(notes, sprintf("weights summed to %.6f; rescaled to sum to 1.",
                              sum(val)))
  }
  keep <- val > 0                        # zero-weight rows contribute nothing
  pts <- pts[keep, , drop = FALSE]; val <- val[keep]
  sw <- sum(val)
  if (!nrow(pts) || !is.finite(sw) || sw <= 0)
    stop("the existing design has no runs.", call. = FALSE)
  list(points = pts, weights = val / sw, n0 = n0, notes = notes)
}

# An existing DATA SET (covariates + a response) reused as an existing design:
# aggregate the observed covariate rows into a support with proportion weights.
# n0 is the number of observations.  The covariate columns are split exactly as
# fit_design() splits them, so the design and the estimates see the same data.
.ui_existing_from_data <- function(data, design_box, response = NULL) {
  sp <- .split_fit_data(data, response)
  X  <- sp$X
  if (ncol(X) != length(design_box))
    stop(sprintf(paste0("the data set has %d covariate column(s) but the model ",
                        "has %d covariate(s)."), ncol(X), length(design_box)),
         call. = FALSE)
  meta <- .parse_design_box(design_box)
  .validate_factor_columns(X, meta$is_factor, meta$nlevels,
                           "the data set's covariates")
  keys <- do.call(paste, c(as.data.frame(X), sep = "\r"))
  uk   <- unique(keys)
  pts  <- X[match(uk, keys), , drop = FALSE]
  cnt  <- as.integer(table(factor(keys, levels = uk)))
  colnames(pts) <- names(design_box)
  df <- as.data.frame(pts); df[["count"]] <- cnt
  list(points = pts, weights = cnt / sum(cnt), n0 = nrow(X),
       counts = cnt, df = df, response = sp$resp_name)
}

# ---- simulation study, POOLED over both stages -----------------------------
# The design is optimised for the COMBINED information of the existing and the
# new runs, so the simulation must estimate theta the way the experiment
# actually will: from both stages at once.  Simulating the new runs alone would
# understate the precision of every design being compared.
#
#   support, counts : the NEW runs (the design under study)
#   existing        : NULL, or list(points, weights, n0, counts) -- the first
#                     stage.  Its responses are unknown, so they are simulated
#                     from `theta` in every replicate.
#   obs             : NULL, or list(data, response) -- an OBSERVED data set used
#                     as the first stage.  Its responses are real, so they are
#                     held FIXED across replicates and only the new runs are
#                     simulated.  Takes precedence over `existing`.
#
# Returns the fields of simulate_design(nsim > 1) that the app consumes, plus
# n_existing / n_new.
.ui_simulate <- function(spec, theta, sigma = 1, support, counts,
                         existing = NULL, obs = NULL, nsim = 1000L, seed = NULL) {
  plan <- .ui_plan(spec)
  k    <- plan$k
  # every link needs the TRUE values here -- even the identity model, whose
  # DESIGN does not depend on them but whose simulated responses do
  if (is.null(theta) || length(theta) != k || any(!is.finite(theta)))
    stop(sprintf("simulation needs %d finite true parameter value(s).", k),
         call. = FALSE)
  theta <- as.numeric(theta)

  rows <- function(X, cnt) {                    # expand a design into runs
    X <- as.matrix(X); storage.mode(X) <- "double"
    Xr <- X[rep(seq_len(nrow(X)), as.integer(cnt)), , drop = FALSE]
    t(apply(Xr, 1L, .model_row, plan = plan))   # -> runs x base_k
  }
  # .model_row returns a length-base_k vector; apply()+t() keeps it a matrix
  fix_dim <- function(M) if (plan$base_k == 1L) matrix(as.numeric(M), ncol = 1L) else M

  Xnew <- fix_dim(rows(support, counts))
  n_new <- nrow(Xnew)

  Xold <- NULL; yold <- NULL
  if (!is.null(obs)) {                          # observed first stage: y is real
    sp   <- .split_fit_data(obs$data, obs$response)
    Xold <- fix_dim(t(apply(sp$X, 1L, .model_row, plan = plan)))
    yold <- .check_response(sp$y, plan, sp$resp_name)
  } else if (!is.null(existing)) {              # design-only first stage: simulate y
    cnt <- existing$counts
    if (is.null(cnt) || sum(cnt) != existing$n0)
      cnt <- .apportion(existing$weights, as.integer(existing$n0))
    Xold <- fix_dim(rows(existing$points, cnt))
  }
  n_old <- if (is.null(Xold)) 0L else nrow(Xold)
  X     <- if (is.null(Xold)) Xnew else rbind(Xold, Xnew)

  if (!is.null(seed)) set.seed(as.integer(seed))
  E <- matrix(NA_real_, nsim, k, dimnames = list(NULL, plan$coef_names))
  conv <- logical(nsim)
  for (s in seq_len(nsim)) {
    y <- c(if (!is.null(yold)) yold                      # observed: held fixed
           else if (n_old) .simulate_y(plan, Xold, theta, sigma),
           .simulate_y(plan, Xnew, theta, sigma))
    fr <- tryCatch(.fit_model(plan, X, y, sigma, want_vcov = FALSE),
                   error = function(e) NULL)
    if (!is.null(fr) && isTRUE(fr$converged)) { E[s, ] <- fr$theta_hat; conv[s] <- TRUE }
  }
  ok <- which(conv)
  if (!length(ok)) stop("no simulation replicate converged.", call. = FALSE)
  list(estimates = E, theta = theta, coef_names = plan$coef_names,
       link = plan$link, N = nrow(X), n_new = n_new, n_existing = n_old,
       nsim = nsim, n_converged = length(ok),
       theta_hat_mean = colMeans(E[ok, , drop = FALSE]),
       bias = colMeans(E[ok, , drop = FALSE]) - theta,
       mse  = colMeans(sweep(E[ok, , drop = FALSE], 2, theta, "-")^2),
       se_empirical = apply(E[ok, , drop = FALSE], 2, stats::sd))
}

# ---- efficiency of a derived design under a DIFFERENT criterion ------------
# Efficiency only means anything against the design that is optimal FOR that
# criterion, so this re-solves for the reference -- with the EXACT same
# step_sequence the original criterion used (it is carried in `spec`).  The
# design's own criterion needs no grid at all (criterion_only), so the
# reference solve is the only grid search in this computation.  Both criteria
# are evaluated with the same existing design, so they live on the same scale.
#   D (p = 0): the criterion is a log       -> eff = exp(crit_ref - crit_design)
#   A (p = 1): the criterion is an average  -> eff = crit_ref / crit_design
# Both are <= 1 because the criterion is minimised.
.ui_efficiency <- function(res, spec, theta = NULL, p_new = 0L, subset_new = NULL,
                           existing = NULL) {
  p_new <- .check_criterion(p_new)
  v <- do.call(verify_optimality,
               c(list(support = res$support, weights = res$weights,
                      criterion_only = TRUE),
                 .ui_solver_args(spec, "verify", theta, p_new, subset_new,
                                 existing)))
  ref <- do.call(optimal_design,
                 .ui_solver_args(spec, "optimal", theta, p_new, subset_new,
                                 existing))
  crit_design <- v$criterion
  crit_ref    <- ref$criterion
  eff <- if (!is.finite(crit_design) || !is.finite(crit_ref)) NA_real_
         else if (p_new == 0L) exp(crit_ref - crit_design)
         else                  crit_ref / crit_design
  list(p = p_new, crit_design = crit_design, crit_ref = crit_ref,
       efficiency = eff, converged = isTRUE(ref$converged),
       max_sensitivity = v$max_sensitivity)
}

#' Launch the owea point-and-click web app.
#'
#' Opens a Shiny application that builds D- and A-optimal designs from
#' point-and-click inputs -- no R code required. It walks through the model
#' (family, covariates, interactions, and the coding of any factor levels),
#' asks whether you already have a design or a data set (either can be reused
#' as a first stage, and a data set can also supply the assumed parameter
#' values via \code{\link{fit_design}}), then the assumed parameters -- typed in
#' or drawn from N(0,1) -- the criterion, the parameters of interest and the
#' design type. It returns the design, its plot, the information matrix, a
#' downloadable table, and can report the design's efficiency under a different
#' criterion. Steps that do not apply are skipped: the linear (normal) model
#' needs no assumed parameters, since its design is not local. This is a
#' convenience wrapper over \code{\link{optimal_design}} /
#' \code{\link{exact_design}}.
#'
#' @param launch.browser open the app in the system browser (default
#'   \code{TRUE}).
#' @param ... further arguments passed to \code{shiny::runApp}.
#' @return Called for its side effect (runs the app); returns nothing useful.
#' @seealso \code{\link{optimal_design}}, \code{\link{exact_design}}.
#' @export
run_owea_app <- function(launch.browser = TRUE, ...) {
  need <- c("shiny", "DT")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("the web app needs package(s) ", paste(miss, collapse = ", "),
         ".\n  install.packages(c(", paste(sprintf('\"%s\"', miss), collapse = ", "),
         "))", call. = FALSE)
  app_dir <- system.file("shiny", "owea-app", package = "owea")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R")))
    stop("could not locate the bundled Shiny app; reinstall 'owea'.", call. = FALSE)
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
