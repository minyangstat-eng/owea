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
# With continuous = TRUE (the app's "search the continuous region" choice) no
# steps are needed or validated: step_sequence is NULL, finest is numeric(0) and
# the returned list carries continuous = TRUE for .ui_solver_args().
.ui_design_box <- function(covariates, continuous = FALSE) {
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
      db[[i]] <- c(cv$lo, cv$hi)
      if (!isTRUE(continuous)) {
        s <- as.numeric(if (!is.null(cv$steps)) cv$steps else cv$step)
        s <- s[is.finite(s)]
        if (!length(s) || any(s <= 0))
          stop(sprintf("continuous '%s' needs a step sequence of positive numbers.",
                       cv$name), call. = FALSE)
        seqs[[length(seqs) + 1L]] <- sort(unique(s), decreasing = TRUE)  # coarsest first
      }
    }
  }
  names(db) <- nm
  if (isTRUE(continuous))
    return(list(design_box = db, step_sequence = NULL, finest = numeric(0),
                continuous = TRUE))
  if (length(seqs) == 0L) {
    step_sequence <- numeric(0); finest <- numeric(0)
  } else {
    K <- max(vapply(seqs, length, integer(1)))
    seqs <- lapply(seqs, function(s)              # pad at the coarse end
                   c(rep(s[1], K - length(s)), s))
    step_sequence <- lapply(seq_len(K), function(j) vapply(seqs, `[`, numeric(1), j))
    finest <- vapply(seqs, function(s) s[length(s)], numeric(1))
  }
  list(design_box = db, step_sequence = step_sequence, finest = finest,
       continuous = FALSE)
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
# With continuous = TRUE there is no grid to warn about: the single count is the
# small coarse grid the continuous search starts from (3 levels per continuous
# covariate times the factor levels).
.ui_grid_sizes <- function(covariates, continuous = FALSE) {
  db <- .ui_design_box(covariates, continuous = continuous)
  fac <- vapply(covariates, function(cv) identical(cv$type, "factor"), logical(1))
  nlev <- prod(vapply(covariates[fac],
                      function(cv) as.numeric(cv$nlevels), numeric(1)))
  if (isTRUE(continuous)) return(nlev * 3^sum(!fac))
  if (length(db$step_sequence) == 0L) return(nlev)               # all-factor
  rng <- vapply(covariates[!fac],
                function(cv) as.numeric(cv$hi) - as.numeric(cv$lo), numeric(1))
  vapply(db$step_sequence, function(by) nlev * prod(round(rng / by) + 1), numeric(1))
}

# Number of design points verify_optimality() will scan for a model spec:
# the full design_box grid at the FINEST step of the step sequence (the step
# .ui_solver_args() passes for the "verify" target), times all factor levels.
# Closed form -- no grid is built.  NA for a continuous-search spec: its
# verification searches the region (multi-start + audit) instead of a grid --
# unless `step` is given, the grid step(s) of an explicit grid check (a scalar
# or one value per continuous covariate), which is then counted for any spec.
.ui_verify_points <- function(spec, step = NULL) {
  db  <- spec$design_box
  fac <- vapply(db, function(b) length(b) == 1L, logical(1))
  nlev <- if (any(fac)) prod(vapply(db[fac], as.numeric, numeric(1))) else 1
  if (!is.null(step)) {
    step <- as.numeric(step)
    if (!length(step) || any(!is.finite(step)) || any(step <= 0)) return(NA_real_)
    if (!any(!fac)) return(nlev)
    rng <- vapply(db[!fac], function(b) b[2] - b[1], numeric(1))
    if (length(step) != 1L && length(step) != length(rng)) return(NA_real_)
    return(nlev * prod(round(rng / step) + 1))
  }
  if (isTRUE(spec$continuous)) return(NA_real_)
  if (!length(spec$finest)) return(nlev)                         # all-factor
  rng <- vapply(db[!fac], function(b) b[2] - b[1], numeric(1))
  nlev * prod(round(rng / spec$finest) + 1)
}

# A suggested grid step for an explicit grid check of a design, as the text the
# app's step box shows: the finest step(s) of a grid spec, or, for a continuous
# search, one fortieth of each continuous covariate's range (2 significant
# digits) -- one value per continuous covariate, comma-separated.
.ui_verify_step_default <- function(spec) {
  if (length(spec$finest)) return(paste(format(spec$finest), collapse = ", "))
  db  <- spec$design_box
  fac <- vapply(db, function(b) length(b) == 1L, logical(1))
  if (!any(!fac)) return("")
  rng <- vapply(db[!fac], function(b) b[2] - b[1], numeric(1))
  paste(format(signif(rng / 40, 2)), collapse = ", ")
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
#   continuous   : TRUE = search the continuous covariates in the continuous
#                  region (optimal_design(continuous = TRUE); no grid steps)
#   n_audit      : random audit points of the continuous search (0 = none);
#                  anything that is not a nonnegative number falls back to the
#                  solver default 20000
#
# Returns a list ready to splice into optimal_design()/exact_design():
#   design_box, step_sequence, link, ncat, f, x, fx, xx, ff, continuous, n_audit
# fx / xx / ff are emitted in the LIST-of-pairs form (robust to > 9 covariates).
.ui_model_spec <- function(covariates, interactions = list(), quadratics = integer(0),
                           link = "identity", ncat = NULL,
                           coding = "zero-sum", continuous = FALSE,
                           n_audit = 20000L) {
  if (length(covariates) == 0L)
    stop("add at least one covariate.", call. = FALSE)
  coding <- match.arg(coding, c("zero-sum", "baseline"))
  db <- .ui_design_box(covariates, continuous = continuous)
  n_audit <- suppressWarnings(as.integer(round(as.numeric(n_audit)[1])))
  if (is.na(n_audit) || n_audit < 0L) n_audit <- 20000L
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
       f = nz(f), x = nz(x), fx = nz(fx), xx = nz(xx), ff = nz(ff),
       continuous = isTRUE(db$continuous),
       n_audit = if (isTRUE(db$continuous)) n_audit else NA_integer_)
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
#   verify_step : "verify" target only -- when given, the check runs over a
#              GRID of the design box at these step(s) (a scalar, or one value
#              per continuous covariate) whatever the spec's own search mode
#   verify_audit : "verify" target only -- when given, the check is the
#              continuous one: a multi-start search of the region plus a
#              RANDOM AUDIT of this many points (>= 1), whatever the spec's mode
#              (verify_step wins if both are given).  With neither, the check
#              runs over the design space the computation used.
#
# Returns a named list to do.call() into.  With no existing design the xi0_* /
# n0 / n1 arguments are OMITTED entirely (passing xi0_points with n0 = 0 would
# make the solver warn that the existing design is ignored).
.ui_solver_args <- function(spec, target = c("optimal", "exact", "verify",
                                             "fit", "simulate"),
                            theta = NULL, p = 0L, subset = NULL,
                            existing = NULL, n = NULL, verify_step = NULL,
                            verify_audit = NULL) {
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

  # the design space: a grid (step_sequence / the finest step for verification)
  # or, for a continuous-search spec, the continuous region itself
  cont_args <- function(a) {              # the continuous search and its audit size
    a$continuous <- TRUE
    if (!is.null(spec$n_audit) && !is.na(spec$n_audit))
      a$n_audit <- as.integer(spec$n_audit)
    a
  }
  if (target %in% c("optimal", "exact")) {
    if (isTRUE(spec$continuous)) args <- cont_args(args)
    else args$step_sequence <- spec$step_sequence
  }
  if (identical(target, "verify")) {
    if (!is.null(verify_step)) {              # an explicit grid check
      st <- as.numeric(verify_step)
      if (!length(st) || any(!is.finite(st)) || any(st <= 0))
        stop("the grid step(s) of the check must be positive numbers.", call. = FALSE)
      args$step       <- st
      args$max_points <- Inf
    } else if (!is.null(verify_audit)) {      # an explicit random audit
      na <- suppressWarnings(as.integer(round(as.numeric(verify_audit)[1])))
      if (is.na(na) || na < 1L)
        stop("the random audit needs a positive number of points.", call. = FALSE)
      args$continuous <- TRUE
      args$n_audit    <- na
    } else if (isTRUE(spec$continuous)) {
      args <- cont_args(args)
    } else {
      args$step       <- spec$finest
      args$max_points <- Inf
    }
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

# ---- R code that reproduces a computation ----------------------------------
# The app's results page shows (and downloads) the solver call it made as a
# plain R script, so a user can rerun and adapt the computation in R.  `args` is
# the do.call() argument list from .ui_solver_args() (plus seed / n for an exact
# design); `verify_args` the "verify" list for the same spec, appended as a
# verify_optimality() call for approximate designs.

# back-tick a name that is not syntactic (covariate names are typed by the user)
.ui_bt <- function(nm) ifelse(make.names(nm) == nm, nm, paste0("`", nm, "`"))

# one argument value as R source (vectors, named lists, matrices as rbind())
.ui_fmt_arg <- function(v, indent = 4L) {
  pad <- strrep(" ", indent)
  # each number on its own (so c(0.2, 1) is not written as c(0.2, 1.0))
  num <- function(z) vapply(as.numeric(z), function(v) format(v, digits = 15), character(1))
  if (is.null(v)) return("NULL")
  if (is.function(v)) return("<function>")
  if (is.matrix(v)) {
    rows <- apply(v, 1, function(r) paste0("c(", paste(num(r), collapse = ", "), ")"))
    return(paste0("rbind(", paste(rows, collapse = paste0(",\n", pad, "      ")), ")"))
  }
  if (is.list(v)) {
    inner <- vapply(v, .ui_fmt_arg, character(1), indent = indent)
    nm <- names(v)
    if (!is.null(nm) && any(nzchar(nm)))
      inner <- paste0(ifelse(nzchar(nm), paste0(.ui_bt(nm), " = "), ""), inner)
    return(paste0("list(", paste(inner, collapse = ", "), ")"))
  }
  s <- if (is.character(v)) paste0('"', v, '"')
       else if (is.logical(v)) ifelse(v, "TRUE", "FALSE")
       else num(v)
  if (length(v) == 1L) s else paste0("c(", paste(s, collapse = ", "), ")")
}

# the call `fn(` + one argument per line + `)`.  `first` / `raw` are ready-made
# "name = expression" lines placed before / after the formatted arguments.
.ui_fmt_call <- function(fn, args, first = character(0), raw = character(0)) {
  keep <- args[!vapply(args, is.null, logical(1))]
  pref <- c("n", "design_box", "step_sequence", "continuous", "n_audit", "link",
            "ncat", "f", "x", "fx", "xx", "ff", "coding", "theta", "p", "subset",
            "xi0_points", "xi0_weights", "n0", "n1", "seed", "step", "max_points",
            "sigma", "nsim", "existing")
  nm   <- names(keep)
  keep <- keep[c(intersect(pref, nm), setdiff(nm, pref))]
  split_nv <- function(s) list(n = trimws(sub("=.*$", "", s)),
                               v = trimws(sub("^[^=]*=", "", s)))
  f1 <- split_nv(first); r1 <- split_nv(raw)
  w  <- max(nchar(c(names(keep), f1$n, r1$n)), 1L)
  line_of <- function(n, v) if (length(n)) paste0(formatC(n, width = -w), " = ", v)
                            else character(0)
  body <- c(line_of(f1$n, f1$v),
            line_of(names(keep), vapply(keep, .ui_fmt_arg, character(1))),
            line_of(r1$n, r1$v))
  body <- paste0("  ", body, c(rep(",", length(body) - 1L), ""))
  c(paste0(fn, "("), body, ")")
}

# `sim` (exact designs only) describes the app's simulation study: theta (the
# true values), sigma (identity link; else NULL), nsim, seed, existing
# (list(points, counts) or NULL), obs (list(data, response) or NULL) and
# designs -- a named list of list(support, counts), the designs compared with.
.ui_r_code <- function(args, target = c("optimal", "exact"), verify_args = NULL,
                       sim = NULL) {
  target <- match.arg(target)
  exact  <- identical(target, "exact")
  fn     <- if (exact) "exact_design" else "optimal_design"
  ver    <- tryCatch(as.character(utils::packageVersion("owea")), error = function(e) "")
  lines  <- c(
    sprintf("## R code generated by the owea app (owea %s): it reproduces the computed design.", ver),
    sprintf("## Edit any argument and rerun; see ?%s for all options.", fn),
    "library(owea)", "")
  if (isTRUE(args$continuous))
    lines <- c(lines,
               "## The continuous search uses random starts and random audit points; a fixed",
               "## seed reproduces this run exactly (any seed gives the same design up to eps0).",
               "set.seed(1)", "")
  call1 <- .ui_fmt_call(fn, args)
  lines <- c(lines, paste0("res <- ", call1[1]), call1[-1],
             if (exact) "print(res)          # runs per support point, per-sample and total criterion"
             else       "print_result(res)   # support points, weights, criterion, max sensitivity")
  if (!exact && !is.null(verify_args)) {
    vc <- .ui_fmt_call("verify_optimality", verify_args,
                       first = c("support = res$support", "weights = res$weights"))
    lines <- c(lines, "",
               "## Check the design against the general equivalence theorem over the design",
               "## space given below (a grid at 'step', or the continuous region)",
               paste0("v <- ", vc[1]), vc[-1],
               "v$is_optimal$value        # TRUE when the max sensitivity is <= tol",
               "v$max_sensitivity",
               "v$efficiency_lower_bound")
  }
  if (exact && !is.null(sim)) lines <- c(lines, "", .ui_r_code_sim(args, sim))
  paste(lines, collapse = "\n")
}

# The simulation-study block of .ui_r_code(): simulate_design() at the exact
# design and at every design the app compared it with (each pooled with the
# first stage, if any), then the mean squared errors side by side.
.ui_r_code_sim <- function(args, sim) {
  mod <- args[intersect(names(args), c("design_box", "link", "ncat", "f", "x", "fx",
                                        "xx", "ff", "coding"))]
  common <- c(mod, list(nsim = as.integer(sim$nsim), seed = as.integer(sim$seed)))
  if (!is.null(sim$sigma)) common$sigma <- as.numeric(sim$sigma)
  raw <- "theta = theta_true"
  out <- c("## Simulation study: simulate responses at the design, refit the model and report",
           "## the estimation error per parameter (mean and median squared error)",
           paste0("theta_true <- ", .ui_fmt_arg(as.numeric(sim$theta)),
                  "   # the true values used in the app"))
  if (!is.null(sim$obs)) {
    d <- sim$obs$data
    if (is.data.frame(d) && nrow(d) <= 500L) {
      csv <- utils::capture.output(utils::write.csv(d, row.names = FALSE))
      out <- c(out,
               "## the observed first stage: its responses are kept, only the new runs are simulated",
               "obs_data <- read.csv(text = c(",
               paste0("  \"", gsub("\"", "\\\\\"", csv), "\"",
                      c(rep(",", length(csv) - 1L), "")),
               "))")
    } else {
      out <- c(out, sprintf(paste0("obs_data <- read.csv(\"your_data.csv\")   # the observed ",
                                   "data set used in the app (%s rows)"),
                            if (is.data.frame(d)) nrow(d) else "?"))
    }
    raw <- c(raw, sprintf("obs = list(data = obs_data, response = %s)",
                          if (is.null(sim$obs$response)) "NULL"
                          else .ui_fmt_arg(sim$obs$response)))
  } else if (!is.null(sim$existing)) {
    common$existing <- list(points = unname(as.matrix(sim$existing$points)),
                            counts = as.integer(sim$existing$counts))
    out <- c(out, "## the first stage (an existing design) is simulated and analysed together with the new runs")
  }
  sc <- .ui_fmt_call("simulate_design", common,
                     first = c("support = res$support", "counts = res$counts"), raw = raw)
  out <- c(out, paste0("sim <- ", sc[1]), sc[-1], "sim$mse", "sim$medse")
  mse_rows <- "`Exact design` = sim$mse"
  for (nm in names(sim$designs)) {
    dsg <- sim$designs[[nm]]
    var <- tolower(gsub("[^A-Za-z0-9]+", "_", nm))
    out <- c(out, "",
             sprintf("## the %s design the app compared with (%d runs)", nm, sum(dsg$counts)),
             paste0(var, "_support <- ", .ui_fmt_arg(unname(as.matrix(dsg$support)))),
             paste0(var, "_counts  <- ", .ui_fmt_arg(as.integer(dsg$counts))))
    sc2 <- .ui_fmt_call("simulate_design", common,
                        first = c(sprintf("support = %s_support", var),
                                  sprintf("counts = %s_counts", var)), raw = raw)
    out <- c(out, paste0("sim_", var, " <- ", sc2[1]), sc2[-1])
    mse_rows <- c(mse_rows, sprintf("%s = sim_%s$mse", .ui_bt(nm), var))
  }
  if (length(sim$designs))
    out <- c(out, "",
             "## mean squared errors side by side; a ratio above 1 means the exact design has",
             "## the smaller error for that parameter",
             paste0("mse <- rbind(", paste(mse_rows, collapse = ", "), ")"),
             "mse",
             "sweep(mse[-1, , drop = FALSE], 2, mse[1, ], \"/\")")
  out
}

# ---- a simple random sample of n runs from the design space -----------------
# The comparison design of the app's simulation study.  After a grid
# computation the runs are a simple random sample (with replacement) from the
# finest-step grid: the grid is a Cartesian product, so drawing a grid index
# independently for each covariate is exactly uniform over the grid points --
# WITHOUT building the grid, whose size grows exponentially with the number
# of covariates as the step shrinks.  (Drawing a continuous point and rounding
# it to the step is almost the same, except that the two ends of each range
# would get half the probability of an interior point.)  After a continuous
# search there is no grid, so each run is drawn uniformly from the continuous
# region.  Factors are uniform over their levels in both cases.  Returns the
# distinct points with their run counts, and the (implied) pool size (NA for
# the continuous draw).
.ui_srs_design <- function(spec, n) {
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L)
    stop("the simple random sample needs n >= 1.", call. = FALSE)
  meta <- .parse_design_box(spec$design_box)
  d <- length(meta$is_factor)
  on_grid <- !isTRUE(spec$continuous) && length(spec$finest) > 0L
  by <- if (on_grid) .expand_stage_step(as.numeric(spec$finest), meta$is_factor,
                                        sum(!meta$is_factor), d) else NULL
  pts <- matrix(0, n, d); pool_size <- if (on_grid) 1 else NA_real_
  for (j in seq_len(d)) {
    if (meta$is_factor[j]) {
      pts[, j] <- sample.int(meta$nlevels[j], n, replace = TRUE)
      if (on_grid) pool_size <- pool_size * meta$nlevels[j]
    } else if (on_grid) {                       # grid values lo + (0..m-1) * step
      m <- round((meta$hi[j] - meta$lo[j]) / by[j]) + 1
      pts[, j] <- meta$lo[j] + (sample.int(m, n, replace = TRUE) - 1) * by[j]
      pool_size <- pool_size * m
    } else {
      pts[, j] <- stats::runif(n, meta$lo[j], meta$hi[j])
    }
  }
  keys <- do.call(paste, c(as.data.frame(pts), sep = "\r"))
  uk   <- unique(keys)
  sup  <- unname(as.matrix(pts[match(uk, keys), , drop = FALSE]))
  colnames(sup) <- names(spec$design_box)
  list(support = sup, counts = as.integer(table(factor(keys, levels = uk))),
       pool_size = pool_size)
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
                         existing = NULL, obs = NULL, nsim = 1000L, seed = NULL)
  .simulate_pooled(.ui_plan(spec), theta, sigma, support, counts, existing, obs,
                   nsim, seed)

# The same from a model "plan" (.build_model_terms()); also the pooled path of
# the exported simulate_design(existing = , obs = ).
#   existing : list(points, counts) or list(points, weights, n0): a first-stage
#              design whose responses are simulated along with the new runs
#   obs      : list(data, response): an observed first stage whose responses
#              are kept fixed (takes precedence over `existing`)
.simulate_pooled <- function(plan, theta, sigma = 1, support, counts,
                             existing = NULL, obs = NULL, nsim = 1000L, seed = NULL) {
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
    has_w <- !is.null(existing$weights) && !is.null(existing$n0)
    if (is.null(cnt)) {
      if (!has_w)
        stop("'existing' needs 'counts' (runs per point), or 'weights' and 'n0'.",
             call. = FALSE)
      cnt <- .apportion(existing$weights, as.integer(existing$n0))
    } else if (has_w && sum(cnt) != existing$n0) {
      cnt <- .apportion(existing$weights, as.integer(existing$n0))
    }
    Xold <- fix_dim(rows(existing$points, as.integer(round(cnt))))
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
  sq <- sweep(E[ok, , drop = FALSE], 2, theta, "-")^2      # squared errors
  list(estimates = E, theta = theta, coef_names = plan$coef_names,
       link = plan$link, N = nrow(X), n_new = n_new, n_existing = n_old,
       nsim = nsim, n_converged = length(ok),
       theta_hat_mean = colMeans(E[ok, , drop = FALSE]),
       bias = colMeans(E[ok, , drop = FALSE]) - theta,
       mse  = colMeans(sq),
       # The MEDIAN squared error, reported alongside the mean.  For a binary
       # or count response a design concentrated on few distinct points can
       # separate in the odd replicate; the fit still "converges" but the
       # coefficients run away, and one such replicate can dominate a mean.
       # The median says what a typical replicate achieved, so the two together
       # show both the average cost and whether a tail is driving it.
       medse = apply(sq, 2, stats::median),
       se_empirical = apply(E[ok, , drop = FALSE], 2, stats::sd))
}

# Share of a parameter's MSE contributed by its single worst replicate, maximised
# over parameters.  Near 1/n_converged for a well-behaved study; close to 1 when
# one separated replicate is carrying the mean.  Used by both branches of the app
# to decide whether to warn about the MSE column.
.ui_sim_dominance <- function(s) {
  if (is.null(s) || inherits(s, "error") || is.null(s$estimates)) return(NA_real_)
  ok <- stats::complete.cases(s$estimates)
  if (!any(ok)) return(NA_real_)
  sq <- sweep(s$estimates[ok, , drop = FALSE], 2, s$theta, "-")^2
  d <- apply(sq, 2, function(z) {
    tot <- sum(z, na.rm = TRUE)
    if (!is.finite(tot) || tot <= 0) 0 else max(z, na.rm = TRUE) / tot
  })
  suppressWarnings(max(d, na.rm = TRUE))
}

# ---- the simulation comparison, as a matrix -------------------------------
# `slots` is a NAMED list, one entry per design being compared, each a numeric
# per-parameter error vector (mse or medse) or NULL -- a design that was not
# run, or one that failed.  Both branches of the app compare the same way, so
# both build their table and their ratio columns from this one matrix and
# cannot disagree about what was compared.
#
# Returns a designs x parameters matrix, or NULL when nothing comparable is
# left.
.ui_sim_matrix <- function(slots, param_names = NULL) {
  if (!length(slots)) return(NULL)
  nms <- names(slots)
  if (is.null(nms)) nms <- paste("design", seq_along(slots))
  keep <- list(); labs <- character(0); k <- NULL
  for (i in seq_along(slots)) {
    v <- slots[[i]]
    if (is.null(v) || inherits(v, "error") || !length(v)) next
    v <- as.numeric(v)
    if (!any(is.finite(v))) next                 # a design that produced nothing
    if (is.null(k)) k <- length(v) else if (length(v) != k) next
    v[!is.finite(v)] <- NA_real_                 # an unusable parameter only
    keep[[length(keep) + 1L]] <- v
    labs <- c(labs, nms[i])
    if (is.null(param_names) && !is.null(names(slots[[i]])))
      param_names <- names(slots[[i]])
  }
  if (!length(keep)) return(NULL)
  M <- do.call(rbind, keep)
  if (is.null(param_names) || length(param_names) != k)
    param_names <- paste0("theta", seq_len(k))
  dimnames(M) <- list(labs, param_names)
  M
}

# ---- the same comparison, as ratios against one reference design ----------
# Row-wise M[i, ] / M[reference, ], for every row but the reference: how many
# times the competing design's error is that of the design the app computed.
# Above 1 the computed design has the smaller error for that parameter; below
# 1 the competitor does.  A ratio is the only honest way to read this table --
# the MSEs themselves live on whatever scale the parameter happens to have.
#
# `reference` is a design NAME, not a position, because .ui_sim_matrix() drops
# the designs that were not run: with only SRS in hand the first row would be
# SRS, and dividing it by itself would silently report 1.
#
# Returns a (designs - 1) x parameters matrix, or NULL when there is nothing
# to compare (the reference absent, or no other design run).
.ui_sim_ratio <- function(M, reference) {
  if (is.null(M) || !is.matrix(M) || nrow(M) < 2L) return(NULL)
  i <- match(reference, rownames(M))
  if (is.na(i)) return(NULL)
  ref <- M[i, ]
  ref[!is.finite(ref) | ref <= 0] <- NA_real_    # no ratio against a zero MSE
  R <- M[-i, , drop = FALSE]
  R <- sweep(R, 2, ref, "/")
  R[!is.finite(R)] <- NA_real_
  R
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
  # ONE efficiency, relative to the TRUE optimum: the ratio to the derived
  # reference design times that design's own certified bound (Becker & Yang,
  # Thm 4.5 / 4.6), so an unconverged reference cannot inflate it
  rb <- if (is.null(ref$efficiency_lower_bound) || !is.finite(ref$efficiency_lower_bound)) 1
        else ref$efficiency_lower_bound
  list(p = p_new, crit_design = crit_design, crit_ref = crit_ref,
       efficiency_lower_bound = if (is.finite(eff)) min(eff, 1) * rb else NA_real_,
       converged = isTRUE(ref$converged),
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
