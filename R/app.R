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
                           link = "identity", ncat = NULL) {
  if (length(covariates) == 0L)
    stop("add at least one covariate.", call. = FALSE)
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
       finest = db$finest, link = link, ncat = ncat,
       f = nz(f), x = nz(x), fx = nz(fx), xx = nz(xx), ff = nz(ff))
}

# Coefficient labels for the spec (order matches the theta the solvers expect).
# Uses model_info_matrix(), which carries "coef_names" for every link and needs
# no theta.
.ui_coef_names <- function(spec) {
  mf <- model_info_matrix(design_box = spec$design_box, link = spec$link,
                          f = spec$f, x = spec$x, fx = spec$fx, xx = spec$xx,
                          ff = spec$ff, ncat = spec$ncat)
  attr(mf, "coef_names")
}

#' Launch the owea point-and-click web app.
#'
#' Opens a Shiny application that builds D- and A-optimal designs from
#' point-and-click inputs -- no R code required. Choose a model family, add
#' covariates (continuous ranges or factor levels), pick interactions by name,
#' set the criterion and parameter values, and get the design, its plot, the
#' information matrix, and a downloadable table. This is a convenience wrapper
#' over \code{\link{optimal_design}} / \code{\link{exact_design}}.
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
