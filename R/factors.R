# ===========================================================================
# factors.R -- support for FACTOR (categorical) covariates.
#
# A covariate is either CONTINUOUS (a c(lo, hi) range, stepped by a grid) or a
# FACTOR with L levels (the integers 1..L).  Factor covariates are carried,
# evaluated and reported as their integer LEVEL: the model functions
# info_vector(x, theta) / info_matrix(x) receive the raw level for a factor
# dimension, and any contrast coding is the user's responsibility inside their
# own model function.  The design space enumerates the levels 1..L.
# ===========================================================================

# Parse a design_box into factor metadata.  Element c(L) (length 1) is a factor
# with L levels; element c(lo, hi) (length 2) is a continuous covariate.
.parse_design_box <- function(design_box) {
  d <- length(design_box)
  is_factor <- logical(d)
  nlevels   <- rep(NA_integer_, d)
  lo <- numeric(d); hi <- numeric(d)
  for (j in seq_len(d)) {
    b <- as.numeric(design_box[[j]])
    if (length(b) == 1L) {
      L <- b
      if (!is.finite(L) || L < 2 || abs(L - round(L)) > 1e-8)
        stop(sprintf(paste0("design_box[[%d]]: a factor must be a single ",
                            "integer >= 2 (the number of levels); got %s."),
                     j, format(L)), call. = FALSE)
      is_factor[j] <- TRUE
      nlevels[j]   <- as.integer(round(L))
      lo[j] <- 1; hi[j] <- as.integer(round(L))
    } else if (length(b) == 2L) {
      if (!all(is.finite(b)) || b[2] < b[1])
        stop(sprintf(paste0("design_box[[%d]]: a continuous covariate must be ",
                            "c(lo, hi) with hi >= lo."), j), call. = FALSE)
      lo[j] <- b[1]; hi[j] <- b[2]
    } else {
      stop(sprintf(paste0("design_box[[%d]] must have length 1 (a factor: the ",
                          "number of levels) or 2 (a continuous covariate: ",
                          "c(lo, hi))."), j), call. = FALSE)
    }
  }
  list(is_factor = is_factor, nlevels = nlevels, lo = lo, hi = hi)
}

# Normalize the candidate_set-path 'factor_levels' argument (one entry per
# candidate_set column) into the same is_factor / nlevels shape the parser
# returns.  NA / 0 / 1 => continuous; an integer >= 2 => a factor with L levels.
.factor_levels_to_meta <- function(factor_levels, ncol) {
  is_factor <- logical(ncol)
  nlevels   <- rep(NA_integer_, ncol)
  if (is.null(factor_levels))
    return(list(is_factor = is_factor, nlevels = nlevels))
  fl <- as.numeric(factor_levels)
  if (length(fl) != ncol)
    stop(sprintf(paste0("'factor_levels' must have one entry per candidate_set ",
                        "column (%d)."), ncol), call. = FALSE)
  for (j in seq_len(ncol)) {
    L <- fl[j]
    if (is.na(L) || L == 0 || L == 1) next                     # continuous
    if (!is.finite(L) || L < 2 || abs(L - round(L)) > 1e-8)
      stop(sprintf(paste0("factor_levels[%d] must be NA/0/1 (continuous) or an ",
                          "integer >= 2 (factor levels); got %s."),
                   j, format(L)), call. = FALSE)
    is_factor[j] <- TRUE
    nlevels[j]   <- as.integer(round(L))
  }
  list(is_factor = is_factor, nlevels = nlevels)
}

# Format one support-point row for printing: factor columns as integer levels
# (no decimals), continuous columns with `cont_fmt`.  is_factor = NULL (or all
# FALSE) reproduces the all-continuous formatting.
.fmt_support_row <- function(row, is_factor, cont_fmt = "%.6f") {
  if (is.null(is_factor)) is_factor <- logical(length(row))
  vapply(seq_along(row), function(j)
    if (isTRUE(is_factor[j])) sprintf("%d", as.integer(round(row[j])))
    else sprintf(cont_fmt, row[j]), character(1))
}

# A valid probe point for inferring k: continuous dims at the box centre, factor
# dims at level 1 (the box centre of a factor, e.g. 1.5, is not a valid level).
.factor_probe_point <- function(lo, hi, is_factor) {
  p <- (lo + hi) / 2
  p[is_factor] <- 1
  p
}

# Round + clamp factor columns of a matrix to integer levels 1..L (continuous
# columns untouched).  Used for representative-point sampling.
.snap_factor_levels <- function(M, is_factor, nlevels) {
  M <- as.matrix(M)
  for (j in which(is_factor)) {
    v <- round(M[, j])
    v[v < 1L]          <- 1L
    v[v > nlevels[j]]  <- nlevels[j]
    M[, j] <- v
  }
  M
}

# Hard check that data supplied for factor covariates (candidate_set rows,
# xi0_points) carry integer levels in 1..L.
.validate_factor_columns <- function(M, is_factor, nlevels, what = "candidate_set") {
  if (!any(is_factor)) return(invisible(TRUE))
  M <- as.matrix(M)
  if (ncol(M) != length(is_factor))
    stop(sprintf("%s has %d columns but the factor specification has %d.",
                 what, ncol(M), length(is_factor)), call. = FALSE)
  for (j in which(is_factor)) {
    v <- M[, j]
    if (any(!is.finite(v)) || any(abs(v - round(v)) > 1e-8) ||
        any(v < 1) || any(v > nlevels[j]))
      stop(sprintf(paste0("%s column %d is a factor with %d levels; entries ",
                          "must be integers in 1..%d."),
                   what, j, nlevels[j], nlevels[j]), call. = FALSE)
  }
  invisible(TRUE)
}
