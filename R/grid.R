# ===========================================================================
# grid.R -- candidate-grid construction and neighbourhood refinement.
# ===========================================================================

#' Build a full rectangular candidate grid.
#'
#' Reconciles the two calling conventions of the original packages: it accepts
#' \code{lower}/\code{upper} corners (information-vector version) and uses the
#' robust point count of the information-matrix version.
#'
#' @param lower,upper numeric vectors: the range of each covariate.
#' @param by grid step: a scalar (same step for every covariate) or a vector
#'   with one step per covariate.
#' @return an \eqn{N \times d} matrix, one row per grid point, the first
#'   dimension varying fastest.
#' @export
make_grid <- function(lower, upper, by) {
  d <- length(lower)
  if (length(upper) != d)
    stop("'lower' and 'upper' must have the same length", call. = FALSE)
  if (length(by) == 1L) by <- rep(by, d)
  if (length(by) != d)
    stop("'by' must have length 1 or length(lower)", call. = FALSE)
  axes <- lapply(seq_len(d), function(j) {
    npt <- round((upper[j] - lower[j]) / by[j]) + 1L
    lower[j] + (seq_len(npt) - 1L) * by[j]
  })
  X <- as.matrix(expand.grid(axes))
  colnames(X) <- if (!is.null(names(lower))) names(lower)
                 else paste0("x", seq_len(d))
  rownames(X) <- NULL
  X
}

#' Discretize a design box into a fixed candidate set.
#'
#' Handles both CONTINUOUS and FACTOR (categorical) covariates, using the same
#' \code{design_box} convention as \code{\link{optimal_design}}: an entry
#' \code{c(lo, hi)} is a continuous covariate discretized by \code{step}, while a
#' single integer \code{L >= 2} is a factor whose levels \code{1..L} are always
#' enumerated in full (factors are not affected by \code{step}).
#'
#' @param design_box list with one entry per covariate: a \code{c(lo, hi)} pair
#'   (continuous) or a single integer \code{L >= 2} (a factor with levels
#'   \code{1..L}). For example \code{list(c(2), c(3), c(0, 6))} is two factors
#'   (2 and 3 levels) and one continuous covariate on \code{[0, 6]}.
#' @param step grid spacing for the CONTINUOUS covariates: a scalar (same
#'   spacing for every continuous covariate), a vector with one entry per
#'   continuous covariate, or a vector with one entry per covariate (factor
#'   entries are ignored). Factor covariates always enumerate all their levels.
#' @return an \eqn{n \times N} numeric matrix of all grid points; factor columns
#'   hold integer levels \code{1..L}.
#' @export
candidate_grid <- function(design_box, step) {
  meta      <- .parse_design_box(design_box)
  is_factor <- meta$is_factor
  by        <- .expand_stage_step(as.numeric(step), is_factor,
                                  sum(!is_factor), length(design_box))
  unname(.factor_make_grid(meta$lo, meta$hi, by, is_factor, meta$nlevels))
}

# Expand one stage's step specification into a full length-N per-covariate `by`
# vector.  `s` may be a scalar (same step for every continuous covariate), a
# vector of length n_continuous, or a vector of length N (one per covariate,
# factor entries ignored).  Factor dims get a placeholder step of 1 (their axes
# always enumerate all levels downstream).
.expand_stage_step <- function(s, is_factor, ncont, N) {
  s  <- as.numeric(s)
  by <- numeric(N)
  if (length(s) == 1L) {
    by[!is_factor] <- s
  } else if (length(s) == ncont) {
    by[!is_factor] <- s
  } else if (length(s) == N) {
    by[!is_factor] <- s[!is_factor]
  } else {
    stop(sprintf(paste0("each grid step must be a scalar, or have length %d ",
                        "(continuous covariates) or %d (all covariates)."),
                 ncont, N), call. = FALSE)
  }
  if (any(!is_factor & (!is.finite(by) | by <= 0)))
    stop("continuous grid steps must be finite and positive.", call. = FALSE)
  by[is_factor] <- 1
  by
}

# Normalize a user step_sequence into a list of per-covariate `by` vectors, one
# per refinement stage (coarsest first).  `step_sequence` may be:
#   * a numeric vector -- each element is one stage's scalar step, applied to
#     every continuous covariate (the classic uniform-step form);
#   * a list           -- each element is one stage's per-covariate steps, a
#     numeric vector of length 1, n_continuous, or N (see .expand_stage_step);
#   * a matrix         -- one row per stage, columns as for the list form.
# Returns a list of length-N numeric `by` vectors (factor entries set to 1).
.normalize_step_sequence <- function(step_sequence, is_factor) {
  N     <- length(is_factor)
  ncont <- sum(!is_factor)
  stages <-
    if (is.list(step_sequence))        step_sequence
    else if (is.matrix(step_sequence)) lapply(seq_len(nrow(step_sequence)),
                                              function(i) step_sequence[i, ])
    else                               as.list(as.numeric(step_sequence))
  lapply(stages, function(s) .expand_stage_step(s, is_factor, ncont, N))
}

# Internal: factor-aware full grid.  Continuous dims are stepped by `by`; factor
# dims (is_factor[j] TRUE) always enumerate their levels 1..nlevels[j].  When
# there are no factors this is exactly make_grid(lower, upper, by).
.factor_make_grid <- function(lower, upper, by, is_factor = NULL, nlevels = NULL) {
  d <- length(lower)
  if (is.null(is_factor) || !any(is_factor)) return(make_grid(lower, upper, by))
  if (length(by) == 1L) by <- rep(by, d)
  lo <- lower; hi <- upper; byv <- by
  for (j in which(is_factor)) { lo[j] <- 1; hi[j] <- nlevels[j]; byv[j] <- 1 }
  make_grid(lo, hi, byv)        # factor dims yield exactly 1..nlevels[j]
}

# Internal: factor-aware refined grid.  Continuous dims use the same
# neighbourhood-refinement logic as refined_grid(); factor dims always enumerate
# ALL their levels at every stage (categorical dims are not refined).  Delegates
# to refined_grid() when there are no factors.
.factor_refined_grid <- function(box_lo, box_hi, prev_support, neighbor_radius,
                                 new_step, is_factor = NULL, nlevels = NULL) {
  if (is.null(is_factor) || !any(is_factor))
    return(refined_grid(box_lo, box_hi, prev_support, neighbor_radius, new_step))
  N    <- length(box_lo)
  if (length(new_step) == 1L)        new_step        <- rep(new_step, N)
  if (length(neighbor_radius) == 1L) neighbor_radius <- rep(neighbor_radius, N)
  Nmax <- round((box_hi - box_lo) / new_step)
  prev_support <- as.matrix(prev_support)
  pieces <- list()
  for (s in seq_len(nrow(prev_support))) {
    x <- prev_support[s, ]
    axes <- lapply(seq_len(N), function(d) {
      if (is_factor[d]) return(seq_len(nlevels[d]))      # all levels, every stage
      jc  <- (x[d] - box_lo[d]) / new_step[d]
      jlo <- max(0,       ceiling(jc - neighbor_radius[d] / new_step[d]))
      jhi <- min(Nmax[d], floor(  jc + neighbor_radius[d] / new_step[d]))
      if (jlo <= jhi) box_lo[d] + (jlo:jhi) * new_step[d] else numeric(0)
    })
    if (all(lengths(axes) > 0))
      pieces[[length(pieces) + 1L]] <- unname(as.matrix(expand.grid(axes)))
    snap <- numeric(N)
    for (d in seq_len(N))
      snap[d] <- if (is_factor[d]) min(max(round(x[d]), 1L), nlevels[d])
                 else box_lo[d] + round((x[d] - box_lo[d]) / new_step[d]) * new_step[d]
    pieces[[length(pieces) + 1L]] <- matrix(snap, nrow = 1L)
  }
  unique(do.call(rbind, pieces))
}

# Internal: refined grid for one stage of appro_opt_seq() (information-vector
# convention).  Keeps the by_cur-spaced lattice points lying within
# +/- 2 * by_prev of any current support point, clipped to [lower, upper].
.refine_grid <- function(centers, by_prev, by_cur, lower, upper) {
  d    <- ncol(centers)
  Mmax <- floor((upper - lower) / by_cur + 1e-9)
  idx  <- vector("list", nrow(centers))
  for (i in seq_len(nrow(centers))) {
    cc <- as.numeric(centers[i, ])
    axes <- lapply(seq_len(d), function(j) {
      m_lo <- ceiling((cc[j] - 2 * by_prev[j] - lower[j]) / by_cur[j] - 1e-9)
      m_hi <- floor  ((cc[j] + 2 * by_prev[j] - lower[j]) / by_cur[j] + 1e-9)
      m_lo <- max(0, m_lo)
      m_hi <- min(Mmax[j], m_hi)
      if (m_hi < m_lo) {
        mm <- round((cc[j] - lower[j]) / by_cur[j])
        min(max(mm, 0), Mmax[j])
      } else m_lo:m_hi
    })
    idx[[i]] <- as.matrix(expand.grid(axes))
  }
  idx <- unique(do.call(rbind, idx))
  X <- sweep(sweep(idx, 2, by_cur, `*`), 2, lower, `+`)
  colnames(X) <- colnames(centers)
  rownames(X) <- NULL
  X
}

# Internal: refined grid for one stage of optimal_design() (information-matrix
# convention).  All points on the global step-`new_step` grid within an
# infinity-radius `neighbor_radius` of some point in `prev_support`, clipped to
# the box, plus the snapped previous support so a warm start is representable.
refined_grid <- function(box_lo, box_hi, prev_support, neighbor_radius, new_step) {
  N    <- length(box_lo)
  if (length(new_step) == 1L)        new_step        <- rep(new_step, N)
  if (length(neighbor_radius) == 1L) neighbor_radius <- rep(neighbor_radius, N)
  Nmax <- round((box_hi - box_lo) / new_step)
  pieces <- list()
  prev_support <- as.matrix(prev_support)
  for (s in seq_len(nrow(prev_support))) {
    x <- prev_support[s, ]
    axes <- lapply(seq_len(N), function(d) {
      jc  <- (x[d] - box_lo[d]) / new_step[d]
      jlo <- max(0,       ceiling(jc - neighbor_radius[d] / new_step[d]))
      jhi <- min(Nmax[d], floor(  jc + neighbor_radius[d] / new_step[d]))
      if (jlo <= jhi) box_lo[d] + (jlo:jhi) * new_step[d] else numeric(0)
    })
    if (all(lengths(axes) > 0))
      pieces[[length(pieces) + 1L]] <- unname(as.matrix(expand.grid(axes)))
    snap <- box_lo + round((x - box_lo) / new_step) * new_step
    pieces[[length(pieces) + 1L]] <- matrix(snap, nrow = 1L)
  }
  unique(do.call(rbind, pieces))
}
