# ===========================================================================
# info.R -- information representations and conversion helpers.
#
# The model enters the algorithm through one per-point function, supplied as
# EITHER an information vector  info_vector(x, theta) -> length-k  (the per-point
# information matrix is f f', "mode 0") OR an information matrix
# info_matrix(x) -> k x k ("mode 1").
# ===========================================================================

# Build the v x k selection matrix wb = dg/dtheta from any of wb / subset /
# grad_g (at most one); default identity (all parameters).
.wb_from <- function(wb, subset, grad_g, theta, k) {
  n_given <- sum(!is.null(wb), !is.null(subset), !is.null(grad_g))
  if (n_given > 1L)
    stop("Supply at most one of 'wb' / 'subset' / 'grad_g'.", call. = FALSE)
  if (!is.null(wb)) {
    wb <- as.matrix(wb); storage.mode(wb) <- "double"
    if (ncol(wb) != k)
      stop("wb must have ncol = k = ", k, ".", call. = FALSE)
    return(wb)
  }
  if (!is.null(subset)) {
    subset <- as.integer(subset)
    if (length(subset) < 1L || any(subset < 1L) || any(subset > k) ||
        anyDuplicated(subset) != 0L)
      stop("'subset' must be distinct parameter indices in 1:", k, ".",
           call. = FALSE)
    S <- matrix(0.0, length(subset), k)
    S[cbind(seq_along(subset), subset)] <- 1.0
    return(S)
  }
  if (!is.null(grad_g)) {
    if (!is.function(grad_g))
      stop("'grad_g' must be a function theta -> v x k matrix.", call. = FALSE)
    return(as.matrix(grad_g(theta)))
  }
  diag(k)
}

# Normalize a user info_matrix to the internal one-argument form info_matrix(x).
# Accepts BOTH info_matrix(x) (theta captured in the closure) and
# info_matrix(x, theta) (theta supplied to optimal_design / DesignProblem); the
# latter is detected by its argument count and bound to the given theta.
.normalize_info_matrix <- function(fn, theta) {
  if (is.null(fn)) return(NULL)
  if (!is.function(fn))
    stop("'info_matrix' must be a function.", call. = FALSE)
  np <- tryCatch(length(formals(fn)), error = function(e) 1L)
  if (length(np) != 1L || is.na(np)) np <- 1L
  if (np >= 2L) {
    if (is.null(theta))
      stop("'theta' is required when 'info_matrix' is defined as ",
           "function(x, theta).", call. = FALSE)
    force(fn); th <- as.numeric(theta)
    function(x) fn(x, th)
  } else {
    fn
  }
}

# TRUE if a user info_vector is written as function(x, theta) (needs theta),
# FALSE if it is function(x) (theta captured in the closure).
.info_vec_needs_theta <- function(fn) {
  np <- tryCatch(length(formals(fn)), error = function(e) 2L)
  if (length(np) != 1L || is.na(np)) np <- 2L
  np >= 2L
}

# Normalize a user info_vector to the internal two-argument form
# info_vector(x, theta). A one-argument function(x) is wrapped so it ignores
# the (unused) theta; a function(x, theta) is used as-is.
.normalize_info_vector <- function(fn, theta) {
  if (is.null(fn)) return(NULL)
  if (!is.function(fn))
    stop("'info_vector' must be a function.", call. = FALSE)
  if (.info_vec_needs_theta(fn)) fn
  else { force(fn); function(x, theta) fn(x) }
}

# Per-point k x k information matrix at x (used for existing designs / merging).
.info_mat_at <- function(x, info_mode, info_vector, info_matrix, theta) {
  if (info_mode == 0L) {
    f <- as.numeric(info_vector(x, theta))
    tcrossprod(f)
  } else {
    as.matrix(info_matrix(x))
  }
}

# The storage column for one candidate point: length-k (mode 0) or length-k^2
# column-major vec of the k x k matrix (mode 1).
.info_col_at <- function(x, info_mode, info_vector, info_matrix, theta) {
  if (info_mode == 0L) as.numeric(info_vector(x, theta))
  else as.numeric(as.matrix(info_matrix(x)))
}

# Build the candidate information matrix (k x N mode 0, or k^2 x N mode 1) by
# evaluating the model at every row of X.  Returns list(info_data, k).
.build_info_data <- function(X, info_mode, info_vector, info_matrix, theta) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  n <- nrow(X)
  if (n == 0L) stop("X has no rows (no design points).", call. = FALSE)
  c1 <- .info_col_at(X[1, ], info_mode, info_vector, info_matrix, theta)
  len <- length(c1)
  out <- matrix(0.0, len, n)
  out[, 1] <- c1
  if (n > 1L)
    for (i in 2:n)
      out[, i] <- .info_col_at(X[i, ], info_mode, info_vector, info_matrix, theta)
  k <- if (info_mode == 0L) len else as.integer(round(sqrt(len)))
  list(info_data = out, k = k)
}

# Existing-design baseline a*I_xi0 and the new-stage fraction b.  When n0 == 0
# this is a zero matrix and b = 1 (single-stage).  xi0_weights are design
# weights (proportions summing to 1), NOT counts -- the existing sample size is
# n0.  They are normalized (with a warning) if they do not sum to 1.
.make_infor0 <- function(xi0_points, xi0_weights, n0, n1,
                         info_mode, info_vector, info_matrix, theta, k) {
  has_existing <- !is.null(xi0_points) && length(xi0_weights) > 0L
  if (n0 <= 0 || !has_existing) {
    if (has_existing && n0 <= 0)
      warning("an existing design (xi0_points/xi0_weights) was supplied but ",
              "n0 = 0, so the existing design is IGNORED. Set n0 to the ",
              "existing sample size to use it.", call. = FALSE)
    return(list(infor0 = matrix(0.0, k, k), b = 1.0))
  }
  a <- n0 / (n0 + n1)
  b <- n1 / (n0 + n1)
  xi0m <- as.matrix(xi0_points); storage.mode(xi0m) <- "double"
  if (nrow(xi0m) != length(xi0_weights))
    stop("length(xi0_weights) must equal nrow(xi0_points).", call. = FALSE)
  sw <- sum(xi0_weights)
  if (sw <= 0)
    stop("xi0_weights must sum to a positive value.", call. = FALSE)
  if (abs(sw - 1) > 1e-6) {
    warning(sprintf(paste0("xi0_weights sum to %g, not 1; normalizing them to ",
                          "proportions. They are design WEIGHTS, not counts -- ",
                          "use n0 for the existing sample size."), sw),
            call. = FALSE)
    xi0_weights <- xi0_weights / sw
  }
  I0 <- matrix(0.0, k, k)
  for (i in seq_len(nrow(xi0m)))
    I0 <- I0 + xi0_weights[i] *
      .info_mat_at(xi0m[i, ], info_mode, info_vector, info_matrix, theta)
  list(infor0 = a * I0, b = b)
}

# Numerical rank of a matrix (singular values above a relative tolerance).
.mat_rank <- function(A, tol = 1e-8) {
  A <- as.matrix(A)
  if (length(A) == 0L || max(abs(A)) == 0) return(0L)
  d <- svd(A, nu = 0, nv = 0)$d
  as.integer(sum(d > tol * max(d)))
}

# Rank r of the per-point information matrix.  For a vector model the per-point
# information is f f' (rank 1).  For a matrix model it is evaluated at a few
# representative points and the largest rank (the structural rank) is taken.
.point_info_rank <- function(info_mode, info_vector, info_matrix, theta, pts) {
  if (info_mode == 0L) return(1L)
  pts <- as.matrix(pts); n <- nrow(pts)
  idx <- unique(round(seq(1, n, length.out = min(7L, n))))
  r <- 1L
  for (i in idx) r <- max(r, .mat_rank(as.matrix(info_matrix(pts[i, ]))))
  as.integer(r)
}

# Minimum number of NEW support points the engine must keep so the COMBINED
# information matrix (a*I_xi0 + new design) is non-singular on the quantity of
# interest.  With k1 = rank(wb) parameters of interest, each new point adding
# rank r, and the existing design already covering up to min(rank(I_xi0), k1) of
# those directions:
#     min_support = max(1, ceil( (k1 - min(rank(infor0), k1)) / r )).
# Single-stage (infor0 = 0) gives ceil(k1 / r); a full-rank existing design
# gives 1.
.min_support_rule <- function(infor0, k, wb, r) {
  k1 <- max(1L, .mat_rank(wb))
  r0 <- .mat_rank(infor0)
  rem <- max(0L, k1 - min(r0, k1))
  max(1L, as.integer(ceiling(rem / max(1L, as.integer(r)))))
}

# Pre-scale candidate information by b: vectors by sqrt(b), matrices by b, so
# the assembled candidate information equals b * I_xi.
.scale_info <- function(info_data, info_mode, b) {
  if (b == 1.0) return(info_data)
  if (info_mode == 0L) info_data * sqrt(b) else info_data * b
}

#' Build the k x N information-vector matrix from a design matrix and a model.
#'
#' @param X design matrix: one row per point, columns are the covariates.
#' @param infor_vec model function \code{infor_vec(x, theta)} returning the
#'   length-k information vector at one point.
#' @param theta parameter vector.
#' @return a \eqn{k \times N} matrix whose column \eqn{n} is \eqn{f(x_n,\theta)}.
#' @export
build_infor_vec_all <- function(X, infor_vec, theta) {
  X <- as.matrix(X)
  N <- nrow(X)
  if (N == 0L) stop("X has no rows (no design points)")
  f1 <- as.numeric(infor_vec(X[1, ], theta))
  k  <- length(f1)
  out <- matrix(0.0, k, N)
  out[, 1] <- f1
  if (N > 1L) for (i in 2:N) out[, i] <- as.numeric(infor_vec(X[i, ], theta))
  out
}

#' Information matrix of a design (information-vector model).
#'
#' Returns \eqn{M = \sum_i w_i f(x_i,\theta) f(x_i,\theta)'}.
#'
#' @param X support points, one row per point.
#' @param weight weights (default equal weights \eqn{1/n}).
#' @param infor_vec information-vector model, either \code{function(x)} or
#'   \code{function(x, theta)}. Omit when a model spec (\code{link} + terms) is
#'   given instead.
#' @param theta parameter vector; required only when \code{infor_vec} is written
#'   as \code{function(x, theta)}. For a model spec that needs \code{theta}
#'   (logit/loglinear), a missing \code{theta} is drawn from \eqn{N(0,1)} with a
#'   warning.
#' @param link,f,x,fx,ff,xx,intercept,coding optional formula-style model spec
#'   (as in \code{\link{optimal_design}}); an alternative to \code{infor_vec}.
#' @param factor_levels optional integer vector marking factor columns of
#'   \code{X}; when omitted with a model spec, factor columns are auto-detected
#'   (columns whose values are all positive integers).
#' @return the \eqn{k \times k} information matrix.
#' @export
infor_matrix <- function(X, weight = NULL, infor_vec = NULL, theta = NULL,
                         link = NULL, f = NULL, x = NULL, fx = NULL, xx = NULL,
                         ff = NULL, intercept = TRUE, coding = "zero-sum",
                         factor_levels = NULL, ncat = NULL) {
  if (!is.null(link) && link %in% c("multinomial", "cumulative"))
    stop("the '", link, "' model has no information-vector form; use ",
         "design_information() to get its design information matrix.",
         call. = FALSE)
  spec <- .resolve_model_spec(link, f, x, fx, xx, intercept, coding,
                              design_box = NULL, candidate_set = X,
                              factor_levels = factor_levels,
                              info_vector = infor_vec, info_matrix = NULL,
                              theta = theta, ff = ff, ncat = ncat)
  infor_vec <- spec$info_vector; theta <- spec$theta
  if (is.null(infor_vec))
    stop("Supply 'infor_vec' or a model spec ('link' + terms).", call. = FALSE)
  if (.info_vec_needs_theta(infor_vec) && is.null(theta))
    stop("'theta' is required when 'infor_vec' is function(x, theta).",
         call. = FALSE)
  infor_vec <- .normalize_info_vector(infor_vec, theta)
  IVA <- build_infor_vec_all(X, infor_vec, theta)
  N   <- ncol(IVA)
  if (is.null(weight)) weight <- rep(1 / N, N)
  weight <- as.numeric(weight)
  if (length(weight) != N)
    stop("length(weight) must equal the number of design points (", N, ")")
  M <- IVA %*% (weight * t(IVA))
  0.5 * (M + t(M))
}

#' Information matrix of a design (information-matrix model).
#'
#' Returns \eqn{M(\xi) = \sum_i w_i\, \mathrm{info\_matrix}(x_i)} from the
#' support points and weights of a design.
#'
#' @param design \eqn{m \times N} matrix of support points (one per row); a
#'   length-m numeric vector is accepted for a 1-D design.
#' @param weights numeric vector of length m.
#' @param info_matrix information-matrix model, either \code{function(x)} or
#'   \code{function(x, theta)}, returning the \eqn{k \times k} matrix. Omit when
#'   a model spec (\code{link} + terms) is given instead.
#' @param theta parameter vector; required only when \code{info_matrix} is
#'   written as \code{function(x, theta)}. For a model spec that needs
#'   \code{theta} (all links but identity), a missing \code{theta} is drawn from
#'   \eqn{N(0,1)} with a warning.
#' @param link,f,x,fx,ff,xx,intercept,coding,ncat optional formula-style model
#'   spec (as in \code{\link{optimal_design}}, including \code{"multinomial"} /
#'   \code{"cumulative"} with \code{ncat}); an alternative to \code{info_matrix}.
#' @param factor_levels optional integer vector marking factor columns of
#'   \code{design}; when omitted with a model spec, factor columns are
#'   auto-detected (columns whose values are all positive integers).
#' @return the \eqn{k \times k} information matrix of the design.
#' @export
design_information <- function(design, weights, info_matrix = NULL, theta = NULL,
                              link = NULL, f = NULL, x = NULL, fx = NULL,
                              xx = NULL, ff = NULL, intercept = TRUE,
                              coding = "zero-sum", factor_levels = NULL,
                              ncat = NULL) {
  design <- as.matrix(design); storage.mode(design) <- "double"
  m <- nrow(design)
  weights <- as.numeric(weights)
  if (length(weights) != m)
    stop("length(weights) (", length(weights),
         ") must equal the number of support points (", m, ").", call. = FALSE)
  if (!is.null(link)) {
    spec <- .resolve_model_spec(link, f, x, fx, xx, intercept, coding,
                                design_box = NULL, candidate_set = design,
                                factor_levels = factor_levels,
                                info_vector = NULL, info_matrix = info_matrix,
                                theta = theta, ff = ff, ncat = ncat)
    theta <- spec$theta
    if (!is.null(spec$info_matrix)) {            # multinomial / cumulative
      info_matrix <- spec$info_matrix
    } else {
      iv <- spec$info_vector
      info_matrix <- if (.info_vec_needs_theta(iv))
                       function(x) tcrossprod(as.numeric(iv(x, theta)))
                     else function(x) tcrossprod(as.numeric(iv(x)))
    }
  }
  if (!is.function(info_matrix))
    stop("Supply 'info_matrix' (a function x -> k x k matrix) or a model spec ",
         "('link' + terms).", call. = FALSE)
  info_matrix <- .normalize_info_matrix(info_matrix, theta)   # accept (x) or (x, theta)
  M <- weights[1] * as.matrix(info_matrix(design[1, ]))
  if (m > 1L)
    for (i in 2:m) M <- M + weights[i] * as.matrix(info_matrix(design[i, ]))
  M
}
