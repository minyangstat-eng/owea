# ===========================================================================
# init-merge.R -- initial-support selection and optional neighbourhood merging.
# ===========================================================================

# Per-point k x k info from a (possibly b-scaled) storage column.
.point_info_R <- function(col, info_mode, k) {
  if (info_mode == 0L) tcrossprod(col)
  else { M <- matrix(col, k, k); 0.5 * (M + t(M)) }
}

# Information matrix of an index set with given weights, from scaled info_data
# (candidate contribution only; add infor0 outside for the combined matrix).
.info_ind_R <- function(idx, w, info_mode, info_data, k) {
  M <- matrix(0.0, k, k)
  for (i in seq_along(idx))
    M <- M + w[i] * .point_info_R(info_data[, idx[i]], info_mode, k)
  M
}

# Deterministic "minmax" start: for each covariate in turn, take the
# not-yet-chosen points with the smallest and largest value (and, if
# with_median, the median).  Returns row indices into X.
.deterministic_support_idx <- function(X, with_median) {
  n <- nrow(X); N <- ncol(X)
  selected <- integer(0)
  for (d in seq_len(N)) {
    rem <- setdiff(seq_len(n), selected)
    if (length(rem) == 0L) break
    rs <- rem[order(X[rem, d])]
    m  <- length(rs)
    picks <- if (with_median) c(rs[1L], rs[(m + 1L) %/% 2L], rs[m])
             else             c(rs[1L], rs[m])
    selected <- union(selected, picks)
  }
  selected
}

# General multiplicative algorithm for optimal design weights (Yu 2010,
# Ann. Statist. 38(3)).  Iterates w_i <- w_i * d_i^lambda / sum(w_j d_j^lambda),
# where d_i = tr(phi'(M) A_i) is the criterion sensitivity, and returns the weight
# vector over all N candidate columns.  Covers, via `wb` (the v x k quantity-of-
# interest matrix) and `pp`:
#   pp == 0 (D / D_A):  phi'(M) = M^-1 wb' (wb M^-1 wb')^-1 wb M^-1 ,  lambda = 1
#   pp == 1 (A / A_A):  phi'(M) = M^-1 wb' wb M^-1 ,                   lambda = 1/2
# and an existing-design offset via `infor0` (M = infor0 + sum_i w_i A_i); the
# stopping bound is the running mu = sum_j w_j d_j (= v for full-D no-offset).
# `wb = NULL` means full-parameter D-optimality (wb = I), which keeps a fast path
# used by the D-optimal warm start (init_method = "MA") and D-optimal solver.
# `info_data` and `infor0` are the already b-scaled inputs from the callers.
.multiplicative_weights <- function(info_data, info_mode, k, infor0,
                                    wb = NULL, pp = 0L,
                                    max_iter = 100L, tol = 1e-8,
                                    details = FALSE) {
  N <- ncol(info_data)
  full_D <- is.null(wb) && pp == 0L             # wb = I, D-optimality: fast path
  if (is.null(wb)) wb <- diag(k)
  wb <- as.matrix(wb)
  lambda <- if (pp == 0L) 1.0 else 0.5
  w <- rep(1 / N, N)
  d <- rep(NA_real_, N); iter <- 0L; singular <- FALSE
  for (it in seq_len(max_iter)) {
    iter <- it
    # M = infor0 + sum_i w_i A_i   (info_data columns are already b-scaled)
    if (info_mode == 0L) M <- infor0 + info_data %*% (t(info_data) * w)
    else { M <- infor0 + matrix(info_data %*% w, k, k); M <- 0.5 * (M + t(M)) }
    Mi <- tryCatch(solve(M), error = function(e) NULL)
    if (is.null(Mi)) { singular <- TRUE; break }
    # criterion gradient C = phi'(M); sensitivity d_i = tr(C A_i)
    if (full_D) {
      C <- Mi
    } else {
      WMi <- wb %*% Mi                          # v x k = wb M^-1
      if (pp == 0L) {
        Si <- tryCatch(solve(WMi %*% t(wb)), error = function(e) NULL)  # (wb M^-1 wb')^-1
        if (is.null(Si)) { singular <- TRUE; break }
        C <- crossprod(WMi, Si %*% WMi)         # k x k
      } else {
        C <- crossprod(WMi)                     # k x k = M^-1 wb' wb M^-1
      }
    }
    if (info_mode == 0L) d <- colSums(info_data * (C %*% info_data))  # f_i' C f_i
    else                 d <- as.numeric(as.numeric(C) %*% info_data)  # tr(C A_i)
    d[d < 0] <- 0
    mu <- sum(w * d)                            # equivalence-theorem bound
    if (!is.finite(mu) || mu <= 0) { singular <- TRUE; break }
    # Stop on the ABSOLUTE equivalence-theorem gap: verify_equiv reports
    # max_d = coeff * (max(d) - mu), so max(d) - mu <= tol guarantees the
    # certified max_d <= tol (coeff <= 1).  A relative rule would overshoot by a
    # factor mu (= v), leaving the design just short of certification.
    if (max(d) - mu <= tol) {
      w <- w * d^lambda / sum(w * d^lambda)     # one more clean-up step
      break
    }
    dl <- d^lambda
    sdl <- sum(w * dl)
    if (!is.finite(sdl) || sdl <= 0) { singular <- TRUE; break }
    w <- w * dl / sdl
  }
  if (details) list(w = w, d = d, iter = iter, singular = singular) else w
}

# Direct multiplicative-algorithm solver on a fixed candidate grid.  Runs the
# general multiplicative algorithm (see .multiplicative_weights) to convergence
# and returns its design -- the stopping rule max_i d_i <= mu IS the equivalence-
# theorem certificate, so no OWEA exchange engine is needed.  Returns the same
# list shape as optimal_design()'s solve_prepared() (support / weights /
# criterion / max_d / converged / iterations, plus a `singular` flag so the caller
# can fall back to OWEA).  Handles D- and A-optimality (`pp`), any quantity of
# interest (`wb`), and an existing design (`infor0`), all via the already b-scaled
# `scaled` information.  The reported criterion / max_d come from the same
# criterion_cpp / verify_equiv_cpp used by the OWEA path.
.ma_solve <- function(X, scaled, info_mode, k, wb, pp, infor0, eps0, max_iter) {
  mw <- .multiplicative_weights(scaled, info_mode, k, infor0, wb = wb, pp = pp,
                                max_iter = max_iter, tol = eps0, details = TRUE)
  w <- mw$w
  if (mw$singular)                              # let the caller fall back to OWEA
    return(list(support = X[0, , drop = FALSE], weights = numeric(0),
                criterion = Inf, max_d = Inf, converged = FALSE,
                iterations = mw$iter, singular = TRUE))
  # Drop numerically-negligible weights and renormalise, but only as far as the
  # design stays certified: for sharply-peaked weight distributions (e.g. subset
  # criteria) an aggressive threshold drops meaningful support and breaks
  # optimality.  Try increasingly-complete supports until the equivalence check
  # passes (or use all positive weights).  No weight re-optimisation is needed --
  # the multiplicative weights are already (near-)optimal.
  wmax <- max(w)
  keep <- which(w > 0); w_sel <- w[keep] / sum(w[keep]); ve <- NULL
  for (thr in c(1e-6, 1e-9, 1e-12, 0)) {
    keep  <- which(w > thr * wmax)
    w_sel <- w[keep] / sum(w[keep])
    opt_infor <- infor0 + .info_ind_R(keep, w_sel, info_mode, scaled, k)
    ve <- verify_equiv_cpp(as.integer(pp), as.matrix(wb), as.integer(info_mode),
                           scaled, opt_infor, as.matrix(infor0))
    if (ve$max_d <= eps0) break
  }
  crit <- criterion_cpp(as.integer(pp), as.integer(keep), as.numeric(w_sel),
                        as.integer(info_mode), scaled, as.matrix(wb),
                        as.matrix(infor0))
  list(support = X[keep, , drop = FALSE], weights = w_sel,
       criterion = crit, max_d = ve$max_d, converged = ve$max_d <= eps0,
       iterations = mw$iter, singular = FALSE)
}

# 1-based initial support indices for the unified engine.  Returns integer(0)
# to let the C++ engine run IBOSS (vector mode default).
.initial_support_idx <- function(X, info_mode, info_data, k, infor0,
                                 init_method = "auto", max_tries = 50L,
                                 ma_max_iter = 100L) {
  if (identical(tolower(init_method), "ma")) init_method <- "MA"
  if (identical(init_method, "auto"))
    init_method <- if (info_mode == 0L) "iboss" else "minmax"
  if (identical(init_method, "iboss")) {
    if (info_mode != 0L)
      stop("init_method = \"iboss\" is only available for information-vector ",
           "input; use \"minmax\" for information-matrix input.", call. = FALSE)
    return(integer(0))                          # C++ runs IBOSS
  }
  if (!init_method %in% c("minmax", "minmaxmedian", "random", "MA"))
    stop("init_method must be one of \"auto\", \"iboss\", \"minmax\", ",
         "\"minmaxmedian\", \"random\", \"MA\".", call. = FALSE)

  n <- nrow(X)
  ok <- function(idx) {
    idx <- unique(idx)
    if (length(idx) < 1L) return(NULL)
    w  <- rep(1 / length(idx), length(idx))
    Ic <- infor0 + .info_ind_R(idx, w, info_mode, info_data, k)
    if (det(Ic) > 1e-20 && rcond(Ic) > 1e-14) idx else NULL
  }

  if (init_method == "MA") {
    w   <- .multiplicative_weights(info_data, info_mode, k, infor0,
                                   max_iter = ma_max_iter)
    ord <- order(w, decreasing = TRUE)
    npk <- min(k + 1L, n)
    idx <- ok(ord[seq_len(npk)])
    if (!is.null(idx)) return(idx)
    # The top-(k+1) set can be rank-deficient (e.g. a pure-corner design for a
    # model with quadratic terms, where every x_i^2 collapses onto the intercept).
    # Rather than fall back to the fragile random search, greedily extend down the
    # MA weight order -- accumulating the (equal-weight) information incrementally
    # -- until the combined matrix is non-singular.  This is deterministic and
    # succeeds whenever the full candidate set has full rank.
    if (npk < n) {
      S <- matrix(0.0, k, k)
      for (j in seq_len(npk)) S <- S + .point_info_R(info_data[, ord[j]], info_mode, k)
      for (j in seq.int(npk + 1L, n)) {
        S  <- S + .point_info_R(info_data[, ord[j]], info_mode, k)
        Ic <- infor0 + S / j
        if (det(Ic) > 1e-20 && rcond(Ic) > 1e-14) return(ord[seq_len(j)])
      }
    }
    # extremely rare (full set singular): fall through to the shared fallbacks.
  }

  if (!init_method %in% c("random", "MA")) {
    idx <- ok(.deterministic_support_idx(X, with_median = (init_method == "minmaxmedian")))
    if (!is.null(idx)) return(idx)
  }
  n_init <- min(k + 1L, n)
  idx <- ok(round(seq(1, n, length.out = n_init)))
  if (!is.null(idx)) return(idx)
  for (t in seq_len(max_tries)) {
    idx <- ok(sort(sample.int(n, n_init, replace = TRUE)))
    if (!is.null(idx)) return(idx)
  }
  stop("Could not find a non-singular starting design after ", max_tries,
       " tries.", call. = FALSE)
}

#' Build an initial support for a design problem (information-matrix model).
#'
#' Convenience wrapper kept for compatibility: returns a starting design whose
#' combined information matrix is non-singular.
#'
#' @param prob a \code{DesignProblem}.
#' @param init_method one of \code{"minmax"} (default), \code{"minmaxmedian"},
#'   \code{"random"}, \code{"iboss"}, \code{"MA"}, \code{"auto"}. \code{"MA"} runs
#'   a multiplicative algorithm (equal-weight start, D-optimality) and keeps the
#'   \code{k + 1} highest-weight points as the initial support.
#' @param max_tries number of random subsets to try in the fallback.
#' @param ma_max_iter maximum number of multiplicative-algorithm iterations when
#'   \code{init_method = "MA"} (default 100); ignored by the other methods.
#' @return list with elements \code{support} and \code{weights}.
#' @export
initial_support <- function(prob, init_method = "minmax", max_tries = 50L,
                            ma_max_iter = 100L) {
  idx <- .initial_support_idx(prob$X, prob$info_mode, prob$info_data, prob$k,
                              prob$infor0, init_method, max_tries, ma_max_iter)
  if (length(idx) == 0L)
    idx <- .initial_support_idx(prob$X, prob$info_mode, prob$info_data, prob$k,
                                prob$infor0, "minmax", max_tries, ma_max_iter)
  list(support = prob$X[idx, , drop = FALSE],
       weights = rep(1 / length(idx), length(idx)))
}

#' Merge support points whose Euclidean distance is below \code{atol}.
#'
#' Iterative weighted-centroid merge: any two points closer than \code{atol}
#' are replaced by their weight-weighted average, with the weights pooled;
#' repeats until no further merge is possible.
#'
#' @param support \eqn{m \times N} matrix of support points (one per row).
#' @param weights numeric vector of length m.
#' @param atol distance threshold below which two points are merged.
#' @return list with the merged \code{support} and \code{weights}.
#' @export
merge_close_points <- function(support, weights, atol = 1e-2) {
  support <- as.matrix(support)
  changed <- TRUE
  while (changed && nrow(support) > 1L) {
    changed <- FALSE
    n <- nrow(support)
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        d <- sqrt(sum((support[i, ] - support[j, ])^2))
        if (d < atol) {
          w_new <- weights[i] + weights[j]
          support[i, ] <- (weights[i] * support[i, ] +
                           weights[j] * support[j, ]) / w_new
          weights[i] <- w_new
          support <- support[-j, , drop = FALSE]
          weights <- weights[-j]
          changed <- TRUE
          break
        }
      }
      if (changed) break
    }
  }
  list(support = support, weights = weights)
}
