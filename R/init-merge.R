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

# 1-based initial support indices for the unified engine.  Returns integer(0)
# to let the C++ engine run IBOSS (vector mode default).
.initial_support_idx <- function(X, info_mode, info_data, k, infor0,
                                 init_method = "auto", max_tries = 50L) {
  if (identical(init_method, "auto"))
    init_method <- if (info_mode == 0L) "iboss" else "minmax"
  if (identical(init_method, "iboss")) {
    if (info_mode != 0L)
      stop("init_method = \"iboss\" is only available for information-vector ",
           "input; use \"minmax\" for information-matrix input.", call. = FALSE)
    return(integer(0))                          # C++ runs IBOSS
  }
  if (!init_method %in% c("minmax", "minmaxmedian", "random"))
    stop("init_method must be one of \"auto\", \"iboss\", \"minmax\", ",
         "\"minmaxmedian\", \"random\".", call. = FALSE)

  n <- nrow(X)
  ok <- function(idx) {
    idx <- unique(idx)
    if (length(idx) < 1L) return(NULL)
    w  <- rep(1 / length(idx), length(idx))
    Ic <- infor0 + .info_ind_R(idx, w, info_mode, info_data, k)
    if (det(Ic) > 1e-20 && rcond(Ic) > 1e-14) idx else NULL
  }

  if (init_method != "random") {
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
#'   \code{"random"}, \code{"iboss"}, \code{"auto"}.
#' @param max_tries number of random subsets to try in the fallback.
#' @return list with elements \code{support} and \code{weights}.
#' @export
initial_support <- function(prob, init_method = "minmax", max_tries = 50L) {
  idx <- .initial_support_idx(prob$X, prob$info_mode, prob$info_data, prob$k,
                              prob$infor0, init_method, max_tries)
  if (length(idx) == 0L)
    idx <- .initial_support_idx(prob$X, prob$info_mode, prob$info_data, prob$k,
                                prob$infor0, "minmax", max_tries)
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
