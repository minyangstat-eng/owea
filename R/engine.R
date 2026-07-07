# ===========================================================================
# engine.R -- internal bridge to the unified C++ solver, plus the optional
# neighbourhood-merge post-processing step.
# ===========================================================================

# Run the unified C++ solver on a prepared, b-scaled candidate set.
# Returns the raw C++ list: index, weight, sensitivity, iter, value.
# min_support defaults to the rank-aware floor (k for single-stage, smaller for
# multistage) so an existing design does not force k spurious new points.
.solve_engine <- function(pp, wb, info_mode, info_data, infor0,
                          init_idx = integer(0),
                          max_iter = 100L, tol = 1e-6, verbose = FALSE,
                          min_support = NULL) {
  infor0 <- as.matrix(infor0)
  k <- if (info_mode == 0L) nrow(as.matrix(info_data))
       else as.integer(round(sqrt(nrow(as.matrix(info_data)))))
  # fallback assumes rank-1 per-point information; callers that know the true
  # per-point rank pass min_support explicitly (see .min_support_rule()).
  if (is.null(min_support))
    min_support <- .min_support_rule(infor0, k, as.matrix(wb), 1L)
  appro_opt_cpp(as.integer(pp), as.matrix(wb), as.integer(info_mode),
                as.matrix(info_data), infor0,
                as.integer(init_idx), as.integer(min_support),
                as.integer(max_iter), as.numeric(tol), isTRUE(verbose))
}

# Candidate-design contribution b * I_xi for an explicit (possibly off-grid)
# support, evaluated through the user's model and scaled by b.
.opt_infor_from_support <- function(support, weights, b,
                                    info_mode, info_vector, info_matrix, theta, k) {
  support <- as.matrix(support)
  M <- matrix(0.0, k, k)
  for (i in seq_len(nrow(support)))
    M <- M + weights[i] *
      .info_mat_at(support[i, ], info_mode, info_vector, info_matrix, theta)
  b * M
}

# Optional neighbourhood merge (default off in the callers).  Merges support
# points within `atol`, recomputes their (off-grid) information through the
# model, and re-optimises the weights on the merged support.  Returns updated
# support / weights / value, or the inputs unchanged if nothing merged.
.apply_merge <- function(support, weights, atol, pp, wb,
                         info_mode, info_vector, info_matrix, theta,
                         infor0, b, k, min_support = NULL) {
  msup <- if (is.null(min_support))
            .min_support_rule(infor0, k, as.matrix(wb), 1L) else min_support
  support <- as.matrix(support)
  if (nrow(support) <= 1L) return(list(support = support, weights = weights,
                                       merged = FALSE))
  mc <- merge_close_points(support, weights, atol)
  if (nrow(mc$support) == nrow(support))
    return(list(support = support, weights = weights, merged = FALSE))

  ms <- as.matrix(mc$support); mw <- mc$weights
  m  <- nrow(ms)
  len <- if (info_mode == 0L) k else k * k
  info_data <- matrix(0.0, len, m)
  for (i in seq_len(m))
    info_data[, i] <- .info_col_at(ms[i, ], info_mode, info_vector,
                                   info_matrix, theta)
  info_data <- .scale_info(info_data, info_mode, b)

  ow <- optimize_weights_cpp(as.integer(pp), as.matrix(wb),
                             as.integer(info_mode), info_data,
                             as.matrix(infor0), as.integer(msup))
  list(support = ms[ow$index, , drop = FALSE], weights = ow$weight,
       value = ow$value, merged = TRUE)
}

#' Directional derivative scan for a given design (equivalence-theorem check).
#'
#' @param prob a \code{DesignProblem}.
#' @param support \eqn{m \times N} matrix of support points.
#' @param weights weights of the design.
#' @return list with \code{x} (the maximiser), \code{d} (max directional
#'   derivative) and \code{index} (its 1-based row in \code{prob$X}).
#' @export
find_best_point <- function(prob, support, weights) {
  opt_infor <- .opt_infor_from_support(support, weights, prob$b,
                                       prob$info_mode, prob$info_vector,
                                       prob$info_matrix, prob$theta, prob$k)
  r <- verify_equiv_cpp(as.integer(prob$p), as.matrix(prob$wb),
                        as.integer(prob$info_mode), as.matrix(prob$info_data),
                        opt_infor, as.matrix(prob$infor0))
  list(x = prob$X[r$index, ], d = r$max_d, index = r$index)
}
