# ===========================================================================
# print.R -- pretty-printing of designs.
# ===========================================================================

#' Pretty-print a design returned by \code{\link{owea}}.
#'
#' @param name a title for the printout.
#' @param res a result list from \code{owea()}.
#' @return \code{res}, invisibly.
#' @export
print_design <- function(name, res) {
  S <- as.matrix(res$support); storage.mode(S) <- "double"
  perm <- do.call(order, lapply(seq_len(ncol(S)), function(j) S[, j]))
  S <- S[perm, , drop = FALSE]; w <- res$weights[perm]
  cat(sprintf("\n=== %s ===\n", name))
  cat(sprintf("converged = %s    iter = %d    max d = %.3e    crit = %.6f\n",
              res$converged, res$iterations, res$max_d, res$criterion))
  cat("support point                    weight\n")
  for (i in seq_len(nrow(S))) {
    pt <- paste(sprintf("%.4f", S[i, ]), collapse = ", ")
    cat(sprintf("  %-30s %.6f\n", pt, w[i]))
  }
  invisible(res)
}

#' Pretty-print a design returned by \code{\link{optimal_design}}.
#'
#' @param res a result list from \code{optimal_design()}.
#' @param title a title for the printout.
#' @return \code{res}, invisibly.
#' @export
print_result <- function(res, title = "Optimal design") {
  # an exact design has its own printout (efficiency %, integer counts)
  if (inherits(res, "exact_design")) return(print(res))
  S <- as.matrix(res$support); storage.mode(S) <- "double"
  perm <- do.call(order, lapply(seq_len(ncol(S)), function(j) S[, j]))
  S <- S[perm, , drop = FALSE]; w <- res$weights[perm]
  cat(sprintf("\n=== %s ===\n", title))
  cat(sprintf("|support| = %d    crit = %.10f    max sensitivity = %.3e  (0 at the optimum)    total time = %.3f s\n",
              nrow(S), res$criterion, res$max_d, res$total_time))
  if (identical(res$method, "continuous"))
    cat(sprintf(paste0("continuous search over the design region (no grid): %d iteration(s); ",
                       "the max sensitivity is over the points examined (multi-start search",
                       "%s)\n"),
                res$iterations,
                if (isTRUE(res$n_audit > 0))
                  sprintf(" + a random audit of %d points", res$n_audit) else ""))
  if (!is.null(res$elfving_bound) && is.finite(res$elfving_bound)) {
    value <- if (isTRUE(res$p == 0L)) exp(res$criterion) else res$criterion
    # is the design's information matrix singular?  Then the sensitivity
    # function could not have certified it, and Elfving's bound is what did.
    singular <- FALSE
    if (!is.null(res$information)) {
      ev <- tryCatch(eigen(as.matrix(res$information), symmetric = TRUE,
                           only.values = TRUE)$values, error = function(e) NULL)
      singular <- !is.null(ev) && min(ev) <= 1e-10 * max(ev)
    }
    cat(sprintf(paste0("c-optimality (one linear combination): Elfving's bound on the optimal ",
                       "value = %.6g; this design = %.6g; gap = %.3g -> efficiency >= %.4f%%",
                       "%s\n"),
                res$elfving_bound, value, res$elfving_gap,
                100 * res$efficiency_lower_bound,
                if (isTRUE(res$converged) && singular)
                  " (certified; the design is singular, which the sensitivity function cannot certify)"
                else if (isTRUE(res$converged)) " (certified)"
                else ""))
  }
  cat("support point                    weight\n")
  for (i in seq_len(nrow(S))) {
    pt <- paste(.fmt_support_row(S[i, ], res$is_factor), collapse = ", ")
    cat(sprintf("  %-32s %.6f\n", pt, w[i]))
  }
  invisible(res)
}
