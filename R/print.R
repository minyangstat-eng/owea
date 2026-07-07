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
  cat("support point                    weight\n")
  for (i in seq_len(nrow(S))) {
    pt <- paste(.fmt_support_row(S[i, ], res$is_factor), collapse = ", ")
    cat(sprintf("  %-32s %.6f\n", pt, w[i]))
  }
  invisible(res)
}
