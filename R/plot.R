# ===========================================================================
# plot.R -- simple base-graphics visualisation of a design.  No extra package
# dependency, so both the console and the Shiny app can call it.
# ===========================================================================

# Core plotter.  `support` is m x d (one design point per row); `amount` is the
# per-point weight (approximate design) or integer count (exact design).
.plot_design <- function(support, amount, is_factor = NULL, cov_names = NULL,
                         amount_label = "weight", main = "Design") {
  support <- as.matrix(support)
  d <- ncol(support)
  amount <- as.numeric(amount)
  if (is.null(cov_names) || length(cov_names) != d)
    cov_names <- if (!is.null(colnames(support))) colnames(support)
                 else paste0("x", seq_len(d))
  if (is.null(is_factor)) is_factor <- rep(FALSE, d)

  if (d == 1L) {
    # stems: covariate value on x, amount as height.
    ord <- order(support[, 1])
    xv <- support[ord, 1]; av <- amount[ord]
    plot(xv, av, type = "h", lwd = 3, col = "#2c7fb8",
         xlab = cov_names[1], ylab = amount_label, main = main,
         ylim = c(0, max(av) * 1.1))
    points(xv, av, pch = 19, col = "#2c7fb8")
    text(xv, av, labels = sprintf("%.3g", av), pos = 3, cex = 0.8)
  } else if (d == 2L) {
    # bubble scatter: point position = (x1, x2), area ~ amount.
    rad <- 3 * sqrt(amount / max(amount)) + 0.8
    pad <- function(v) { r <- range(v); d <- diff(r); if (d == 0) d <- 1
                         c(r[1] - 0.1 * d, r[2] + 0.1 * d) }
    plot(support[, 1], support[, 2], cex = rad, pch = 21,
         bg = "#7fcdbb", col = "#2c7fb8", lwd = 1.5,
         xlim = pad(support[, 1]), ylim = pad(support[, 2]),
         xlab = cov_names[1], ylab = cov_names[2], main = main)
    text(support[, 1], support[, 2], labels = sprintf("%.3g", amount),
         cex = 0.75, pos = 3)
  } else {
    # >= 3 covariates: one bar per support point, labelled by its coordinates.
    labs <- apply(support, 1, function(r)
      paste(sprintf("%.3g", r), collapse = ","))
    ord <- order(-amount)
    op <- par(mar = c(8, 4, 4, 2) + 0.1); on.exit(par(op), add = TRUE)
    barplot(amount[ord], names.arg = labs[ord], las = 2, col = "#7fcdbb",
            border = "#2c7fb8", ylab = amount_label, main = main, cex.names = 0.7)
    mtext(paste0("support point (", paste(cov_names, collapse = ", "), ")"),
          side = 1, line = 6, cex = 0.8)
  }
  invisible(NULL)
}

#' Plot a design returned by \code{\link{optimal_design}} or
#' \code{\link{exact_design}}.
#'
#' Draws a quick visual of the design: a stem plot of weight/count vs the single
#' covariate (one covariate), a bubble scatter with area proportional to
#' weight/count (two covariates), or a labelled bar per support point (three or
#' more). Uses base graphics only.
#'
#' @param res a result from \code{\link{optimal_design}} (uses \code{weights})
#'   or \code{\link{exact_design}} (uses integer \code{counts}).
#' @param cov_names optional covariate labels (defaults to the support column
#'   names, or \code{x1, x2, ...}).
#' @param main plot title.
#' @return \code{NULL}, invisibly (called for the plot).
#' @seealso \code{\link{optimal_design}}, \code{\link{exact_design}}.
#' @export
plot_design <- function(res, cov_names = NULL, main = "Design") {
  if (!is.null(res$counts)) { amount <- res$counts; lab <- "count" }
  else                      { amount <- res$weights; lab <- "weight" }
  .plot_design(res$support, amount, is_factor = res$is_factor,
               cov_names = cov_names, amount_label = lab, main = main)
}
