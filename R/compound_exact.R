# ===========================================================================
# compound_exact.R -- exact (integer-run) COMPOUND designs.
#
# The compound counterpart of exact.R, and it follows the same four steps:
#
#   1. solve the approximate compound optimum (the reference design);
#   2. round its weights to integer counts by largest-remainder apportionment;
#   3. repair the total to exactly n using the COMPOUND sensitivity function --
#      peel a run off the least informative support point, or add one to the
#      most informative candidate;
#   4. improve by random exchanges, accepting only a strict improvement.
#
# The one structural difference from exact.R is the direction: Phi_p is
# MINIMISED there, whereas Psi_alpha is MAXIMISED here, so every comparison is
# reversed and the efficiency is a ratio the other way up
# (Psi_alpha(exact) / Psi_alpha(approximate), in (0, 1]).
#
# Speed matters in step 4, so the per-candidate information is built ONCE and
# every criterion/sensitivity evaluation goes straight to the C++ core by index
# -- exactly the arrangement exact.R uses.
# ===========================================================================

#' Exact (integer-run) compound optimal design.
#'
#' The exact-design counterpart of \code{\link{compound_design}}: allocates a
#' fixed number of runs \code{n} over design points so that the compound
#' criterion \eqn{\Psi_\alpha} is as large as possible. Supply the components
#' and the design region exactly as for \code{\link{compound_design}}, plus
#' \code{n}.
#'
#' @param n integer sample size --- the number of runs to allocate.
#' @param components,alpha,efficiency,psi_star,reference_bound the compound
#'   criterion, exactly as in \code{\link{compound_design}}.
#' @param design_box,step_sequence,candidate_set,factor_levels the design
#'   region, as in \code{\link{compound_design}}.
#' @param xi0_points,xi0_weights,n0,n1 an existing design to augment. When an
#'   existing design is given and \code{n1} is not set explicitly it defaults to
#'   \code{n}, so the design is optimised for the balance actually realised.
#' @param max_exchange number of random exchanges attempted (default 1000).
#' @param seed optional integer seed, for reproducibility.
#' @param snap_support for the \code{design_box} path: if \code{TRUE} (default)
#'   the approximate support is snapped to the finest grid; if \code{FALSE} the
#'   off-grid support points are appended to the candidate set so the exact
#'   design may use them.
#' @param max_iter,eps0,accept_tol,verbose passed to the internal
#'   \code{\link{compound_design}} call for the approximate reference.
#' @return a list of class \code{"compound_exact_design"} with \code{support},
#'   \code{counts} (integers summing to \code{n}), \code{weights}
#'   (\code{counts / n}), \code{criterion} (the compound value of the exact
#'   design), \code{criterion_approx} (the approximate optimum it is measured
#'   against), \code{efficiency_exact_lower_bound} (a guaranteed lower bound
#'   on the exact design's efficiency relative to the TRUE compound optimum:
#'   their ratio times the approximate design's own bound
#'   \code{criterion_bound}),
#'   \code{psi}, \code{psi_star}, \code{component_criterion} and
#'   \code{component_criterion_star} (the per-component values on the scale
#'   \code{\link{optimal_design}} reports, see \code{\link{compound_design}}),
#'   \code{efficiency_lower_bound} (per component, relative to the TRUE
#'   component optimum),
#'   \code{cross_efficiency}, \code{alpha}, \code{at}, \code{information},
#'   \code{n}, \code{exchanges}, \code{n_candidates} (the size of the candidate
#'   set the exchanges searched) and \code{approx} (the full approximate
#'   result).
#' @details The reported \code{efficiency_exact_lower_bound} is a LOWER bound on the exact
#'   design's efficiency: it compares against the approximate optimum, which is
#'   at least as good as any exact design, so the true efficiency is at least
#'   this value.
#'
#'   With a multistage \code{step_sequence} the exchange step searches a
#'   \emph{coarse} grid over the whole box (the coarsest step in the sequence)
#'   together with a \emph{fine} neighbourhood of the approximate support at
#'   the finest step --- the same candidate set \code{\link{exact_design}}
#'   builds. The finest grid is never materialised over the whole region, whose
#'   size grows as the reciprocal of the final step raised to the number of
#'   continuous covariates. A single-step \code{step_sequence} is unaffected:
#'   coarse and fine coincide and the candidate set is the full grid, as
#'   before. Supply \code{candidate_set} to control the points exactly.
#' @seealso \code{\link{compound_design}}, \code{\link{compound_criterion}},
#'   \code{\link{exact_design}} (the single-criterion version).
#' @examples
#' \dontrun{
#' f1 <- function(x, th) { q <- c(1, x[1]); e <- sum(q * th)
#'                         (exp(e / 2) / (1 + exp(e))) * q }
#' f2 <- function(x, th) { q <- c(1, x[1], x[1]^2); e <- sum(q * th)
#'                         (exp(e / 2) / (1 + exp(e))) * q }
#' r <- compound_exact_design(
#'   n = 30,
#'   components = list(list(info_vector = f1, theta = c(0.5, 1), p = 0),
#'                     list(info_vector = f2, theta = c(0.5, 1, -0.5), p = 0)),
#'   design_box = list(c(-2, 2)), step_sequence = c(0.5, 0.1), seed = 1)
#' print(r)
#' }
#' @export
compound_exact_design <- function(n, components, alpha = NULL,
                                  design_box = NULL, step_sequence = NULL,
                                  candidate_set = NULL, factor_levels = NULL,
                                  efficiency = TRUE, psi_star = NULL,
                                  reference_bound = NULL,
                                  xi0_points = NULL, xi0_weights = numeric(0),
                                  n0 = 0, n1 = 1,
                                  max_exchange = 1000L, seed = NULL,
                                  snap_support = TRUE,
                                  max_iter = 100L, eps0 = 1e-6,
                                  accept_tol = 1e-9, verbose = FALSE) {
  t_start <- proc.time()[3]

  n <- as.integer(round(n))
  if (length(n) != 1L || is.na(n) || n < 1L)
    stop("'n' must be a positive integer sample size.", call. = FALSE)

  # As in exact_design(): with an existing design, the new stage should weigh as
  # the n runs actually being added, so n1 defaults to n.
  if (missing(n1) && !is.null(xi0_points) && n0 > 0) {
    n1 <- n
    if (isTRUE(verbose))
      cat(sprintf(paste0("  compound_exact_design(): n1 not supplied; using ",
                         "n1 = n = %d for the n0 : n1 balance.\n"), n))
  }

  # ---- 1. the approximate compound optimum (the reference) ----------------
  ap <- compound_design(components = components, alpha = alpha,
                        design_box = design_box, step_sequence = step_sequence,
                        candidate_set = candidate_set,
                        factor_levels = factor_levels,
                        efficiency = efficiency, psi_star = psi_star,
                        reference_bound = reference_bound,
                        xi0_points = xi0_points, xi0_weights = xi0_weights,
                        n0 = n0, n1 = n1,
                        max_iter = max_iter, eps0 = eps0,
                        accept_tol = accept_tol, verbose = verbose)

  # ---- 2. rebuild the model quantities on ONE candidate set ---------------
  su <- .cmp_setup(components, design_box, candidate_set, factor_levels,
                   xi0_points, xi0_weights, n0, n1)
  comps <- su$comps
  J  <- length(comps)
  ms <- max(vapply(comps, function(cc) cc$min_support, numeric(1)))
  if (n < ms)
    stop(sprintf(paste0("n = %d is below the %d run(s) needed for a ",
                        "non-singular information matrix in every component."),
                 n, ms), call. = FALSE)

  use_set <- !is.null(candidate_set)
  if (use_set) {
    X <- as.matrix(candidate_set); storage.mode(X) <- "double"
  } else {
    # Candidate set: a coarse full-box grid (global reach for the add /
    # exchange steps) UNION a fine neighbourhood of the approximate support at
    # the finest step -- the same construction exact_design() uses, and for the
    # same reason: materialising the finest grid over the WHOLE box costs the
    # cube (or worse) of the final step, and every one of those points then
    # carries per-component information.  A step sequence therefore stops being
    # something you pay for here.  With a single step coarse == fine, so this
    # reduces exactly to the full grid and nothing changes.
    stage_by <- .normalize_step_sequence(
      if (is.null(step_sequence)) numeric(0) else step_sequence, su$is_factor)
    if (!length(stage_by)) {
      X <- .factor_make_grid(su$meta$lo, su$meta$hi,
                             rep(1, length(su$is_factor)), su$is_factor,
                             su$nlevels)
    } else {
      coarse <- do.call(pmax, stage_by)
      fine   <- do.call(pmin, stage_by)
      # per-covariate neighbourhood radius: the second-finest step, or twice
      # the finest when that covariate has only one distinct step
      radius <- vapply(seq_along(fine), function(d) {
        u <- sort(unique(vapply(stage_by, `[`, numeric(1), d)))
        if (length(u) >= 2L) u[2] else 2 * u[1]
      }, numeric(1))
      X <- unique(rbind(
        .factor_make_grid(su$meta$lo, su$meta$hi, coarse, su$is_factor,
                          su$nlevels),
        .factor_refined_grid(su$meta$lo, su$meta$hi, ap$support, radius, fine,
                             su$is_factor, su$nlevels)))
    }
    if (!isTRUE(snap_support)) {
      # keep the approximate support exactly where it is, off-grid and all
      off <- as.matrix(ap$support)
      keep <- vapply(seq_len(nrow(off)), function(i) {
        d <- rowSums((X - matrix(off[i, ], nrow(X), ncol(X), byrow = TRUE))^2)
        min(d) > 1e-20
      }, logical(1))
      if (any(keep)) X <- rbind(X, off[keep, , drop = FALSE])
    }
  }
  Ncand <- nrow(X)

  # per-component candidate information, built ONCE
  data <- lapply(comps, function(cc) .cmp_data(cc, X))
  im <- as.integer(vapply(comps, function(cc) cc$info_mode, numeric(1)))
  wb <- lapply(comps, function(cc) cc$wb)
  i0 <- lapply(comps, function(cc) cc$infor0)
  pp <- as.integer(vapply(comps, function(cc) cc$p, numeric(1)))
  at <- as.numeric(ap$at)

  # Psi_alpha and the compound sensitivity at an integer allocation
  crit_of <- function(cnt) {
    idx <- which(cnt > 0L)
    if (!length(idx)) return(-Inf)
    v <- tryCatch(compound_criterion_cpp(data, im, wb, i0, pp, at,
                                         as.integer(idx),
                                         as.numeric(cnt[idx] / sum(cnt))),
                  error = function(e) NA_real_)
    if (!is.finite(v)) -Inf else v
  }
  dir_of <- function(cnt) {
    idx <- which(cnt > 0L)
    compound_dirderiv_cpp(data, im, wb, i0, pp, at, as.integer(idx),
                          as.numeric(cnt[idx] / sum(cnt)))
  }

  # ---- 3. round the approximate weights to integer counts -----------------
  a_idx <- .nearest_idx(X, ap$support)
  aw <- .combine_idx_w(a_idx, as.numeric(ap$weights))
  cnt <- integer(Ncand)
  cnt[aw$idx] <- .apportion(aw$w / sum(aw$w), n)

  # ---- 4. repair the total to exactly n using the sensitivity -------------
  total <- sum(cnt)
  while (total > n) {                          # drop the least informative run
    d   <- dir_of(cnt)
    sup <- which(cnt > 0L)
    ok  <- sup[order(d[sup])]                  # smallest sensitivity first
    done <- FALSE
    for (i in ok) {
      if (cnt[i] == 1L && sum(cnt > 0L) - 1L < ms) next
      prop <- cnt; prop[i] <- prop[i] - 1L
      if (is.finite(crit_of(prop))) { cnt <- prop; done <- TRUE; break }
    }
    if (!done)
      stop("compound_exact_design(): cannot reduce the design to n without ",
           "making a component singular.", call. = FALSE)
    total <- total - 1L
  }
  while (total < n) {                          # add the most informative run
    cnt[which.max(dir_of(cnt))] <- cnt[which.max(dir_of(cnt))] + 1L
    total <- total + 1L
  }

  # ---- 5. random exchanges (accept iff STRICTLY better) -------------------
  # Psi_alpha is maximised, so the acceptance test is reversed relative to
  # exact_design(), which minimises Phi_p.
  if (!is.null(seed)) set.seed(as.integer(seed))
  cur <- crit_of(cnt); accepted <- 0L
  if (max_exchange > 0L) for (t in seq_len(as.integer(max_exchange))) {
    sup <- which(cnt > 0L)
    i <- sup[sample.int(length(sup), 1L)]
    j <- sample.int(Ncand, 1L)
    if (j == i) next
    if (cnt[i] == 1L && sum(cnt > 0L) - 1L < ms) next
    prop <- cnt; prop[i] <- prop[i] - 1L; prop[j] <- prop[j] + 1L
    nc <- crit_of(prop)
    if (is.finite(nc) && nc > cur + 1e-12) {
      cnt <- prop; cur <- nc; accepted <- accepted + 1L
    }
  }

  # ---- 6. assemble the result ---------------------------------------------
  idx <- sort(which(cnt > 0L))
  support <- X[idx, , drop = FALSE]
  counts  <- as.integer(cnt[idx])
  w       <- counts / n
  crit_exact <- crit_of(cnt)

  psi <- compound_psi_cpp(data, im, wb, i0, pp, at, as.integer(idx),
                          as.numeric(w))
  nms <- names(ap$psi)
  eff <- if (isTRUE(ap$efficiency_weighted)) psi / ap$psi_star
         else rep(NA_real_, J)

  Mlist <- lapply(seq_len(J), function(j) {
    cc <- comps[[j]]
    M <- cc$infor0 + .opt_infor_from_support(support, w, cc$b, cc$info_mode,
                                             cc$info_vector, cc$info_matrix,
                                             cc$theta_use, cc$k)
    dimnames(M) <- list(cc$coef_names, cc$coef_names)
    M
  })

  # the cross-efficiency table, with THIS exact design as the last row
  cross <- NULL
  if (!is.null(ap$reference_designs) && isTRUE(ap$efficiency_weighted)) {
    rows <- lapply(ap$reference_designs, function(d)
      .cmp_eval(comps, at, d$support, d$weights)$psi / ap$psi_star)
    rows[[J + 1L]] <- as.numeric(eff)
    cross <- do.call(rbind, rows)
    dimnames(cross) <- list(c(paste("optimal for", nms), "THIS DESIGN"), nms)
    # lower bounds relative to each component's TRUE optimum
    rb <- as.numeric(ap$reference_bound)
    cross <- sweep(pmin(cross, 1), 2, ifelse(is.finite(rb), rb, 1), `*`)
  }

  out <- list(support = support, counts = counts, weights = w,
              criterion = crit_exact,
              criterion_approx = ap$criterion,
              # relative to the TRUE compound optimum: the ratio to the approximate
              # compound design times that design's certified bound (eq. 27)
              efficiency_exact_lower_bound = .cmp_certify(crit_exact / ap$criterion,
                                                          ap$criterion_bound),
              psi = stats::setNames(as.numeric(psi), nms),
              psi_star = ap$psi_star,
              # the same values on the scale optimal_design() reports
              component_criterion      = stats::setNames(.cmp_single_scale(psi, ap$p), nms),
              component_criterion_star = stats::setNames(.cmp_single_scale(ap$psi_star, ap$p), nms),
              # ONE efficiency per component, relative to the TRUE component
              # optimum (ratio x the reference design's bound, Thm 4.5 / 4.6)
              efficiency_lower_bound = stats::setNames(.cmp_certify(eff, ap$reference_bound), nms),
              reference_bound = ap$reference_bound,
              cross_efficiency = cross,
              alpha = ap$alpha, at = ap$at,
              p = ap$p,
              information = stats::setNames(Mlist, nms),
              n = n, n0 = as.integer(n0), exchanges = accepted,
              n_candidates = Ncand,        # the set the exchanges searched
              efficiency_weighted = ap$efficiency_weighted,
              is_factor = su$is_factor,
              approx = ap,
              total_time = proc.time()[3] - t_start)
  class(out) <- "compound_exact_design"
  out
}

#' @param x a \code{"compound_exact_design"} result.
#' @param ... ignored.
#' @rdname compound_exact_design
#' @export
print.compound_exact_design <- function(x, ...) {
  cat("Exact compound optimal design\n")
  cat(sprintf("  runs       : %d over %d support point(s)\n",
              x$n, nrow(x$support)))
  cat(sprintf("  criterion  : Psi_alpha = %.8f%s\n", x$criterion,
              if (isTRUE(x$efficiency_weighted))
                "   (weighted average efficiency, in (0,1])" else ""))
  cat(sprintf("  efficiency : >= %.4f%% of the compound optimum   (approximate compound design: %.8f)\n",
              100 * x$efficiency_exact_lower_bound, x$criterion_approx))
  cat(sprintf("  exchanges  : %d accepted\n", x$exchanges))
  cat("\n  per-component criterion values, on the scale optimal_design() reports",
      "\n  (D: log det Sigma / v;  A: tr(Sigma) / v;  smaller is better)\n")
  tab <- data.frame(alpha = round(as.numeric(x$alpha), 4),
                    p = ifelse(x$p == 0, "D", "A"),
                    criterion = signif(as.numeric(x$component_criterion), 7),
                    stringsAsFactors = FALSE)
  if (isTRUE(x$efficiency_weighted)) {
    tab$optimal    <- signif(as.numeric(x$component_criterion_star), 7)
    tab$efficiency_lower_bound <- round(as.numeric(x$efficiency_lower_bound), 6)
  }
  rownames(tab) <- names(x$psi)
  print(tab)
  if (!is.null(x$cross_efficiency)) {
    cat("\n  efficiency of each design under every component",
        "\n  (rows: design;  columns: component;  lower bounds relative to each",
        "\n   component's true optimum)\n")
    print(round(x$cross_efficiency, 3))
  }
  cat("\n  design (support points and runs):\n")
  d <- data.frame(x$support, count = x$counts)
  names(d) <- c(paste0("x", seq_len(ncol(x$support))), "count")
  print(d, row.names = FALSE)
  invisible(x)
}
