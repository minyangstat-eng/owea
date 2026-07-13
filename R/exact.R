# ===========================================================================
# exact.R -- exact (integer) optimal designs built on the approximate optimum.
#
# Given a sample size N, exact_design() converts the approximate Phi_p-optimal
# design (continuous weights) into a highly efficient EXACT design: an integer
# allocation n_1, ..., n_m of N runs over design points, with sum n_i = N.
#
# The construction follows Yang, Biedermann & Tang (2013, JASA 108:1411-1420):
#   1. Round the approximate weights to integer counts (Hamilton / largest
#      remainder apportionment).
#   2. Adjust the total to exactly N using the sensitivity (directional
#      derivative) function: when too large, peel a replication off the most
#      INEFFICIENT support point (smallest sensitivity); when too small, add a
#      replication to the most INFORMATIVE candidate (largest sensitivity).
#   3. Improve by random exchanges: move one run from a support point to
#      another candidate, accept iff the criterion strictly improves.
#   4. Report the efficiency of the exact design relative to the approximate
#      optimum (the exact design is "projected" onto the approximate design
#      space, i.e. evaluated as an approximate design with weights n_i / N).
# ===========================================================================

# -- internal helpers -------------------------------------------------------

# Hamilton (largest-remainder) apportionment: round the weights w (summing to
# 1) to integers summing to N.  Floor everyone, then hand the N - sum(floor)
# leftover units to the largest fractional remainders.  Ties break by lowest
# index (stable order), so the result is deterministic.
.apportion <- function(w, N) {
  raw  <- w * N
  base <- floor(raw)
  rem  <- N - sum(base)
  if (rem > 0) {
    ord <- order(raw - base, decreasing = TRUE)   # stable: ties -> lower index
    take <- ord[seq_len(rem)]
    base[take] <- base[take] + 1
  }
  as.integer(base)
}

# Combine duplicate indices, summing their weights (used when several
# approximate support points snap to the same grid column).
.combine_idx_w <- function(idx, w) {
  u <- sort(unique(idx))
  list(idx = u, w = vapply(u, function(j) sum(w[idx == j]), numeric(1)))
}

#' Construct an exact (integer) Phi_p-optimal design for a fixed sample size.
#'
#' The exact-design counterpart of \code{\link{optimal_design}}. Supply the
#' model and design region exactly as for \code{optimal_design()} plus an
#' integer sample size \code{n}; \code{exact_design()} computes the approximate
#' optimum internally, rounds it to an integer allocation, repairs the total to
#' exactly \code{n} and improves it by random exchanges, and reports the
#' efficiency of the resulting exact design relative to the approximate optimum.
#'
#' @param n integer sample size \eqn{N} (the number of runs to allocate).
#' @param design_box,step_sequence design region and grid steps, as in
#'   \code{\link{optimal_design}}; \code{design_box} may include FACTOR
#'   covariates (a single integer \code{L >= 2}) alongside continuous
#'   \code{c(lo, hi)} entries. Omit when \code{candidate_set} is given.
#' @param factor_levels for the \code{candidate_set} path, per-column factor
#'   markers, as in \code{\link{optimal_design}}. Factor covariates are passed
#'   to the model as integer levels and printed as integers.
#' @param candidate_set optional \eqn{n \times N} matrix of candidate points
#'   (one per row); the exact design is found on this fixed set.
#' @param info_vector,info_matrix the model (supply exactly one), as in
#'   \code{\link{optimal_design}}.
#' @param theta parameter values (required for a \code{function(x, theta)}
#'   model, \code{grad_g}, or an existing design).
#' @param link,f,x,fx,ff,xx,intercept,coding,ncat formula-style model spec, as in
#'   \code{\link{optimal_design}} (an alternative to
#'   \code{info_vector}/\code{info_matrix} for linear / logistic / log-linear /
#'   multinomial / cumulative models; \code{ff} = factor \eqn{\times} factor
#'   interactions, \code{ncat} = number of response categories for
#'   \code{"multinomial"}/\code{"cumulative"}). When a spec needs \code{theta}
#'   and none is given, it is drawn from \eqn{N(0,1)} once, with a warning, and
#'   returned in \code{theta}. The result gains \code{coef_names} / \code{link}
#'   (see \code{\link{model_summary}}).
#' @param p criterion: 0 = D-optimality or 1 = A-optimality (only these two are
#'   supported).
#' @param wb,subset,grad_g quantity of interest (at most one); default identity.
#' @param xi0_points,xi0_weights,n0,n1 existing-design support, weights and the
#'   existing / new sample sizes (multistage). Here \code{n} is the EXACT size
#'   of the NEW stage. \code{xi0_weights} are proportions (summing to 1, not
#'   counts; normalized with a warning otherwise) and \code{n0} is the existing
#'   sample size; the existing : new balance is \code{n0 : n1}. When an existing
#'   design is given (\code{xi0_points} with \code{n0 > 0}) and \code{n1} is not
#'   supplied, \code{n1} defaults to \code{n} so the new stage is weighted as
#'   the \code{n} runs being allocated; supplying \code{n1 != n} is allowed but
#'   warns, since the design is then optimized for a different existing : new
#'   balance than the \code{n} runs realized.
#' @param max_exchange maximum number of random exchanges in the improvement
#'   step (default 1000).
#' @param seed optional integer seed making the random exchanges reproducible.
#' @param snap_support (\code{design_box} path only) if \code{TRUE} (default),
#'   snap the approximate support to the candidate grid and re-optimise the
#'   reference weights on it, keeping a clean grid; if \code{FALSE}, append the
#'   (possibly off-grid) approximate support to the grid so it is represented
#'   exactly. The candidate grid is a coarse full-box grid (at the coarsest
#'   \code{step_sequence} step) plus a fine neighbourhood of the approximate
#'   support at the finest step -- fine resolution where it matters without
#'   materialising the full finest grid over the whole box.
#' @param init_method,auto_warm_start,max_iter,eps0,accept_tol passed through to
#'   the internal \code{\link{optimal_design}} call.
#' @param solver,ma_max_iter passed through to the internal
#'   \code{\link{optimal_design}} call: \code{solver = "MA"} (alias
#'   \code{"multiplicative"}) solves the reference APPROXIMATE design with the
#'   multiplicative algorithm instead of the OWEA exchange engine (see
#'   \code{\link{optimal_design}}); the exact-design rounding and exchange steps
#'   are unchanged. Default \code{solver = "owea"}.
#' @param check_global,global_step,global_max_points passed through to the
#'   internal \code{\link{optimal_design}} call (\code{design_box} path):
#'   verify the approximate design over a fine grid spanning the whole box and
#'   report its global optimality. See \code{\link{optimal_design}}.
#' @param merge,merge_factor,merge_atol passed through to the internal
#'   \code{\link{optimal_design}} call (merging of the approximate design).
#' @param verbose print progress.
#' @return an object of class \code{"exact_design"}: a list with \code{support}
#'   (the exact support points, one per row; zero-count points dropped),
#'   \code{counts} (the integer replications \eqn{n_i}, summing to \code{n}),
#'   \code{weights} (\eqn{n_i / n}), \code{n}, \code{criterion} (the exact
#'   design's normalised per-sample Phi_p value), \code{criterion_total} (the
#'   Phi_p value of the TOTAL information matrix -- the criterion for all
#'   \code{n} runs at once), \code{information} (the counts-based combined TOTAL
#'   Fisher information \eqn{\sum_i n_i M(x_i)}, plus \eqn{n_0 I_{\xi_0}} when an
#'   existing design is supplied), \code{efficiency} (a guaranteed LOWER
#'   BOUND on the exact design's efficiency, in \eqn{(0, 1]}: it is measured
#'   against the approximate optimum, which is at least as good as any exact
#'   design, so the true efficiency is at least this), \code{criterion_approx},
#'   \code{approx}
#'   (the approximate design used, including its \code{global_max_d} /
#'   \code{global_check} when \code{check_global = TRUE}), \code{p},
#'   \code{candidate_set}, \code{n_exchange_accepted}, \code{converged}
#'   (inherited from the approximate solve) and \code{total_time}.
#' @seealso \code{\link{optimal_design}} for the approximate design.
#' @export
exact_design <- function(n,
                         design_box = NULL, step_sequence = NULL,
                         candidate_set = NULL,
                         info_vector = NULL, info_matrix = NULL,
                         theta = NULL, p = 0L,
                         link = NULL, f = NULL, x = NULL, fx = NULL,
                         xx = NULL, ff = NULL, intercept = TRUE,
                         coding = "zero-sum", ncat = NULL,
                         wb = NULL, subset = NULL, grad_g = NULL,
                         xi0_points = NULL, xi0_weights = numeric(0),
                         n0 = 0, n1 = 1,
                         factor_levels = NULL,
                         max_exchange = 1000L, seed = NULL,
                         snap_support = TRUE,
                         init_method = "auto", auto_warm_start = TRUE,
                         solver = "owea", ma_max_iter = 100L,
                         check_global = FALSE, global_step = NULL,
                         global_max_points = 1e6,
                         max_iter = 100L, eps0 = 1e-6, accept_tol = 1e-9,
                         merge = FALSE, merge_factor = 1.5, merge_atol = NULL,
                         verbose = FALSE) {
  t_start <- proc.time()[3]

  n <- as.integer(round(n))
  if (length(n) != 1L || is.na(n) || n < 1L)
    stop("'n' must be a positive integer sample size.", call. = FALSE)
  p <- .check_criterion(p)

  # When an existing design is given and n1 was not set explicitly, default the
  # new-stage size n1 to the exact sample size n (the existing : new balance is
  # n0 : n1, so the new stage should weigh as the n runs actually being added).
  if (missing(n1) && !is.null(xi0_points) && n0 > 0) {
    n1 <- n
    if (isTRUE(verbose))
      cat(sprintf(paste0("  exact_design(): n1 not supplied; using n1 = n = %d ",
                        "for the existing : new (n0 : n1) balance.\n"), n))
  } else if (!missing(n1) && !is.null(xi0_points) && n0 > 0 && n1 != n) {
    warning(sprintf(paste0("n1 = %d differs from the exact sample size n = %d. ",
                          "The design is optimized (and its criterion reported) ",
                          "for an existing : new balance of n0 : n1 = %d : %d, ",
                          "but %d runs are then allocated -- so the integer ",
                          "design is optimal for a DIFFERENT balance than the ",
                          "one realized. Set n1 = n (or omit n1) unless this is ",
                          "intentional."),
                   n1, n, n0, n1, n), call. = FALSE)
  }

  # Resolve a formula-style model spec ONCE here (so theta is drawn once and
  # reported), then delegate to optimal_design() with the built info_vector.
  spec <- .resolve_model_spec(link, f, x, fx, xx, intercept, coding,
                              design_box, candidate_set, factor_levels,
                              info_vector, info_matrix, theta, ff = ff,
                              ncat = ncat)
  info_vector <- spec$info_vector; info_matrix <- spec$info_matrix
  theta <- spec$theta
  if (spec$spec_given && is.null(factor_levels)) factor_levels <- spec$factor_levels
  coef_names <- spec$coef_names; model_link <- spec$link

  if (is.null(info_matrix) && is.null(info_vector))
    stop("Supply 'info_matrix', 'info_vector', or a model spec ('link' + terms).",
         call. = FALSE)
  if (!is.null(info_matrix) && !is.null(info_vector))
    stop("Supply only one of 'info_matrix' / 'info_vector'.", call. = FALSE)
  info_mode <- if (is.null(info_vector)) 1L else 0L

  use_set <- !is.null(candidate_set)
  if (!use_set && is.null(design_box))
    stop("Supply 'candidate_set', or 'design_box' (plus 'step_sequence' for any ",
         "continuous covariates).", call. = FALSE)
  # step_sequence may be omitted for an all-factor design_box.
  if (!use_set && is.null(step_sequence)) step_sequence <- numeric(0)

  # ---- 1. approximate optimum (the reference design) -------------------
  # Muffle only optimal_design()'s informational "multistage convergence is
  # LOCAL" warning; genuine non-convergence warnings still surface.
  appr <- withCallingHandlers(
    optimal_design(design_box = design_box, step_sequence = step_sequence,
                   info_vector = info_vector, info_matrix = info_matrix,
                   theta = theta, p = p, wb = wb, subset = subset,
                   grad_g = grad_g, xi0_points = xi0_points,
                   xi0_weights = xi0_weights, n0 = n0, n1 = n1,
                   candidate_set = candidate_set,
                   factor_levels = factor_levels,
                   merge = merge,
                   merge_factor = merge_factor, merge_atol = merge_atol,
                   init_method = init_method, auto_warm_start = auto_warm_start,
                   solver = solver, ma_max_iter = ma_max_iter,
                   check_global = check_global, global_step = global_step,
                   global_max_points = global_max_points,
                   max_iter = max_iter, eps0 = eps0, accept_tol = accept_tol,
                   verbose = verbose),
    warning = function(w) {
      # muffle informational warnings that exact_design re-emits from its own
      # setup (the existing-design ones below) or that are not actionable here.
      msg <- conditionMessage(w)
      if (grepl("multistage convergence is LOCAL", msg, fixed = TRUE) ||
          grepl("existing design is IGNORED",      msg, fixed = TRUE) ||
          grepl("normalizing them to proportions", msg, fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  if (!isTRUE(appr$converged))
    warning("exact_design(): the approximate design did not converge; the ",
            "efficiency reference is not a certified optimum.", call. = FALSE)

  # ---- 2. setup: rebuild the model quantities (mirrors optimal_design) --
  if (info_mode == 0L) {
    if (.info_vec_needs_theta(info_vector) && is.null(theta))
      stop("'theta' is required when 'info_vector' is function(x, theta).",
           call. = FALSE)
    info_vector <- .normalize_info_vector(info_vector, theta)
  }
  info_matrix <- .normalize_info_matrix(info_matrix, theta)

  # factor (categorical) covariate metadata (mirrors optimal_design; factors
  # are carried in level space and passed to the model as integer levels).
  meta <- if (use_set)
            .factor_levels_to_meta(factor_levels, ncol(as.matrix(candidate_set)))
          else .parse_design_box(design_box)
  is_factor <- meta$is_factor; nlevels <- meta$nlevels
  if (use_set)
    .validate_factor_columns(candidate_set, is_factor, nlevels, "candidate_set")
  if (!is.null(xi0_points))
    .validate_factor_columns(xi0_points, is_factor, nlevels, "xi0_points")

  probe <- if (use_set) as.numeric(as.matrix(candidate_set)[1, ])
           else .factor_probe_point(meta$lo, meta$hi, is_factor)
  k <- if (info_mode == 1L) nrow(as.matrix(info_matrix(probe)))
       else length(as.numeric(info_vector(probe, theta)))
  theta_use <- if (is.null(theta)) rep(0.0, k) else as.numeric(theta)

  wb_use <- .wb_from(wb, subset, grad_g, theta_use, k)
  i0 <- .make_infor0(xi0_points, xi0_weights, n0, n1,
                     info_mode, info_vector, info_matrix, theta_use, k)
  infor0 <- i0$infor0; b <- i0$b

  scaled_of <- function(X)
    .scale_info(.build_info_data(X, info_mode, info_vector, info_matrix,
                                 theta_use)$info_data, info_mode, b)

  # ---- 3. fix one candidate set X and the approximate design on it ------
  if (use_set) {
    X <- as.matrix(candidate_set); storage.mode(X) <- "double"
    a_idx <- .nearest_idx(X, appr$support)
    cw <- .combine_idx_w(a_idx, appr$weights)
    a_idx <- cw$idx; a_w <- cw$w
  } else {
    box_lo <- meta$lo
    box_hi <- meta$hi
    # Candidate set: a coarse full-box grid (global reach for the add / exchange
    # steps) UNION a fine neighbourhood of the approximate support at the finest
    # step.  This gives fine resolution WHERE IT MATTERS without materialising
    # the full finest grid over the whole box -- which for a multi-D box can be
    # millions of points whose model evaluation dominates the run time.  Factor
    # dims always enumerate all their levels (no refinement).
    nstage <- if (is.list(step_sequence)) length(step_sequence)
              else if (is.matrix(step_sequence)) nrow(step_sequence)
              else length(step_sequence)
    if (nstage == 0L) {
      if (any(!is_factor))
        stop("'step_sequence' must contain at least one grid step for ",
             "continuous covariates.", call. = FALSE)
      X0 <- .factor_make_grid(box_lo, box_hi, 1, is_factor, nlevels)
    } else {
      # per-covariate coarsest / finest steps across stages, plus a per-covariate
      # neighbourhood radius (the second-finest step, or twice the finest when a
      # covariate has only one distinct step).  Reduces to the single-scale
      # behaviour when every stage uses one scalar step.
      stage_by <- .normalize_step_sequence(step_sequence, is_factor)
      coarse   <- do.call(pmax, stage_by)
      fine     <- do.call(pmin, stage_by)
      radius   <- vapply(seq_along(fine), function(d) {
        u <- sort(unique(vapply(stage_by, `[`, numeric(1), d)))
        if (length(u) >= 2L) u[2] else 2 * u[1]
      }, numeric(1))
      X0 <- unique(rbind(
        .factor_make_grid(box_lo, box_hi, coarse, is_factor, nlevels),
        .factor_refined_grid(box_lo, box_hi, appr$support, radius, fine,
                             is_factor, nlevels)))
    }
    if (isTRUE(snap_support)) {
      X <- X0
      a_idx <- sort(unique(.nearest_idx(X, appr$support)))
    } else {
      X <- unique(rbind(X0, as.matrix(appr$support)))
      cw <- .combine_idx_w(.nearest_idx(X, appr$support), appr$weights)
      a_idx <- cw$idx; a_w <- cw$w
    }
  }
  storage.mode(X) <- "double"
  scaled <- scaled_of(X)
  Ncand  <- ncol(scaled)

  # per-point information rank and the minimum number of support points that
  # keeps the COMBINED information matrix non-singular on the quantity of
  # interest (re-using optimal_design()'s rule).
  r_pt <- .point_info_rank(info_mode, info_vector, info_matrix, theta_use,
                           X[unique(round(seq(1, nrow(X),
                              length.out = min(7L, nrow(X))))), , drop = FALSE])
  ms <- .min_support_rule(infor0, k, wb_use, r_pt)

  # criterion / sensitivity helpers on the fixed candidate set --------------
  crit_of <- function(cnt) {
    idx <- which(cnt > 0L)
    criterion_cpp(as.integer(p), as.integer(idx),
                  as.numeric(cnt[idx] / sum(cnt)), info_mode, scaled,
                  wb_use, infor0)
  }
  # full sensitivity (directional-derivative) vector of the design `cnt`.
  dir_of <- function(cnt) {
    total <- sum(cnt); idx <- which(cnt > 0L)
    oi <- .info_ind_R(idx, cnt[idx] / total, info_mode, scaled, k)
    as.numeric(directional_deriv_cpp(as.integer(p), wb_use, info_mode, scaled,
                                     oi, infor0))
  }

  # re-optimise the reference weights on the (snapped) support so the efficiency
  # reference is the approximate optimum PROJECTED onto the candidate set.
  if (!use_set && isTRUE(snap_support)) {
    ow <- optimize_weights_cpp(as.integer(p), wb_use, info_mode,
                               scaled[, a_idx, drop = FALSE], infor0, ms)
    a_idx <- a_idx[ow$index]; a_w <- ow$weight
    crit_approx <- ow$value
  } else {
    crit_approx <- criterion_cpp(as.integer(p), as.integer(a_idx),
                                 as.numeric(a_w), info_mode, scaled,
                                 wb_use, infor0)
  }

  if (n < ms)
    stop(sprintf(paste0("n = %d is below the minimum support (%d) needed for a ",
                        "non-singular information matrix on the quantity of ",
                        "interest; increase n."), n, ms), call. = FALSE)
  if (n == ms)
    warning(sprintf(paste0("n = %d equals the minimum support (%d): the exact ",
                           "design is forced to one run per support point."),
                    n, ms), call. = FALSE)

  # ---- 4. round the approximate weights to integer counts --------------
  cnt <- integer(Ncand)
  cnt[a_idx] <- .apportion(a_w / sum(a_w), n)

  # ---- 5. adjust the total to exactly n via the sensitivity function ----
  total <- sum(cnt)
  while (total > n) {                                  # remove most inefficient
    d   <- dir_of(cnt)
    sup <- which(cnt > 0L)
    ok  <- sup[order(d[sup])]                          # smallest sensitivity first
    done <- FALSE
    for (i in ok) {
      drops_point <- cnt[i] == 1L
      if (drops_point && sum(cnt > 0L) - 1L < ms) next
      prop <- cnt; prop[i] <- prop[i] - 1L
      if (is.finite(crit_of(prop))) { cnt <- prop; done <- TRUE; break }
    }
    if (!done) stop("exact_design(): cannot reduce the design to n without ",
                    "making it singular.", call. = FALSE)
    total <- total - 1L
  }
  while (total < n) {                                  # add most informative
    d <- dir_of(cnt)
    j <- which.max(d)
    cnt[j] <- cnt[j] + 1L
    total  <- total + 1L
  }

  # ---- 6. random exchanges (accept iff strictly better) ----------------
  if (!is.null(seed)) set.seed(as.integer(seed))
  cur <- crit_of(cnt); accepted <- 0L
  if (max_exchange > 0L) for (t in seq_len(as.integer(max_exchange))) {
    sup <- which(cnt > 0L)
    i <- sup[sample.int(length(sup), 1L)]
    j <- sample.int(Ncand, 1L)
    if (j == i) next
    if (cnt[i] == 1L && sum(cnt > 0L) - 1L < ms) next  # would break min support
    prop <- cnt; prop[i] <- prop[i] - 1L; prop[j] <- prop[j] + 1L
    nc <- crit_of(prop)
    if (is.finite(nc) && nc < cur - 1e-12) {
      cnt <- prop; cur <- nc; accepted <- accepted + 1L
    }
  }

  # ---- 7. assemble the result ------------------------------------------
  idx <- sort(which(cnt > 0L))
  crit_exact <- crit_of(cnt)
  efficiency <- if (p == 0L) exp(crit_approx - crit_exact)
                else          crit_approx / crit_exact

  # counts-based combined TOTAL Fisher information (existing design + realised
  # runs): (n0 + n1) * infor0 = n0 * I_xi0 (0 if none), plus sum_i n_i * M(x_i).
  info_counts <- (n0 + n1) * infor0 +
    .opt_infor_from_support(X[idx, , drop = FALSE], cnt[idx], 1,
                            info_mode, info_vector, info_matrix, theta_use, k)
  dimnames(info_counts) <- list(coef_names, coef_names)
  # total-sample criterion: Phi_p evaluated at the TOTAL information matrix
  # (same normalised convention as `criterion`, but for all n runs at once).
  crit_total <- criterion_cpp(as.integer(p), 1L, 1.0, 1L,
                              matrix(as.numeric(info_counts), ncol = 1L),
                              wb_use, matrix(0.0, k, k))

  out <- list(
    support      = X[idx, , drop = FALSE],
    counts       = as.integer(cnt[idx]),
    weights      = cnt[idx] / n,
    n            = n,
    criterion    = crit_exact,
    criterion_total = crit_total,
    information  = info_counts,
    efficiency   = efficiency,
    criterion_approx = crit_approx,
    approx       = list(support = appr$support, weights = appr$weights,
                        criterion = appr$criterion, max_d = appr$max_d,
                        converged = appr$converged,
                        global_max_d = appr$global_max_d,
                        global_check = appr$global_check),
    p            = as.integer(p),
    candidate_set = X,
    is_factor    = is_factor,
    theta        = theta,
    coef_names   = coef_names,
    link         = model_link,
    n_exchange_accepted = accepted,
    converged    = appr$converged,
    total_time   = proc.time()[3] - t_start)
  class(out) <- "exact_design"
  out
}

#' Print an exact design.
#'
#' @param x an \code{"exact_design"} object from \code{\link{exact_design}}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.exact_design <- function(x, ...) {
  S <- as.matrix(x$support); storage.mode(S) <- "double"
  perm <- do.call(order, lapply(seq_len(ncol(S)), function(j) S[, j]))
  S <- S[perm, , drop = FALSE]; cnt <- x$counts[perm]
  cat(sprintf("\n=== Exact design (n = %d) ===\n", x$n))
  cat(sprintf("|support| = %d    efficiency >= %.2f%%    approx. max sensitivity = %.3e  (0 at the optimum)    total time = %.3f s\n",
              nrow(S), 100 * x$efficiency, x$approx$max_d, x$total_time))
  cat(sprintf("criterion (per sample) = %.6f    criterion (total sample) = %.6f\n",
              x$criterion, x$criterion_total))
  cat("support point                    count\n")
  for (i in seq_len(nrow(S))) {
    pt <- paste(.fmt_support_row(S[i, ], x$is_factor), collapse = ", ")
    cat(sprintf("  %-32s %d\n", pt, cnt[i]))
  }
  invisible(x)
}
