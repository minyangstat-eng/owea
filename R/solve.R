# ===========================================================================
# solve.R -- user-facing solvers built on the unified C++ engine.
# ===========================================================================

# Nearest row index (1-based) in X for each row of `pts` (warm-start snapping).
.nearest_idx <- function(X, pts) {
  X <- as.matrix(X); pts <- as.matrix(pts)
  vapply(seq_len(nrow(pts)), function(i) {
    d <- rowSums((X - matrix(pts[i, ], nrow(X), ncol(X), byrow = TRUE))^2)
    which.min(d)
  }, integer(1))
}

#' Run OWEA on a fixed candidate set.
#'
#' @param prob a \code{\link{DesignProblem}}.
#' @param eps0 stopping threshold on the directional derivative.
#' @param max_outer maximum number of outer iterations.
#' @param verbose print per-iteration progress.
#' @param init_weight accepted for compatibility (unused by the C++ engine).
#' @param merge if \code{TRUE}, merge neighbouring support points within
#'   \code{merge_atol} after convergence (default \code{FALSE}).
#' @param merge_atol distance threshold for merging (used only if
#'   \code{merge = TRUE}).
#' @param init_method initial-support strategy: \code{"minmax"} (default),
#'   \code{"minmaxmedian"}, \code{"random"}, \code{"iboss"}, \code{"auto"}.
#' @param full_scan_every accepted for compatibility (the engine scans the whole
#'   grid every iteration via one batched product).
#' @param warm_support,warm_weights optional warm-start design (snapped to the
#'   candidate grid).
#' @return list with \code{support}, \code{weights}, \code{criterion},
#'   \code{max_d}, \code{iterations}, \code{converged}.
#' @export
owea <- function(prob, eps0 = 1e-6, max_outer = 2000L, verbose = FALSE,
                 init_weight = 0.01, merge = FALSE, merge_atol = 1e-2,
                 init_method = "minmax", full_scan_every = "auto",
                 warm_support = NULL, warm_weights = NULL) {
  if (!inherits(prob, "DesignProblem"))
    stop("'prob' must be a DesignProblem (see DesignProblem()).", call. = FALSE)

  init_idx <- integer(0)
  if (!is.null(warm_support)) {
    init_idx <- unique(.nearest_idx(prob$X, warm_support))
    w  <- rep(1 / length(init_idx), length(init_idx))
    Ic <- prob$infor0 + .info_ind_R(init_idx, w, prob$info_mode,
                                    prob$info_data, prob$k)
    if (!(det(Ic) > 1e-20 && rcond(Ic) > 1e-14)) init_idx <- integer(0)
  }
  if (length(init_idx) == 0L)
    init_idx <- .initial_support_idx(prob$X, prob$info_mode, prob$info_data,
                                     prob$k, prob$infor0, init_method)

  res <- .solve_engine(prob$p, prob$wb, prob$info_mode, prob$info_data,
                       prob$infor0, init_idx, max_outer, eps0, verbose,
                       min_support = prob$min_support)

  support   <- prob$X[res$index, , drop = FALSE]
  weights   <- res$weight
  criterion <- res$value
  max_d     <- res$sensitivity
  converged <- res$sensitivity <= eps0

  if (isTRUE(merge)) {
    mg <- .apply_merge(support, weights, merge_atol, prob$p, prob$wb,
                       prob$info_mode, prob$info_vector, prob$info_matrix,
                       prob$theta, prob$infor0, prob$b, prob$k,
                       min_support = prob$min_support)
    if (isTRUE(mg$merged)) {
      support <- mg$support; weights <- mg$weights; criterion <- mg$value
      fb <- find_best_point(prob, support, weights)
      max_d <- fb$d; converged <- max_d <= eps0
    }
  }

  if (!converged)
    warning(sprintf("owea(): the returned design did NOT converge (max_d = %.3e > eps0 = %g); it is not optimal. Increase max_outer, coarsen the grid, or supply a warm start.",
                    max_d, eps0), call. = FALSE)

  list(support = support, weights = weights, criterion = criterion,
       max_d = max_d, iterations = res$iter, converged = converged)
}

#' Find a Phi_p-optimal design on a continuous design box.
#'
#' The recommended high-level entry point. Supply the model as EITHER
#' \code{info_vector(x, theta)} (fast path) OR \code{info_matrix(x)} (general
#' path). Candidate points come EITHER from a continuous \code{design_box} +
#' \code{step_sequence} (multistage grid refinement) OR from a fixed
#' \code{candidate_set} you provide. Supports existing designs (\code{xi0_*},
#' \code{n0}, \code{n1}) and optional neighbourhood merging.
#'
#' @param design_box list with one entry per covariate. A continuous covariate
#'   is a \code{c(lo, hi)} pair; a FACTOR (categorical) covariate is a single
#'   positive integer \code{L >= 2} giving its number of levels (the levels are
#'   the integers \code{1..L}). For example \code{list(c(2), c(3), c(0, 1))} is
#'   two factors (2 and 3 levels) and one continuous covariate on \code{[0,1]}.
#'   A factor covariate is passed to \code{info_vector} / \code{info_matrix} as
#'   its raw integer level; any contrast coding is done inside your model
#'   function. Reported support points show factor covariates as integer levels.
#'   Omit when \code{candidate_set} is supplied.
#' @param step_sequence grid steps for the multistage refinement, coarsest
#'   first; applies only to continuous covariates (factor covariates always
#'   enumerate all their levels). Either a numeric vector (each element is one
#'   stage's step, the SAME for every continuous covariate) OR, to give each
#'   covariate its own scale, a list of per-stage step vectors (or a matrix with
#'   one row per stage). Each per-stage vector may have length 1 (uniform),
#'   length equal to the number of continuous covariates, or length equal to the
#'   number of covariates (factor entries ignored). For example, with
#'   \code{design_box = list(c(0, 60), c(0, 6), c(0, 0.6))},
#'   \code{step_sequence = list(c(1, 0.1, 0.01), c(0.5, 0.05, 0.005))} runs two
#'   stages with a different step per covariate. May be omitted (or
#'   \code{numeric(0)}) when every covariate is a factor. Omit when
#'   \code{candidate_set} is supplied.
#' @param factor_levels for the \code{candidate_set} path only: an integer
#'   vector with one entry per \code{candidate_set} column marking factor
#'   columns (the number of levels \code{L >= 2}); \code{NA}/\code{0}/\code{1}
#'   marks a continuous column. Marked columns must hold integer levels in
#'   \code{1..L}; they are passed to the model as-is and printed as integers.
#'   Default \code{NULL} (all columns continuous).
#' @param candidate_set optional \eqn{n \times N} matrix of candidate design
#'   points (one row per point; need not be a regular grid). When given, the
#'   design is found on this fixed set and \code{design_box} / \code{step_sequence}
#'   are ignored.
#' @param merge_atol distance threshold for merging with \code{candidate_set}
#'   (used only if \code{merge = TRUE}); default \code{merge_factor} times the
#'   smallest positive per-coordinate gap in \code{candidate_set}.
#' @param info_vector information-vector model, either \code{function(x)} or
#'   \code{function(x, theta)} returning the length-k vector \eqn{f}.
#' @param info_matrix information-matrix model, either \code{function(x)} or
#'   \code{function(x, theta)} returning the \eqn{k \times k} matrix. Supply
#'   exactly one of \code{info_vector} / \code{info_matrix}. \code{theta} is
#'   required only for the two-argument form (or for \code{grad_g} / existing
#'   designs).
#' @param theta parameter values (required for \code{info_vector}).
#' @param link,f,x,fx,ff,xx,intercept,coding,ncat formula-style model spec, an
#'   alternative to \code{info_vector}/\code{info_matrix}. \code{link} is
#'   \code{"identity"} (linear), \code{"logit"} (logistic), \code{"loglinear"}
#'   (Poisson log-link), \code{"multinomial"} (baseline-category logit) or
#'   \code{"cumulative"} (proportional-odds ordinal logit). \code{f} = main
#'   factor effects (factor indices), \code{x} = main continuous effects
#'   (continuous indices), \code{fx} = factor \eqn{\times} continuous
#'   interactions (two-digit codes like \code{c(11, 23)} or a list of index
#'   pairs), \code{ff} = factor \eqn{\times} factor interactions (two-digit
#'   codes like \code{c(12, 13)} or index pairs; the two factors must differ),
#'   \code{xx} = continuous quadratic / interaction terms (same encoding, equal
#'   digits = quadratic). \code{intercept} (default \code{TRUE}) and
#'   \code{coding} (\code{"zero-sum"} default, or \code{"baseline"}) control the
#'   intercept and factor contrast coding. \code{ncat} is the number of response
#'   categories \eqn{J \ge 2} for \code{"multinomial"} / \code{"cumulative"}
#'   (their \code{theta} is stacked -- \eqn{(\beta_1,\dots,\beta_{J-1})} for
#'   multinomial, \eqn{(\alpha_1,\dots,\alpha_{J-1},\beta)} with increasing
#'   thresholds for cumulative). When \code{theta} is required (every link but
#'   identity) but missing, each parameter is drawn from \eqn{N(0,1)} with a
#'   warning. The result gains \code{coef_names} / \code{link}; see
#'   \code{\link{model_summary}} and \code{\link{model_info_vector}}.
#' @param p criterion: 0 = D-optimality or 1 = A-optimality (only these two are
#'   supported).
#' @param wb,subset,grad_g quantity of interest: a \eqn{v \times k} matrix
#'   \code{wb}, OR a parameter \code{subset}, OR a function \code{grad_g} (at
#'   most one); default identity (all parameters).
#' @param xi0_points,xi0_weights,n0,n1 existing-design support, weights and the
#'   existing / new sample sizes (multistage). \code{xi0_weights} are design
#'   WEIGHTS (proportions summing to 1, NOT counts); they are normalized with a
#'   warning otherwise. The existing : new balance is the ratio \code{n0 : n1}.
#'   \code{n0 = 0} is single-stage and the existing design is ignored (a warning
#'   is issued if \code{xi0_points} was supplied with \code{n0 = 0}).
#' @param merge if \code{TRUE}, merge neighbouring support points of the FINAL
#'   design once, after the multistage refinement finishes (not at every stage,
#'   which could fuse distinct optimal points at a coarse resolution). Default
#'   \code{FALSE}.
#' @param merge_factor merge tolerance multiplier: points within
#'   \code{merge_factor * step} (the final/finest step) are merged. Used only if
#'   \code{merge = TRUE}.
#' @param init_method initial-support strategy (\code{"auto"} = IBOSS for
#'   vector input, minmax for matrix input).
#' @param auto_warm_start if \code{TRUE} (default) and a \code{candidate_set}
#'   solve fails to converge from a cold start, automatically retry warm-started
#'   from a quick coarse multistage solve over the candidate set's bounding box.
#' @param check_global if \code{TRUE} (only meaningful for the \code{design_box}
#'   path), after the multistage solve converges, verify the design over a fine
#'   grid spanning the WHOLE box (the equivalence theorem) and report
#'   \code{global_max_d} / \code{global_check}. Default \code{FALSE}.
#' @param global_step grid step for the \code{check_global} verification grid;
#'   a scalar or a per-covariate vector (as in \code{step_sequence}). Default the
#'   finest per-covariate step across \code{step_sequence}. Note this
#'   materializes a full grid over the whole box, which may be large.
#' @param global_max_points safety cap (default \code{1e6}) on the
#'   \code{check_global} verification grid. If the grid would exceed it, in an
#'   interactive session you are asked whether to proceed; otherwise the global
#'   check is skipped with a warning (\code{global_check = NA}). Raise it, set a
#'   coarser \code{global_step}, or verify on a \code{candidate_set} instead.
#' @param max_iter maximum outer iterations per stage.
#' @param eps0 stopping threshold on the directional derivative.
#' @param accept_tol a refinement stage is kept only if it converges and the
#'   criterion does not worsen by more than this.
#' @param verbose print stage-by-stage progress.
#' @return list with \code{support}, \code{weights}, \code{criterion},
#'   \code{max_d}, \code{information} (the resulting per-observation information
#'   matrix \eqn{M(\xi,\theta)}; when an existing design is supplied this is the
#'   COMBINED matrix \eqn{a\,I_{\xi_0} + b\,M(\xi)} with \eqn{a = n_0/(n_0+n_1)},
#'   \eqn{b = n_1/(n_0+n_1)}), \code{converged}, \code{times}, \code{grid_sizes},
#'   \code{total_time}, \code{box_lo}, \code{box_hi}, \code{p}, and
#'   \code{global_max_d} / \code{global_check} (the whole-box equivalence-theorem
#'   maximum and whether it is \eqn{\le} \code{eps0}; \code{NA} when not checked).
#'   For the multistage (\code{design_box}) path a \code{converged} design is
#'   optimal only over the refined neighbourhood grids unless
#'   \code{check_global = TRUE} certifies it over the whole box.
#' @export
optimal_design <- function(design_box = NULL, step_sequence = NULL,
                           info_vector = NULL, info_matrix = NULL,
                           theta = NULL, p = 0L,
                           link = NULL, f = NULL, x = NULL, fx = NULL,
                           xx = NULL, ff = NULL, intercept = TRUE,
                           coding = "zero-sum", ncat = NULL,
                           wb = NULL, subset = NULL, grad_g = NULL,
                           xi0_points = NULL, xi0_weights = numeric(0),
                           n0 = 0, n1 = 1, candidate_set = NULL,
                           factor_levels = NULL,
                           merge = FALSE, merge_factor = 1.5, merge_atol = NULL,
                           init_method = "auto", auto_warm_start = TRUE,
                           check_global = FALSE, global_step = NULL,
                           global_max_points = 1e6,
                           max_iter = 100L, eps0 = 1e-6,
                           accept_tol = 1e-9, verbose = FALSE) {
  p <- .check_criterion(p)

  # a formula-style model spec ('link' + f/x/fx/ff/xx) is a third way to specify
  # the model; it builds info_vector (and draws theta ~ N(0,1) for logit/loglinear
  # when theta is missing).
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
  # step_sequence may be omitted for an all-factor design_box (no continuous
  # covariates to step over); the continuous-covariate case is validated below.
  if (!use_set && is.null(step_sequence)) step_sequence <- numeric(0)
  if (info_mode == 0L) {
    if (.info_vec_needs_theta(info_vector) && is.null(theta))
      stop("'theta' is required when 'info_vector' is function(x, theta).",
           call. = FALSE)
    info_vector <- .normalize_info_vector(info_vector, theta)
  }
  # accept info_matrix(x) or info_matrix(x, theta)
  info_matrix <- .normalize_info_matrix(info_matrix, theta)

  # ---- factor (categorical) covariate metadata ----------------------------
  # Factor covariates are carried in LEVEL space (integers 1..L); the model
  # receives the raw level and does any contrast coding itself.
  meta <- if (use_set)
            .factor_levels_to_meta(factor_levels, ncol(as.matrix(candidate_set)))
          else .parse_design_box(design_box)
  is_factor <- meta$is_factor; nlevels <- meta$nlevels
  if (use_set)
    .validate_factor_columns(candidate_set, is_factor, nlevels, "candidate_set")
  if (!is.null(xi0_points))
    .validate_factor_columns(xi0_points, is_factor, nlevels, "xi0_points")
  if (any(is_factor) && isTRUE(merge)) {
    warning("optimal_design(): merging is not supported with factor covariates ",
            "(centroid merging would create invalid factor levels); proceeding ",
            "with merge = FALSE.", call. = FALSE)
    merge <- FALSE
  }

  # infer k from a probe point (a candidate point, or the box centre); factor
  # dims use a valid level (1) rather than the non-integer box centre.
  probe <- if (use_set) as.numeric(as.matrix(candidate_set)[1, ])
           else .factor_probe_point(meta$lo, meta$hi, is_factor)
  k <- if (info_mode == 1L) nrow(as.matrix(info_matrix(probe)))
       else length(as.numeric(info_vector(probe, theta)))
  theta_use <- if (is.null(theta)) rep(0.0, k) else as.numeric(theta)

  wb_use <- .wb_from(wb, subset, grad_g, theta_use, k)
  i0 <- .make_infor0(xi0_points, xi0_weights, n0, n1,
                     info_mode, info_vector, info_matrix, theta_use, k)
  infor0 <- i0$infor0; b <- i0$b

  # rank-aware minimum support: ceil((rank(wb) - existing coverage) / per-point
  # information rank).  Evaluate the per-point rank at representative points.
  samp <- if (use_set) {
            Xs <- as.matrix(candidate_set)
            Xs[unique(round(seq(1, nrow(Xs), length.out = min(7L, nrow(Xs))))), ,
               drop = FALSE]
          } else {
            bl <- meta$lo; bh <- meta$hi
            t(vapply(c(0.5, 0.25, 0.75, 0.1, 0.9),
                     function(a) bl + a * (bh - bl), numeric(length(bl))))
          }
  if (any(is_factor)) samp <- .snap_factor_levels(samp, is_factor, nlevels)
  r_pt <- .point_info_rank(info_mode, info_vector, info_matrix, theta_use, samp)
  ms   <- .min_support_rule(infor0, k, wb_use, r_pt)

  # scaled candidate information for a grid X (built once per grid)
  scaled_of <- function(X)
    .scale_info(.build_info_data(X, info_mode, info_vector, info_matrix,
                                 theta_use)$info_data, info_mode, b)

  # solve on a grid whose scaled information is already built (no merging here;
  # merging is applied once at the end via finalize_merge()).
  solve_prepared <- function(X, scaled, init_idx) {
    res <- .solve_engine(p, wb_use, info_mode, scaled, infor0,
                         init_idx, max_iter, eps0, FALSE, min_support = ms)
    list(support = X[res$index, , drop = FALSE], weights = res$weight,
         criterion = res$value, max_d = res$sensitivity,
         converged = res$sensitivity <= eps0, iterations = res$iter)
  }
  solve_on_grid <- function(X, init_idx) solve_prepared(X, scaled_of(X), init_idx)

  # auto-warm-start retry: if a cold solve does not converge, warm-start it from
  # a coarse multistage solve over [blo, bhi] (makes a hard first stage /
  # candidate set converge instead of getting stuck on a rank-deficient design).
  robust_solve <- function(X, scaled, blo, bhi) {
    init <- .initial_support_idx(X, info_mode, scaled, k, infor0, init_method)
    cand <- solve_prepared(X, scaled, init)
    if (!cand$converged && isTRUE(auto_warm_start)) {
      if (verbose)
        cat("  solve did not converge; warm-starting from a coarse multistage solve...\n")
      sc   <- max(bhi - blo) / 6
      warm <- multistage(blo, bhi,
                         .normalize_step_sequence(c(sc, sc / 2, sc / 4), is_factor),
                         robust_first = FALSE)
      init2 <- unique(.nearest_idx(X, warm$support))
      w2 <- rep(1 / length(init2), length(init2))
      Ic <- infor0 + .info_ind_R(init2, w2, info_mode, scaled, k)
      if (det(Ic) > 1e-20 && rcond(Ic) > 1e-14) {
        cand2 <- solve_prepared(X, scaled, init2)
        if (cand2$converged || cand2$criterion < cand$criterion) cand <- cand2
      }
    }
    cand
  }

  # merge neighbouring support points ONCE (at the final stage / candidate set),
  # re-optimise the weights, and re-verify max_d on the given scaled grid.
  finalize_merge <- function(cand, scaled_grid, atol) {
    if (!isTRUE(merge)) return(cand)
    mg <- .apply_merge(cand$support, cand$weights, atol, p, wb_use, info_mode,
                       info_vector, info_matrix, theta_use, infor0, b, k,
                       min_support = ms)
    if (!isTRUE(mg$merged)) return(cand)
    oi <- .opt_infor_from_support(mg$support, mg$weights, b, info_mode,
                                  info_vector, info_matrix, theta_use, k)
    ve <- verify_equiv_cpp(as.integer(p), wb_use, info_mode, scaled_grid, oi, infor0)
    cand$support <- mg$support; cand$weights <- mg$weights
    cand$criterion <- mg$value; cand$max_d <- ve$max_d
    cand$converged <- ve$max_d <= eps0
    cand
  }

  # multistage refinement over a box (NO merging in the stages).  Tracks the grid
  # of the last accepted stage so merging can be applied once, at the end.
  # `steps` is a list of per-covariate `by` vectors, one per stage (coarsest
  # first); see .normalize_step_sequence().
  multistage <- function(blo, bhi, steps, robust_first = FALSE) {
    res <- NULL; res_X <- NULL; res_step <- NULL
    times <- numeric(0); gsz <- integer(0)
    fmt_step <- function(s) paste(format(s[!is_factor]), collapse = ",")
    for (i in seq_along(steps)) {
      step <- steps[[i]]
      if (i == 1L) {
        X  <- .factor_make_grid(blo, bhi, step, is_factor, nlevels)
        sX <- scaled_of(X)
        tt <- system.time(
          cand <- if (robust_first) robust_solve(X, sX, blo, bhi)
                  else solve_prepared(X, sX,
                         .initial_support_idx(X, info_mode, sX, k, infor0, init_method)))[3]
        res <- cand; res_X <- X; res_step <- step
        times <- c(times, tt); gsz <- c(gsz, nrow(X))
        if (verbose)
          cat(sprintf("  step=(%s)  |X|=%-10d time=%8.3f s  |S|=%-2d crit=%.6f conv=%s\n",
                      fmt_step(step), nrow(X), tt, nrow(res$support), res$criterion, res$converged))
      } else {
        Xpts <- .factor_refined_grid(blo, bhi, res$support, steps[[i - 1L]], step,
                                     is_factor, nlevels)
        init <- unique(.nearest_idx(Xpts, res$support))
        tt   <- system.time(cand <- solve_on_grid(Xpts, init))[3]
        times <- c(times, tt); gsz <- c(gsz, nrow(Xpts))
        if (cand$converged && cand$criterion <= res$criterion + accept_tol) {
          res <- cand; res_X <- Xpts; res_step <- step
          if (verbose)
            cat(sprintf("  step=(%s)  |X|=%-10d time=%8.3f s  |S|=%-2d crit=%.6f\n",
                        fmt_step(step), nrow(Xpts), tt, nrow(res$support), res$criterion))
        } else if (verbose) {
          cat(sprintf("  step=(%s)  |X|=%-10d time=%8.3f s  -> rejected (crit=%.6f, converged=%s); keeping previous\n",
                      fmt_step(step), nrow(Xpts), tt, cand$criterion, cand$converged))
        }
      }
    }
    res$times <- times; res$grid_sizes <- gsz
    res$final_X <- res_X; res$final_step <- res_step
    res
  }

  # max directional derivative of a design over a fine grid spanning the box
  global_md <- function(support, weights, blo, bhi, step) {
    Xv <- .factor_make_grid(blo, bhi, step, is_factor, nlevels)
    oi <- .opt_infor_from_support(support, weights, b, info_mode,
                                  info_vector, info_matrix, theta_use, k)
    ve <- verify_equiv_cpp(as.integer(p), wb_use, info_mode, scaled_of(Xv),
                           oi, infor0)
    list(max_d = ve$max_d, npoints = nrow(Xv))
  }

  # ---- fixed candidate set: a single solve (+ auto warm-start fallback) ----
  if (use_set) {
    X <- as.matrix(candidate_set); storage.mode(X) <- "double"
    atol <- if (isTRUE(merge)) {
      if (!is.null(merge_atol)) merge_atol else {
        gaps <- apply(X, 2, function(col) {
          u <- sort(unique(col)); d <- diff(u); d <- d[d > 0]
          if (length(d)) min(d) else NA_real_
        })
        gaps <- gaps[is.finite(gaps)]
        merge_factor * (if (length(gaps)) min(gaps) else 1e-2)
      }
    } else NA_real_
    scaledX <- scaled_of(X)
    ttot <- system.time({
      cand <- robust_solve(X, scaledX, apply(X, 2, min), apply(X, 2, max))
      cand <- finalize_merge(cand, scaledX, atol)
    })[3]
    if (verbose)
      cat(sprintf("  candidate set |X|=%-10d time=%8.3f s  |S|=%-2d crit=%.6f converged=%s\n",
                  nrow(X), ttot, nrow(cand$support), cand$criterion, cand$converged))

    info_out <- infor0 + .opt_infor_from_support(cand$support, cand$weights, b,
                                                 info_mode, info_vector,
                                                 info_matrix, theta_use, k)
    dimnames(info_out) <- list(coef_names, coef_names)
    out <- list(support = cand$support, weights = cand$weights,
                criterion = cand$criterion, max_d = cand$max_d,
                information = info_out,
                converged = cand$converged, times = ttot,
                grid_sizes = nrow(X), total_time = ttot,
                box_lo = apply(X, 2, min), box_hi = apply(X, 2, max), p = p,
                is_factor = is_factor, theta = theta,
                coef_names = coef_names, link = model_link,
                global_max_d = if (cand$converged) cand$max_d else NA_real_,
                global_check = cand$converged)
    if (!cand$converged)
      warning(sprintf("optimal_design(): the returned design did NOT converge (max_d = %.3e > eps0 = %g); it is not optimal. Increase max_iter, coarsen the grid, or supply a warm start.",
                      cand$max_d, eps0), call. = FALSE)
    return(out)
  }

  # ---- design box (continuous and/or factor): multistage grid refinement ----
  box_lo <- meta$lo
  box_hi <- meta$hi
  # step_sequence applies only to continuous covariates.  Normalize it into a
  # list of per-covariate `by` vectors (one per stage); a scalar per stage is
  # the classic uniform-step form, while a list/matrix lets each covariate use
  # its own scale.  An all-factor design space is fully enumerated in a single
  # stage (refinement does nothing to categorical dims), so collapse it.
  nstage <- if (is.list(step_sequence)) length(step_sequence)
            else if (is.matrix(step_sequence)) nrow(step_sequence)
            else length(step_sequence)
  if (nstage == 0L) {
    if (any(!is_factor))
      stop("'step_sequence' must contain at least one grid step for continuous ",
           "covariates.", call. = FALSE)
    stage_by <- list(rep(1, length(is_factor)))  # all-factor: one stage, unused
  } else {
    stage_by <- .normalize_step_sequence(step_sequence, is_factor)
    if (all(is_factor)) stage_by <- stage_by[1]  # factors: one stage suffices
  }
  res <- multistage(box_lo, box_hi, stage_by, robust_first = auto_warm_start)
  # merge neighbouring support points ONLY at the end, on the final (finest)
  # grid -- not at every stage (which could fuse distinct optimal points early).
  # With per-covariate steps the merge tolerance is a single Euclidean radius,
  # so use the smallest continuous step of the final stage.
  if (isTRUE(merge))
    res <- finalize_merge(res, scaled_of(res$final_X),
                          merge_factor * min(res$final_step[!is_factor]))

  info_out <- infor0 + .opt_infor_from_support(res$support, res$weights, b,
                                               info_mode, info_vector,
                                               info_matrix, theta_use, k)
  dimnames(info_out) <- list(coef_names, coef_names)
  out <- list(support = res$support, weights = res$weights,
              criterion = res$criterion, max_d = res$max_d,
              information = info_out,
              converged = res$converged, times = res$times,
              grid_sizes = res$grid_sizes, total_time = sum(res$times),
              box_lo = box_lo, box_hi = box_hi, p = p,
              is_factor = is_factor, theta = theta,
              coef_names = coef_names, link = model_link,
              global_max_d = NA_real_, global_check = NA)

  if (isTRUE(res$converged)) {
    if (isTRUE(check_global)) {
      # default verification step: the finest per-covariate step across stages.
      gstep <- if (is.null(global_step)) do.call(pmin, stage_by)
               else .expand_stage_step(as.numeric(global_step), is_factor,
                                       sum(!is_factor), length(is_factor))
      gstep_lab <- paste(format(gstep[!is_factor]), collapse = ",")
      npts  <- prod(ifelse(is_factor, nlevels,               # grid size, no build
                           round((box_hi - box_lo) / gstep) + 1))
      do_it <- TRUE
      if (npts > global_max_points) {
        msg <- sprintf("optimal_design(): the whole-box global check at step (%s) would evaluate %.0f design points (> global_max_points = %.0f).",
                       gstep_lab, npts, global_max_points)
        do_it <- if (interactive())
                   isTRUE(utils::askYesNo(paste(msg, "Evaluate all of them anyway?"),
                                          default = FALSE))
                 else FALSE
        if (!do_it)
          warning(paste(msg, "Global check SKIPPED (global_check = NA). Raise global_max_points, set a coarser global_step, or verify on a fixed candidate_set."),
                  call. = FALSE)
      }
      if (do_it) {
        g <- global_md(res$support, res$weights, box_lo, box_hi, gstep)
        out$global_max_d <- g$max_d
        out$global_check <- (g$max_d <= eps0)
        if (verbose)
          cat(sprintf("  global check (step (%s), |X|=%d): max_d = %.3e -> %s\n",
                      gstep_lab, g$npoints, g$max_d,
                      if (out$global_check) "GLOBAL optimum" else "LOCAL only"))
        if (!out$global_check)
          warning(sprintf("optimal_design(): the design is optimal over the refined grids but NOT over the whole box (global max_d = %.3e at step (%s) > eps0 = %g) -- it is only LOCALLY optimal. Try a finer/longer step_sequence, a smaller global_step, or solve on a fixed candidate_set warm-started from this design.",
                          g$max_d, gstep_lab, eps0), call. = FALSE)
      }
    } else {
      warning("optimal_design(): multistage convergence is LOCAL -- the equivalence theorem was verified only over the refined neighbourhood grids, not the whole design box. Pass check_global = TRUE to verify global optimality.",
              call. = FALSE)
    }
  } else {
    warning(sprintf("optimal_design(): the returned design did NOT converge (max_d = %.3e > eps0 = %g); it is not optimal. Increase max_iter or adjust the step_sequence.",
                    res$max_d, eps0), call. = FALSE)
  }
  out
}

#' Solve for an optimal design from an information-vector model (fast path).
#'
#' Compatibility entry point of the information-vector package. Candidate points
#' may be given as a precomputed \eqn{k \times N} matrix \code{infor_vec_all},
#' or as a design matrix \code{X} with a per-point function \code{infor_vec} and
#' \code{theta}.
#'
#' @param pp criterion: 0 = D-optimality or 1 = A-optimality (only these two are
#'   supported).
#' @param wb \eqn{v \times k} selection matrix \eqn{dg/d\theta}.
#' @param infor_vec_all \eqn{k \times N} information-vector matrix.
#' @param X \eqn{N \times d} design matrix (one row per candidate point).
#' @param infor_vec model function \code{infor_vec(x, theta)}.
#' @param theta parameter vector.
#' @param max_iter maximum outer iterations.
#' @param tol convergence tolerance on the directional derivative.
#' @param verbose print per-iteration progress.
#' @return list with \code{index}, \code{weight}, \code{points} (if \code{X}
#'   given), \code{sensitivity}, \code{iter}, \code{value}, \code{design}.
#' @export
appro_opt <- function(pp, wb, infor_vec_all = NULL,
                      X = NULL, infor_vec = NULL, theta = NULL,
                      max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  pp <- .check_criterion(pp)
  if (is.null(infor_vec_all)) {
    if (is.null(X) || is.null(infor_vec) || is.null(theta))
      stop("Supply either 'infor_vec_all', or all of 'X', 'infor_vec' and ",
           "'theta'.", call. = FALSE)
    infor_vec_all <- build_infor_vec_all(X, infor_vec, theta)
  }
  infor_vec_all <- as.matrix(infor_vec_all); storage.mode(infor_vec_all) <- "double"
  wb <- as.matrix(wb); storage.mode(wb) <- "double"
  k  <- nrow(infor_vec_all)
  if (ncol(wb) != k) stop("wb must have ncol = nrow(infor_vec_all) = ", k, ".")
  if (k < 2) stop("need at least 2 parameters (nrow(infor_vec_all) >= 2).")

  res <- .solve_engine(pp, wb, 0L, infor_vec_all, matrix(0.0, k, k),
                       integer(0), max_iter, tol, verbose)
  res$design <- cbind(index = res$index, weight = res$weight)
  if (!is.null(X)) res$points <- as.matrix(X)[res$index, , drop = FALSE]
  res
}

#' Sequential (multi-resolution) optimal design from an information-vector model.
#'
#' @param pp,wb as in \code{\link{appro_opt}}.
#' @param lower,upper range of each covariate.
#' @param by_seq list of grid steps, coarsest first.
#' @param infor_vec model function \code{infor_vec(x, theta)}.
#' @param theta parameter vector.
#' @param max_iter,tol,verbose as in \code{\link{appro_opt}}.
#' @return list with \code{points}, \code{weight}, \code{design},
#'   \code{sensitivity}, \code{iter}, \code{value}.
#' @export
appro_opt_seq <- function(pp, wb, lower, upper, by_seq, infor_vec, theta,
                          max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  d <- length(lower)
  if (length(upper) != d)
    stop("'lower' and 'upper' must have the same length", call. = FALSE)
  if (!is.list(by_seq)) by_seq <- as.list(by_seq)
  by_seq <- lapply(by_seq, function(b) if (length(b) == 1L) rep(b, d) else b)
  nst <- length(by_seq)
  if (nst < 1L) stop("'by_seq' must contain at least one grid size", call. = FALSE)

  X   <- make_grid(lower, upper, by_seq[[1]])
  if (verbose)
    cat(sprintf("stage 1: step (%s) -> %d grid points\n",
                paste(format(by_seq[[1]]), collapse = ", "), nrow(X)))
  res <- appro_opt(pp, wb, X = X, infor_vec = infor_vec, theta = theta,
                   max_iter = max_iter, tol = tol)
  pts <- X[res$index, , drop = FALSE]
  if (verbose)
    cat(sprintf("         -> %d support points, criterion %.6f\n",
                nrow(pts), res$value))

  for (s in seq_len(nst)[-1]) {
    X   <- .refine_grid(pts, by_seq[[s - 1]], by_seq[[s]], lower, upper)
    if (verbose)
      cat(sprintf("stage %d: step (%s) -> %d grid points\n",
                  s, paste(format(by_seq[[s]]), collapse = ", "), nrow(X)))
    res <- appro_opt(pp, wb, X = X, infor_vec = infor_vec, theta = theta,
                     max_iter = max_iter, tol = tol)
    pts <- X[res$index, , drop = FALSE]
    if (verbose)
      cat(sprintf("         -> %d support points, criterion %.6f\n",
                  nrow(pts), res$value))
  }

  list(points = pts, weight = res$weight,
       design = cbind(pts, weight = res$weight),
       sensitivity = res$sensitivity, iter = res$iter, value = res$value)
}
