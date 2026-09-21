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
#'   \code{"minmaxmedian"}, \code{"random"}, \code{"iboss"}, \code{"MA"},
#'   \code{"auto"}. \code{"MA"} runs a multiplicative algorithm (equal-weight
#'   start, D-optimality) and keeps the \code{k + 1} highest-weight candidate
#'   points as the starting support -- often a much better warm start in higher
#'   dimensions.
#' @param full_scan_every accepted for compatibility (the engine scans the whole
#'   grid every iteration via one batched product).
#' @param warm_support,warm_weights optional warm-start design (snapped to the
#'   candidate grid).
#' @param ma_max_iter maximum number of multiplicative-algorithm iterations used
#'   by \code{init_method = "MA"} (default 100); ignored by the other methods.
#' @return list with \code{support}, \code{weights}, \code{criterion},
#'   \code{max_d}, \code{iterations}, \code{converged}.
#' @export
owea <- function(prob, eps0 = 1e-6, max_outer = 2000L, verbose = FALSE,
                 init_weight = 0.01, merge = FALSE, merge_atol = 1e-2,
                 init_method = "minmax", full_scan_every = "auto",
                 warm_support = NULL, warm_weights = NULL,
                 ma_max_iter = 100L) {
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
                                     prob$k, prob$infor0, init_method,
                                     ma_max_iter = ma_max_iter)

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
#'   vector input, minmax for matrix input). Other choices: \code{"minmax"},
#'   \code{"minmaxmedian"}, \code{"random"}, \code{"iboss"}, and \code{"MA"} (a
#'   multiplicative algorithm that keeps the \code{k + 1} highest-weight points as
#'   the warm-start support -- often faster in higher dimensions).
#' @param solver design algorithm: \code{"owea"} (default) runs the OWEA exchange
#'   engine (initialised per \code{init_method}); \code{"MA"} (alias
#'   \code{"multiplicative"}) instead runs the general multiplicative algorithm
#'   (Yu 2010) as a \emph{direct} solver and returns its design, skipping the
#'   exchange engine -- typically far faster for larger problems. \code{"MA"}
#'   handles both D-optimality (\code{p = 0}) and A-optimality (\code{p = 1}), any
#'   quantity of interest (\code{subset} / \code{grad_g} / \code{wb}) and an
#'   existing design (\code{n0 > 0}). For a criterion outside its
#'   guaranteed-convergent class (e.g. c-optimality / rank-1 \code{wb}), any
#'   per-grid solve that fails to certify optimality automatically falls back to
#'   the OWEA engine, so the result is always correct. When \code{solver = "MA"},
#'   \code{init_method} and \code{auto_warm_start} are unused. \code{"owea"} is
#'   generally faster, but \code{"MA"} has the advantage for a large dimension of
#'   the information matrix (many parameters / a big grid).
#' @param ma_max_iter maximum number of multiplicative-algorithm iterations, used
#'   both by \code{init_method = "MA"} and \code{solver = "MA"} (default 100; the
#'   direct solver raises the effective cap to at least 10000 so it can converge).
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
#' @param engine \code{"classic"} (default) runs the original OWEA engine.
#'   \code{"active-set"} replaces the weight step by an active-set Newton method:
#'   the step is cut at the first weight that would reach zero and accepted
#'   under an Armijo decrease of the criterion, all points the step pushes to
#'   zero are pruned in one batch (never a point whose directional derivative is
#'   still positive, and never below the rank-aware minimum support or past the
#'   estimability of the quantity of interest), and \code{add_per_iter}
#'   candidates are added per exchange iteration. Applies to every criterion
#'   and \code{wb}/\code{subset}/\code{grad_g}; both engines converge to the
#'   same design, the active-set engine typically in far fewer Newton solves.
#' @param add_per_iter number of violating candidate points added per exchange
#'   iteration when \code{engine = "active-set"} (at most 20\% of the current
#'   support per iteration); ignored by the classic engine, which adds one.
#' @param continuous if \code{TRUE} (\code{design_box} path only), search the
#'   continuous covariates in the CONTINUOUS design region instead of on a grid;
#'   no \code{step_sequence} is needed. The solver starts from the OWEA solution
#'   on a small coarse grid (\code{init_levels} levels per continuous covariate
#'   times all levels of every factor), then alternates (i) polishing -- with
#'   the weights fixed, the criterion is minimised over the continuous
#'   coordinates of all support points by L-BFGS-B, using the closed-form
#'   gradient \eqn{-(w_i/v)\,\partial\,\mathrm{tr}(P\,I(x))/\partial x} at
#'   \eqn{x_i}, where \eqn{P} is the matrix behind the directional derivative
#'   and \eqn{v} the number of parameters of interest -- (ii) the package's
#'   Newton weight step on the current support, (iii) merging of support points
#'   closer than \code{merge_tol}, and (iv) adding the maximiser of the
#'   directional derivative over the region, found by multi-start L-BFGS-B
#'   (started from the support points, the box vertices and \code{n_starts}
#'   random points, for every level combination of the factors). It stops when
#'   that maximum is \eqn{\le} \code{eps0} and a random audit of \code{n_audit}
#'   points (its best points polished locally) finds no violation. The
#'   certificate (\code{max_d}, \code{efficiency_lower_bound}) is the largest
#'   directional derivative over the points examined, not an exhaustive scan.
#'   Factor covariates keep their levels. On this path \code{step_sequence} and
#'   \code{merge} are ignored and \code{check_global} needs \code{global_step}.
#'   Default \code{FALSE}: the grid path is unchanged.
#' @param info_jacobian optional derivative of the per-point information with
#'   respect to the continuous covariates, used by \code{continuous = TRUE}:
#'   \code{function(x, theta)} (or \code{function(x)}) returning a
#'   \eqn{k \times n_c} matrix for \code{info_vector} (\eqn{\partial f/\partial
#'   x}) or a \eqn{k^2 \times n_c} matrix for \code{info_matrix} (the
#'   column-major \eqn{\partial\,\mathrm{vec}\,I(x)/\partial x}), with
#'   \eqn{n_c} the number of continuous covariates in \code{design_box} order.
#'   Default \code{NULL}: a formula-style model uses its analytic derivative
#'   and a user-supplied model function is differentiated by finite differences
#'   (central inside the box, one-sided at its boundary).
#' @param init_levels number of equally spaced levels per continuous covariate
#'   of the coarse grid the continuous search starts from (default 3: the two
#'   ends and the midpoint; raised automatically if the grid would have fewer
#'   than \eqn{k + 1} points).
#' @param init_points optional matrix of starting support points for the
#'   continuous search (one row per point, one column per covariate, factor
#'   columns holding levels), replacing the coarse grid.
#' @param n_starts number of random starts of the multi-start maximisation of
#'   the directional derivative (default 30, shared across the level
#'   combinations of the factors), on top of the support points and the box
#'   vertices.
#' @param n_audit number of random points of the final audit (default 20000;
#'   0 skips it). It costs one model evaluation per point: about a second per
#'   million points for a formula-style model, one to two seconds per 100,000
#'   points for a user-supplied model function.
#' @param merge_tol continuous search only: two support points with identical
#'   factor levels whose continuous coordinates, each scaled by the range of
#'   its covariate, lie within this Euclidean distance are merged (default
#'   \code{1e-3}).
#' @param verbose print stage-by-stage progress.
#' @return list with \code{support}, \code{weights}, \code{criterion},
#'   \code{max_d}, \code{information} (the resulting per-observation information
#'   matrix \eqn{M(\xi,\theta)}; when an existing design is supplied this is the
#'   COMBINED matrix \eqn{a\,I_{\xi_0} + b\,M(\xi)} with \eqn{a = n_0/(n_0+n_1)},
#'   \eqn{b = n_1/(n_0+n_1)}), \code{converged}, \code{times} (elapsed seconds
#'   in the solver itself, one entry per stage --- including any stage whose
#'   refinement was computed and then rejected), \code{grid_sizes},
#'   \code{total_time} (elapsed seconds for the WHOLE call, as
#'   \code{\link{exact_design}} and \code{compound_design} also report it: grid
#'   construction, the model evaluation over it and any \code{check_global}
#'   verification included, so it exceeds \code{sum(times)} --- several-fold on
#'   a large grid, where building the candidate set dominates the solve),
#'   \code{box_lo}, \code{box_hi}, \code{p}, \code{efficiency_lower_bound} (a
#'   guaranteed lower bound on the design's efficiency relative to the TRUE
#'   optimum over the design space it was checked on, computed from the design
#'   alone: \eqn{1/(1 + \mathrm{max\_d}/v)} for D and
#'   \eqn{1 - \mathrm{max\_d}/\mathrm{criterion}} for A, with \eqn{v} the number
#'   of parameters of interest --- Theorems 4.5 and 4.6 of Becker and Yang,
#'   \emph{Post Hoc Control Group Selection via Constrained Optimal Design};
#'   it equals 1 at a certified optimum and stays valid when \code{converged}
#'   is \code{FALSE}), and
#'   \code{global_max_d} / \code{global_check} (the whole-box equivalence-theorem
#'   maximum and whether it is \eqn{\le} \code{eps0}; \code{NA} when not checked;
#'   when the whole-box check runs, \code{efficiency_lower_bound} is recomputed
#'   from \code{global_max_d}).
#'   For the multistage (\code{design_box}) path a \code{converged} design is
#'   optimal only over the refined neighbourhood grids unless
#'   \code{check_global = TRUE} certifies it over the whole box.
#'   With \code{continuous = TRUE} the result also carries
#'   \code{method = "continuous"}, \code{iterations}, \code{history} (one row
#'   per iteration: time, support size, criterion, max sensitivity),
#'   \code{maximiser} (the point attaining \code{max_d}), \code{n_audit},
#'   \code{audit_max_d} and \code{jacobian} (\code{"analytic"},
#'   \code{"finite differences"} or \code{"user"}); \code{global_max_d} /
#'   \code{global_check} then repeat \code{max_d} / \code{converged}, since the
#'   region-wide search is the check, \code{grid_sizes} is the size of the
#'   coarse starting grid and \code{times} holds its solve time and the time of
#'   the continuous search.
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
                           solver = "owea", ma_max_iter = 100L,
                           check_global = FALSE, global_step = NULL,
                           global_max_points = 1e6,
                           max_iter = 100L, eps0 = 1e-6,
                           engine = c("classic", "active-set"), add_per_iter = 1L,
                           accept_tol = 1e-9,
                           continuous = FALSE, info_jacobian = NULL,
                           init_levels = 3L, init_points = NULL,
                           n_starts = 30L, n_audit = 20000L, merge_tol = 1e-3,
                           verbose = FALSE) {
  # Wall clock for the WHOLE call, as exact_design() and compound_design() also
  # report it: the grid construction, the model evaluation over it and any
  # check_global verification are part of what the call costs, and timing only
  # the solver understated a large-grid run several-fold.  The per-stage solver
  # times stay in `times`.
  t_start <- proc.time()[3]
  p <- .check_criterion(p)
  engine <- match.arg(engine)
  engine_mode <- if (identical(engine, "active-set")) 1L else 0L
  add_per_iter <- suppressWarnings(as.integer(add_per_iter))[1]
  if (is.na(add_per_iter) || add_per_iter < 1L)
    stop("'add_per_iter' must be a positive integer.", call. = FALSE)

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

  # c-optimality -- ONE linear combination c'theta (a one-row wb / grad_g, or a
  # one-parameter subset).  Its optimal design is often singular, which the
  # sensitivity function cannot certify; Elfving's bound can (see elfving.R).
  # The design's value V = c' M^- c is exp(criterion) for p = 0, criterion for
  # p = 1 (both normalised by v = 1).
  is_c    <- nrow(wb_use) == 1L
  cvec    <- if (is_c) as.numeric(wb_use) else NULL
  c_value <- function(crit) if (p == 0L) exp(crit) else crit
  certify_c <- function(cand, scaled) {              # grid / candidate-set paths
    if (!is_c || is.null(cand$weights) || !length(cand$weights)) return(cand)
    M <- infor0 + .opt_infor_from_support(cand$support, cand$weights, b, info_mode,
                                          info_vector, info_matrix, theta_use, k)
    cert <- .c_certificate(cvec, info_mode, scaled, infor0, k,
                           c_value(cand$criterion), eps0, M = M)
    cand$elfving   <- cert
    cand$converged <- isTRUE(cand$converged) || isTRUE(cert$certified)
    cand
  }

  # solver: "owea" (OWEA exchange engine, default) or "MA" (the general
  # multiplicative algorithm as a direct solver).  MA covers D- and A-optimality,
  # any quantity of interest (subset / grad_g / wb) and an existing design; for a
  # criterion outside its guaranteed-convergent class (e.g. c-optimality / rank-1
  # wb) each per-grid MA solve that fails to certify falls back to the engine.
  solver <- tolower(as.character(solver)[1])
  use_ma <- solver %in% c("ma", "multiplicative")
  if (!use_ma && !identical(solver, "owea"))
    stop("solver must be \"owea\" or \"MA\".", call. = FALSE)

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
                         init_idx, max_iter, eps0, FALSE, min_support = ms,
                         engine_mode = engine_mode, add_per_iter = add_per_iter)
    certify_c(list(support = X[res$index, , drop = FALSE], weights = res$weight,
                   criterion = res$value, max_d = res$sensitivity,
                   converged = res$sensitivity <= eps0, iterations = res$iter),
              scaled)
  }
  solve_on_grid <- function(X, init_idx) solve_prepared(X, scaled_of(X), init_idx)

  # solver = "MA": run the general multiplicative algorithm to convergence and
  # return its design directly (no OWEA engine).  Handles D-/A-optimality, any
  # quantity of interest (wb_use) and an existing design (infor0).
  ma_cap <- max(as.integer(ma_max_iter), 10000L)
  ma_solve_on_grid <- function(X, scaled)
    .ma_solve(X, scaled, info_mode, k, wb_use, p, infor0, eps0, ma_cap)
  # try MA on a grid; if it is singular or fails to certify (e.g. c-optimality),
  # fall back to the OWEA engine so the result is always correct.
  ma_or <- function(X, scaled, owea_thunk) {
    r <- ma_solve_on_grid(X, scaled)
    if (!isTRUE(r$singular) && isTRUE(r$converged)) return(r)
    if (verbose) cat("  solver=MA did not certify on this grid; using OWEA.\n")
    owea_thunk()
  }

  # auto-warm-start retry: if a cold solve does not converge, warm-start it from
  # a coarse multistage solve over [blo, bhi] (makes a hard first stage /
  # candidate set converge instead of getting stuck on a rank-deficient design).
  robust_solve <- function(X, scaled, blo, bhi) {
    init <- .initial_support_idx(X, info_mode, scaled, k, infor0, init_method,
                                 ma_max_iter = ma_max_iter)
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
    certify_c(cand, scaled_grid)
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
          cand <- if (use_ma)
                    ma_or(X, sX, function()
                      if (robust_first) robust_solve(X, sX, blo, bhi)
                      else solve_prepared(X, sX,
                             .initial_support_idx(X, info_mode, sX, k, infor0,
                                                  init_method, ma_max_iter = ma_max_iter)))
                  else if (robust_first) robust_solve(X, sX, blo, bhi)
                  else solve_prepared(X, sX,
                         .initial_support_idx(X, info_mode, sX, k, infor0, init_method,
                                              ma_max_iter = ma_max_iter)))[3]
        res <- cand; res_X <- X; res_step <- step
        times <- c(times, tt); gsz <- c(gsz, nrow(X))
        if (verbose)
          cat(sprintf("  step=(%s)  |X|=%-10d time=%8.3f s  |S|=%-2d crit=%.6f conv=%s\n",
                      fmt_step(step), nrow(X), tt, nrow(res$support), res$criterion, res$converged))
      } else {
        Xpts <- .factor_refined_grid(blo, bhi, res$support, steps[[i - 1L]], step,
                                     is_factor, nlevels)
        init <- unique(.nearest_idx(Xpts, res$support))
        tt   <- system.time(
          cand <- if (use_ma)
                    ma_or(Xpts, scaled_of(Xpts), function() solve_on_grid(Xpts, init))
                  else solve_on_grid(Xpts, init))[3]
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
      cand <- if (use_ma)
                ma_or(X, scaledX, function()
                  robust_solve(X, scaledX, apply(X, 2, min), apply(X, 2, max)))
              else robust_solve(X, scaledX, apply(X, 2, min), apply(X, 2, max))
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
                grid_sizes = nrow(X), total_time = NA_real_,  # stamped below
                box_lo = apply(X, 2, min), box_hi = apply(X, 2, max), p = p,
                is_factor = is_factor, theta = theta,
                coef_names = coef_names, link = model_link,
                # certified lower bound on the efficiency vs the true optimum
                # over this candidate set (Becker & Yang, Thm 4.5 / 4.6); for
                # c-optimality Elfving's bound, which also covers singular designs
                efficiency_lower_bound =
                  if (is_c && !is.null(cand$elfving) && is.finite(cand$elfving$efficiency))
                    cand$elfving$efficiency
                  else .efficiency_bound(p, cand$max_d, cand$criterion, nrow(wb_use)),
                global_max_d = if (cand$converged) cand$max_d else NA_real_,
                global_check = cand$converged)
    if (is_c && !is.null(cand$elfving))
      out[c("elfving_bound", "elfving_gap", "elfving_h")] <-
        list(cand$elfving$bound, cand$elfving$gap, cand$elfving$h)
    if (!cand$converged)
      warning(sprintf("optimal_design(): the returned design did NOT converge (max_d = %.3e > eps0 = %g); it is not optimal. Increase max_iter, coarsen the grid, or supply a warm start.",
                      cand$max_d, eps0), call. = FALSE)
    out$total_time <- proc.time()[3] - t_start
    return(out)
  }

  # ---- design box, continuous = TRUE: grid-free search over the region -------
  # (see continuous.R).  Only the continuous covariates are searched; an
  # all-factor box has nothing to search and is enumerated as usual below.
  if (isTRUE(continuous) && any(!is_factor)) {
    box_lo <- meta$lo; box_hi <- meta$hi
    ev <- .cont_evaluator(info_mode, info_vector, info_matrix, theta_use, k, b,
                          is_factor, box_lo, box_hi, info_jacobian)
    if (!is.null(init_points)) {
      # the user's starting support
      X0 <- as.matrix(init_points); storage.mode(X0) <- "double"
      if (ncol(X0) != length(is_factor))
        stop(sprintf("'init_points' must have one column per covariate (%d).",
                     length(is_factor)), call. = FALSE)
      .validate_factor_columns(X0, is_factor, nlevels, "init_points")
      cont <- which(!is_factor)
      if (any(sweep(X0[, cont, drop = FALSE], 2, box_lo[cont], "<")) ||
          any(sweep(X0[, cont, drop = FALSE], 2, box_hi[cont], ">")))
        stop("'init_points' must lie inside the design box.", call. = FALSE)
      init_sup <- X0; init_w <- rep(1 / nrow(X0), nrow(X0))
      t_init <- 0; n_grid <- nrow(X0)
    } else {
      # the OWEA engine on a small coarse grid: init_levels levels per continuous
      # covariate (raised if that gives fewer than k + 1 points) x factor levels
      L <- max(2L, as.integer(init_levels))
      nfac <- if (any(is_factor)) prod(nlevels[is_factor]) else 1
      while (L^sum(!is_factor) * nfac < k + 1 && L < 50L) L <- L + 1L
      X0 <- .cont_initial_grid(box_lo, box_hi, is_factor, nlevels, L)
      sX <- scaled_of(X0)
      t_init <- system.time(
        g0 <- solve_prepared(X0, sX,
                             .initial_support_idx(X0, info_mode, sX, k, infor0,
                                                  init_method,
                                                  ma_max_iter = ma_max_iter)))[3]
      init_sup <- g0$support; init_w <- g0$weights; n_grid <- nrow(X0)
      if (verbose)
        cat(sprintf(paste0("  coarse start grid: %d levels per continuous covariate, ",
                           "%d points  |S|=%-2d crit=%.6f conv=%s\n"),
                    L, nrow(X0), nrow(g0$support), g0$criterion, g0$converged))
    }
    cs <- .cont_solve(ev, p, wb_use, infor0, ms, box_lo, box_hi, is_factor, nlevels,
                      init_sup, init_w, eps0, max_iter, n_starts, n_audit,
                      merge_tol, verbose)
    c_cert <- NULL
    if (is_c) {                                      # Elfving's bound over the region
      c_cert <- .c_certificate_cont(ev, cvec, infor0, box_lo, box_hi, is_factor,
                                    nlevels, cs$support, c_value(cs$criterion), eps0)
      cs$converged <- isTRUE(cs$converged) || isTRUE(c_cert$certified)
    }
    info_out <- infor0 + .opt_infor_from_support(cs$support, cs$weights, b,
                                                 info_mode, info_vector,
                                                 info_matrix, theta_use, k)
    dimnames(info_out) <- list(coef_names, coef_names)
    out <- list(support = cs$support, weights = cs$weights,
                criterion = cs$criterion, max_d = cs$max_d,
                information = info_out, converged = cs$converged,
                times = c(t_init, cs$time), grid_sizes = n_grid,
                total_time = NA_real_,   # stamped below
                box_lo = box_lo, box_hi = box_hi, p = p,
                is_factor = is_factor, theta = theta,
                coef_names = coef_names, link = model_link,
                # the bound is over the points examined by the search + audit
                # (Elfving's bound for c-optimality)
                efficiency_lower_bound =
                  if (!is.null(c_cert) && is.finite(c_cert$efficiency)) c_cert$efficiency
                  else .efficiency_bound(p, cs$max_d, cs$criterion, nrow(wb_use)),
                elfving_bound = if (is.null(c_cert)) NULL else c_cert$bound,
                elfving_gap   = if (is.null(c_cert)) NULL else c_cert$gap,
                elfving_h     = if (is.null(c_cert)) NULL else c_cert$h,
                global_max_d = cs$max_d, global_check = cs$converged,
                method = "continuous", iterations = cs$iterations,
                history = cs$history, maximiser = cs$maximiser,
                n_audit = as.integer(n_audit),
                audit_max_d = if (is.null(cs$audit)) NA_real_ else cs$audit$max_d,
                jacobian = if (!is.null(info_jacobian)) "user"
                           else if (ev$analytic) "analytic" else "finite differences")
    if (isTRUE(check_global)) {
      # an optional grid scan on top of the search, at a step the user gives
      if (is.null(global_step)) {
        warning("optimal_design(): with continuous = TRUE, check_global = TRUE needs 'global_step' (there is no step_sequence to take it from); the grid check was skipped.",
                call. = FALSE)
      } else {
        gstep <- .expand_stage_step(as.numeric(global_step), is_factor,
                                    sum(!is_factor), length(is_factor))
        npts  <- prod(ifelse(is_factor, nlevels, round((box_hi - box_lo) / gstep) + 1))
        if (npts > global_max_points) {
          warning(sprintf("optimal_design(): the whole-box grid check at step (%s) would evaluate %.0f design points (> global_max_points = %.0f); it was skipped.",
                          paste(format(gstep[!is_factor]), collapse = ","), npts,
                          global_max_points), call. = FALSE)
        } else {
          g <- global_md(cs$support, cs$weights, box_lo, box_hi, gstep)
          out$grid_check_max_d <- g$max_d
          out$global_max_d     <- max(g$max_d, cs$max_d)
          out$global_check     <- out$global_max_d <= eps0
          out$efficiency_lower_bound <- .efficiency_bound(p, out$global_max_d,
                                                          cs$criterion, nrow(wb_use))
        }
      }
    }
    if (!cs$converged)
      warning(sprintf("optimal_design(): the continuous search did NOT converge (max_d = %.3e > eps0 = %g); the design is not certified optimal. Increase max_iter or n_starts, or supply init_points.",
                      cs$max_d, eps0), call. = FALSE)
    out$total_time <- proc.time()[3] - t_start
    return(out)
  } else if (isTRUE(continuous)) {
    message("optimal_design(): continuous = TRUE, but every covariate is a factor; ",
            "the levels are enumerated as usual.")
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
              grid_sizes = res$grid_sizes, total_time = NA_real_,  # stamped below
              box_lo = box_lo, box_hi = box_hi, p = p,
              is_factor = is_factor, theta = theta,
              coef_names = coef_names, link = model_link,
              # certified lower bound on the efficiency vs the optimum over the
              # grids the design was checked on (replaced by the whole-box bound
              # below when check_global runs); Elfving's bound for c-optimality
              efficiency_lower_bound =
                if (is_c && !is.null(res$elfving) && is.finite(res$elfving$efficiency))
                  res$elfving$efficiency
                else .efficiency_bound(p, res$max_d, res$criterion, nrow(wb_use)),
              global_max_d = NA_real_, global_check = NA)
  if (is_c && !is.null(res$elfving))
    out[c("elfving_bound", "elfving_gap", "elfving_h")] <-
      list(res$elfving$bound, res$elfving$gap, res$elfving$h)

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
        out$efficiency_lower_bound <- .efficiency_bound(p, g$max_d, res$criterion,
                                                        nrow(wb_use))
        if (is_c) {                                  # Elfving's bound on the whole-box grid
          Xv <- .factor_make_grid(box_lo, box_hi, gstep, is_factor, nlevels)
          gc <- certify_c(list(support = res$support, weights = res$weights,
                               criterion = res$criterion, converged = FALSE),
                          scaled_of(Xv))
          out$global_check <- isTRUE(gc$converged)
          if (is.finite(gc$elfving$efficiency))
            out$efficiency_lower_bound <- gc$elfving$efficiency
          out[c("elfving_bound", "elfving_gap", "elfving_h")] <-
            list(gc$elfving$bound, gc$elfving$gap, gc$elfving$h)
        }
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
  # stamped last, so the check_global verification above is counted too
  out$total_time <- proc.time()[3] - t_start
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
