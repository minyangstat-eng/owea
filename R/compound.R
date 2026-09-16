# ===========================================================================
# compound.R -- compound (multi-criterion / multi-model) optimal designs.
#
# ADDITIVE FILE: nothing here modifies the single-criterion code path.  It
# reuses the existing internal helpers (.resolve_model_spec, .build_info_data,
# .wb_from, .make_infor0, the grid builders, .initial_support_idx, ...) exactly
# as optimal_design() does, and drives the separate C++ core in
# src/owea_compound.cpp.
#
# The criterion (see theory/compound-criterion.tex):
#
#     Psi_alpha(xi) = sum_j alpha_j Psi_j(xi) / Psi_j^*   (efficiency = TRUE)
#                   = sum_j alpha_j Psi_j(xi)             (efficiency = FALSE)
#
# where component j carries its OWN model -- its own per-point information
# I_j(x, theta_j), parameter dimension k_j, existing-design information
# M_{0,j}, quantity of interest g_j (dimension v_j) and type (D or A) -- and
#
#     Psi_j = det(S_j^{-1})^{1/v_j}   (p = 0, D)      S_j = G_j M_j^{-1} G_j'
#     Psi_j = v_j / tr(S_j)           (p = 1, A)      M_j = a M_{0,j} + b M_{xi,j}
#
# Both are "per parameter": positively homogeneous of degree one in M_j and
# equal to 1 at S_j = I, so components of different dimension are comparable
# and Psi_alpha is a weighted average of efficiencies lying in (0, 1].
# ===========================================================================

# ---- component specification ---------------------------------------------
# A component is a list holding the same model arguments optimal_design()
# accepts, plus its criterion p and its quantity of interest.  Everything is
# resolved here, once, into the pieces the C++ core needs.
.cmp_resolve_one <- function(cmp, jj, design_box, candidate_set, factor_levels) {
  if (!is.list(cmp))
    stop(sprintf("components[[%d]] must be a list.", jj), call. = FALSE)
  gv <- function(nm, default = NULL) if (is.null(cmp[[nm]])) default else cmp[[nm]]

  p <- .check_criterion(gv("p", 0L))

  spec <- .resolve_model_spec(gv("link"), gv("f"), gv("x"), gv("fx"), gv("xx"),
                             gv("intercept", TRUE), gv("coding", "zero-sum"),
                             design_box, candidate_set, factor_levels,
                             gv("info_vector"), gv("info_matrix"), gv("theta"),
                             ff = gv("ff"), ncat = gv("ncat"))
  info_vector <- spec$info_vector
  info_matrix <- spec$info_matrix
  theta       <- spec$theta

  if (is.null(info_matrix) && is.null(info_vector))
    stop(sprintf(paste0("components[[%d]]: supply 'info_vector', 'info_matrix', ",
                        "or a model spec ('link' + terms)."), jj), call. = FALSE)
  if (!is.null(info_matrix) && !is.null(info_vector))
    stop(sprintf("components[[%d]]: supply only one of 'info_vector' / 'info_matrix'.",
                 jj), call. = FALSE)

  info_mode <- if (is.null(info_vector)) 1L else 0L
  if (info_mode == 0L && .info_vec_needs_theta(info_vector) && is.null(theta))
    stop(sprintf(paste0("components[[%d]]: 'theta' is required when 'info_vector' ",
                        "is function(x, theta)."), jj), call. = FALSE)
  info_vector <- .normalize_info_vector(info_vector, theta)
  info_matrix <- .normalize_info_matrix(info_matrix, theta)

  list(p = p, info_mode = info_mode,
       info_vector = info_vector, info_matrix = info_matrix,
       theta = theta, coef_names = spec$coef_names, link = spec$link,
       spec_given = spec$spec_given, factor_levels = spec$factor_levels,
       wb_in = gv("wb"), subset = gv("subset"), grad_g = gv("grad_g"),
       name = gv("name", paste0("component ", jj)))
}

# k for a component, from a probe point in the design space.
.cmp_k <- function(cc, probe) {
  if (cc$info_mode == 1L) nrow(as.matrix(cc$info_matrix(probe)))
  else length(as.numeric(cc$info_vector(probe, cc$theta)))
}

# Per-component wb / infor0 / b / min_support, given the existing design.
.cmp_finalize <- function(cc, probe, xi0_points, xi0_weights, n0, n1, samp) {
  k <- .cmp_k(cc, probe)
  theta_use <- if (is.null(cc$theta)) rep(0.0, k) else as.numeric(cc$theta)
  wb <- .wb_from(cc$wb_in, cc$subset, cc$grad_g, theta_use, k)
  i0 <- .make_infor0(xi0_points, xi0_weights, n0, n1,
                     cc$info_mode, cc$info_vector, cc$info_matrix, theta_use, k)
  r_pt <- .point_info_rank(cc$info_mode, cc$info_vector, cc$info_matrix,
                           theta_use, samp)
  cc$k <- k; cc$theta_use <- theta_use; cc$wb <- wb
  cc$infor0 <- i0$infor0; cc$b <- i0$b
  cc$min_support <- .min_support_rule(i0$infor0, k, wb, r_pt)
  cc
}

# b-scaled candidate information for one component on the grid X.
.cmp_data <- function(cc, X)
  .scale_info(.build_info_data(X, cc$info_mode, cc$info_vector, cc$info_matrix,
                               cc$theta_use)$info_data,
              cc$info_mode, cc$b)

# ---- shared setup ----------------------------------------------------------
# Resolve every component's model and the shared design-space metadata, and
# finalize each component (wb, infor0, b, k, min_support) against the existing
# design.  Used by compound_design() and compound_exact_design() alike, so the
# two cannot disagree about the models or the design region.
.cmp_setup <- function(components, design_box, candidate_set, factor_levels,
                       xi0_points, xi0_weights, n0, n1) {
  J <- length(components)
  use_set <- !is.null(candidate_set)

  comps <- lapply(seq_len(J), function(j)
    .cmp_resolve_one(components[[j]], j, design_box, candidate_set, factor_levels))

  fl <- factor_levels
  if (is.null(fl)) {
    for (cc in comps) if (cc$spec_given && !is.null(cc$factor_levels)) {
      fl <- cc$factor_levels; break
    }
  }
  meta <- if (use_set) .factor_levels_to_meta(fl, ncol(as.matrix(candidate_set)))
          else .parse_design_box(design_box)
  is_factor <- meta$is_factor; nlevels <- meta$nlevels
  if (use_set) .validate_factor_columns(candidate_set, is_factor, nlevels,
                                        "candidate_set")
  if (!is.null(xi0_points))
    .validate_factor_columns(xi0_points, is_factor, nlevels, "xi0_points")

  probe <- if (use_set) as.numeric(as.matrix(candidate_set)[1, ])
           else .factor_probe_point(meta$lo, meta$hi, is_factor)
  samp <- if (use_set) {
            Xs <- as.matrix(candidate_set)
            Xs[unique(round(seq(1, nrow(Xs), length.out = min(7L, nrow(Xs))))), ,
               drop = FALSE]
          } else {
            t(vapply(c(0.5, 0.25, 0.75, 0.1, 0.9),
                     function(a) meta$lo + a * (meta$hi - meta$lo),
                     numeric(length(meta$lo))))
          }
  if (any(is_factor)) samp <- .snap_factor_levels(samp, is_factor, nlevels)

  comps <- lapply(comps, .cmp_finalize, probe = probe,
                  xi0_points = xi0_points, xi0_weights = xi0_weights,
                  n0 = n0, n1 = n1, samp = samp)

  list(comps = comps, meta = meta, is_factor = is_factor, nlevels = nlevels)
}

# ---- evaluate an ARBITRARY design ------------------------------------------
# Evaluate the compound criterion at a design (support + weights) that need not
# lie on any grid, optionally scanning a candidate set for the sensitivity.
#
# The support is appended to the scanning grid and indexed at the end, so the
# design's own information is built from the EXACT support points rather than
# from nearest grid neighbours; the sensitivity is then read off the leading
# nrow(X) entries.  (Snapping the support to the grid first would perturb both
# the criterion and the offset for an off-grid design.)
.cmp_eval <- function(comps, at, support, weights, X = NULL) {
  support <- as.matrix(support); storage.mode(support) <- "double"
  m <- nrow(support)
  if (length(weights) != m)
    stop("'weights' must have one value per support point.", call. = FALSE)

  Xall <- if (is.null(X)) support else rbind(as.matrix(X), support)
  idx  <- if (is.null(X)) seq_len(m) else nrow(as.matrix(X)) + seq_len(m)
  data <- lapply(comps, function(cc) .cmp_data(cc, Xall))

  im <- as.integer(vapply(comps, function(cc) cc$info_mode, numeric(1)))
  wb <- lapply(comps, function(cc) cc$wb)
  i0 <- lapply(comps, function(cc) cc$infor0)
  pp <- as.integer(vapply(comps, function(cc) cc$p, numeric(1)))

  psi <- compound_psi_cpp(data, im, wb, i0, pp, as.numeric(at),
                          as.integer(idx), as.numeric(weights))
  val <- compound_criterion_cpp(data, im, wb, i0, pp, as.numeric(at),
                                as.integer(idx), as.numeric(weights))
  out <- list(criterion = val, psi = psi)
  if (!is.null(X)) {
    d <- compound_dirderiv_cpp(data, im, wb, i0, pp, as.numeric(at),
                               as.integer(idx), as.numeric(weights))
    v <- compound_verify_cpp(data, im, wb, i0, pp, as.numeric(at),
                             as.integer(idx), as.numeric(weights))
    out$sensitivity <- d[seq_len(nrow(as.matrix(X)))]   # grid points only
    out$max_d       <- max(out$sensitivity)
    out$index       <- which.max(out$sensitivity)
    out$euler       <- v$euler - v$value
  }
  out
}

# ---- initial support ------------------------------------------------------
# The compound engine needs an EXPLICIT starting support (unlike the
# single-criterion engine it has no built-in IBOSS), and that support must make
# EVERY component's information matrix non-singular -- a start adequate for a
# three-parameter model can easily be singular for a four-parameter one.  So
# take the union of the per-component "minmax" starts, which is non-singular
# for each component and therefore (being a superset) for all of them.  Zero
# weights are pruned by the weight optimiser, so an over-large start is safe.
.cmp_init_idx <- function(comps, data, X) {
  idx <- integer(0)
  for (j in seq_along(comps)) {
    ij <- tryCatch(
      .initial_support_idx(X, comps[[j]]$info_mode, data[[j]], comps[[j]]$k,
                           comps[[j]]$infor0, "minmax"),
      error = function(e) integer(0))
    idx <- c(idx, ij)
  }
  idx <- unique(idx)
  ms  <- max(vapply(comps, function(cc) cc$min_support, numeric(1)))
  if (length(idx) < ms) {
    # deterministic spread over the candidate set as a last resort
    extra <- unique(round(seq(1, nrow(X), length.out = min(nrow(X), ms + 2L))))
    idx <- unique(c(idx, extra))
  }
  as.integer(idx)
}

# ---- one solve on a fixed candidate set -----------------------------------
# `comps` is the finalized component list, `at` the weights actually used
# (alpha, or alpha / psi_star).  Returns the raw C++ result plus the grid.
.cmp_solve_set <- function(comps, at, X, data = NULL, init_idx = NULL,
                           max_iter = 100L, eps0 = 1e-6, verbose = FALSE) {
  if (is.null(data)) data <- lapply(comps, function(cc) .cmp_data(cc, X))
  if (is.null(init_idx) || length(init_idx) == 0L)
    init_idx <- .cmp_init_idx(comps, data, X)
  ms <- max(vapply(comps, function(cc) cc$min_support, numeric(1)))
  res <- compound_appro_opt_cpp(
    data,
    as.integer(vapply(comps, function(cc) cc$info_mode, numeric(1))),
    lapply(comps, function(cc) cc$wb),
    lapply(comps, function(cc) cc$infor0),
    as.integer(vapply(comps, function(cc) cc$p, numeric(1))),
    as.numeric(at), as.integer(init_idx), as.integer(ms),
    as.integer(max_iter), as.numeric(eps0), isTRUE(verbose))
  res$support   <- X[res$index, , drop = FALSE]
  res$weights   <- res$weight            # name it as the rest of the package does
  res$converged <- res$sensitivity <= eps0
  res
}

# ---- multistage solve over a design box -----------------------------------
.cmp_solve_box <- function(comps, at, box_lo, box_hi, stage_by, is_factor,
                           nlevels, max_iter, eps0, accept_tol, verbose) {
  res <- NULL; times <- numeric(0); gsz <- integer(0); prev_step <- NULL
  fmt <- function(s) paste(format(s[!is_factor]), collapse = ",")
  for (i in seq_along(stage_by)) {
    step <- stage_by[[i]]
    if (i == 1L) {
      X <- .factor_make_grid(box_lo, box_hi, step, is_factor, nlevels)
      init <- NULL
    } else {
      X <- .factor_refined_grid(box_lo, box_hi, res$support, prev_step, step,
                                is_factor, nlevels)
      init <- unique(.nearest_idx(X, res$support))
    }
    tt <- system.time(
      cand <- .cmp_solve_set(comps, at, X, init_idx = init,
                             max_iter = max_iter, eps0 = eps0,
                             verbose = FALSE))[3]
    times <- c(times, tt); gsz <- c(gsz, nrow(X))
    # MAXIMISATION: keep a refinement only if it converges and does not lose.
    keep <- is.null(res) || (cand$converged && cand$value >= res$value - accept_tol)
    if (keep) { res <- cand; prev_step <- step }
    if (verbose)
      cat(sprintf("  step=(%s)  |X|=%-9d time=%7.3f s  |S|=%-2d Psi=%.8f %s\n",
                  fmt(step), nrow(X), tt, nrow(cand$support), cand$value,
                  if (keep) "" else "-> rejected, keeping previous"))
  }
  res$times <- times; res$grid_sizes <- gsz; res$final_step <- prev_step
  res
}

# ---- reference solve with the single-criterion engine ---------------------
# The reference value Psi_j^* must be the TRUE optimum of component j alone:
# an under-estimated Psi_j^* inflates every efficiency reported for that
# component above 1.  The compound engine's own single-component solve can
# stall (it starts from a minmax support and, until 0.4.0, could prune into a
# singular information matrix), so compound_design() also solves each
# component with optimal_design() -- which has its own starting rules and warm
# starts -- and keeps whichever design scores higher under the component.
# Returns NULL if optimal_design() cannot be run for this component.
.cmp_reference_single <- function(cc, design_box, step_sequence, candidate_set,
                                  factor_levels, xi0_points, xi0_weights, n0, n1,
                                  max_iter, eps0, accept_tol) {
  args <- list(p = cc$p, wb = cc$wb, theta = cc$theta_use,
               xi0_points = xi0_points, xi0_weights = xi0_weights,
               n0 = n0, n1 = n1, max_iter = max_iter, eps0 = eps0,
               accept_tol = accept_tol, verbose = FALSE)
  if (cc$info_mode == 0L) args$info_vector <- cc$info_vector
  else                    args$info_matrix <- cc$info_matrix
  if (!is.null(candidate_set)) {
    args$candidate_set <- candidate_set
    args$factor_levels <- factor_levels
  } else {
    args$design_box    <- design_box
    args$step_sequence <- step_sequence
  }
  r <- tryCatch(suppressWarnings(do.call(optimal_design, args)),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  list(support = r$support, weights = r$weights,
       converged = isTRUE(r$converged), max_d = r$max_d,
       efficiency_bound = r$efficiency_lower_bound)
}

# ONE efficiency per component: the ratio to the derived reference design times
# that reference's certified bound (a guaranteed lower bound relative to the
# TRUE component optimum).  Where no bound is available (psi_star supplied
# without reference_bound) the plain ratio is returned.
.cmp_certify <- function(eff, bound) {
  eff <- as.numeric(eff); bound <- as.numeric(bound)
  ifelse(is.finite(bound), pmin(eff, 1) * bound, eff)
}

# ---- certified bound for a reference design, on the compound Psi scale ----
# For a single component solved by the compound engine, the reported maximum
# sensitivity is  sens = (Psi/v) * E  (D)  or  Psi^2 * E_A  (A), with E / E_A the
# optimality gaps of Becker & Yang (Def. 4.3 / eq. 25) that owea reports as
# max_d.  Their Theorems 4.5 / 4.6 then give, for Psi_j = det(S^-1)^(1/v) (D)
# or v / tr(S) (A):
#     D:  Psi_j / Psi_j* >= 1 / (1 + E/v) = Psi_j / (Psi_j + sens)
#     A:  Psi_j / Psi_j* >= 1 - E_A / Phi_1 = 1 - sens / Psi_j
.cmp_ref_bound <- function(p, psi, sens) {
  if (!is.finite(psi) || !is.finite(sens) || psi <= 0) return(NA_real_)
  s <- max(as.numeric(sens), 0)
  if (as.integer(p) == 0L) psi / (psi + s) else max(0, 1 - s / psi)
}

# ---- per-component criterion on the single-criterion scale ---------------
# The compound criterion is OPTIMISED through Psi_j (see the header), but the
# per-component values are REPORTED on the scale optimal_design() uses, so the
# two entry points agree number for number:
#     D (p = 0):  log det(S_j) / v_j  = -log(Psi_j)
#     A (p = 1):  tr(S_j) / v_j       =  1 / Psi_j
# (smaller is better on this scale).  The efficiencies are unchanged:
# exp(crit* - crit) for D and crit* / crit for A both equal Psi_j / Psi_j*.
.cmp_single_scale <- function(psi, p) {
  psi <- as.numeric(psi); p <- as.numeric(p)
  ifelse(p == 0, -log(psi), 1 / psi)
}

#' Compound (multi-criterion, multi-model) optimal design.
#'
#' Finds one approximate design that is good simultaneously for several
#' criteria, each of which may belong to a \emph{different model}. The
#' components share the design region and the design; everything else --- the
#' model, the parameter dimension, the quantity of interest, the criterion
#' type --- is component-specific.
#'
#' The criterion maximised is
#' \deqn{\Psi_\alpha(\xi)=\sum_j \alpha_j \Psi_j(\xi)/\Psi_j^{*},}
#' the weighted average of the component \emph{efficiencies}, where
#' \eqn{\Psi_j^{*}} is the value attained by the design optimal for component
#' \eqn{j} alone (computed automatically over the same design space and the
#' same stage structure). Each component criterion is stated per parameter,
#' \eqn{\Psi_j=\det(S_j^{-1})^{1/v_j}} for \code{p = 0} (D) and
#' \eqn{\Psi_j=v_j/\mathrm{tr}(S_j)} for \code{p = 1} (A), with
#' \eqn{S_j=G_j M_j^{-1}G_j'}; so every component lies in \eqn{(0,1]} after
#' normalisation and \eqn{\Psi_\alpha} is directly readable as an average
#' efficiency. Set \code{efficiency = FALSE} to weight the raw \eqn{\Psi_j}
#' instead (rarely what is wanted --- see the note in \code{Details}).
#'
#' @param components a list of components, one per criterion. Each element is
#'   itself a list holding the model exactly as \code{\link{optimal_design}}
#'   takes it --- \code{info_vector} or \code{info_matrix} plus \code{theta},
#'   OR a formula-style spec (\code{link}, \code{f}, \code{x}, \code{fx},
#'   \code{xx}, \code{ff}, \code{intercept}, \code{coding}, \code{ncat}) ---
#'   plus \code{p} (0 = D, 1 = A; default 0), optionally one of \code{subset} /
#'   \code{grad_g} / \code{wb} for the quantity of interest, and optionally a
#'   \code{name} used in the printout. Different components may use entirely
#'   different models, with different numbers of parameters.
#' @param alpha nonnegative component weights, recycled to
#'   \code{length(components)} and normalised to sum to 1. Default: equal
#'   weights.
#' @param design_box,step_sequence continuous design region and grid steps, as
#'   in \code{\link{optimal_design}}. Supply these OR \code{candidate_set}.
#' @param candidate_set optional \eqn{n \times d} matrix of candidate points.
#' @param factor_levels factor specification for the \code{candidate_set} path.
#' @param efficiency if \code{TRUE} (default) each component is divided by its
#'   own optimum \eqn{\Psi_j^{*}}, so \code{alpha} weights efficiencies. If
#'   \code{FALSE} the raw \eqn{\Psi_j} are weighted.
#' @param reference_bound optional numeric vector, one value per component: the
#'   certified efficiency bound of the design each supplied \code{psi_star} came
#'   from (the \code{reference_bound} element of an earlier result). Used only
#'   together with \code{psi_star}; without it the certified efficiencies are
#'   \code{NA} when \code{psi_star} is supplied.
#' @param psi_star optional numeric vector of reference values
#'   \eqn{\Psi_j^{*}}; by default they are computed by \code{length(components)}
#'   single-component runs over the same design space (each one an ordinary
#'   single-criterion optimal design problem).
#' @param xi0_points,xi0_weights,n0,n1 an existing design to augment, as in
#'   \code{\link{optimal_design}}. The runs are shared by all components; the
#'   information each model extracts from them is not.
#' @param max_iter maximum exchange iterations.
#' @param eps0 stopping threshold on the compound sensitivity function.
#' @param accept_tol a refinement stage is kept only if it converges and does
#'   not lose more than this in \eqn{\Psi_\alpha}.
#' @param verbose print progress.
#' @details
#'   \strong{Why efficiency weighting.} The raw \eqn{\Psi_j} of different
#'   models are not comparable --- a logistic model with \eqn{k=3} and a Poisson
#'   model with \eqn{k=4} produce criterion values on unrelated scales --- so
#'   with \code{efficiency = FALSE} the \code{alpha} absorb that mismatch
#'   instead of expressing a preference. The two settings give different
#'   optima unless all \eqn{\Psi_j^{*}} coincide.
#'
#'   \strong{Reference values.} Each \eqn{\Psi_j^{*}} is itself an optimal
#'   design problem and is solved here over the same design space and with the
#'   same \code{n0}/\code{n1}/\code{xi0_points}, which is what keeps every
#'   efficiency in \eqn{(0,1]}. They are returned in \code{psi_star} and should
#'   be reported with the design: they determine the effective weights.
#'
#'   \strong{Optimality.} \eqn{\Psi_\alpha} is concave, so the returned design
#'   is a global optimum when \code{max_d} \eqn{\le} \code{eps0}, by the
#'   compound equivalence theorem: at the optimum
#'   \eqn{\sum_j [\,b\,\mathrm{tr}(K_j I_j(x)) + a\,\mathrm{tr}(K_j M_{0,j})\,]
#'   \le \Psi_\alpha} for every \eqn{x}, with equality on the support.
#'   The optimal information matrices are unique; the design attaining them
#'   need not be.
#' @return a list with \code{support}, \code{weights}, \code{criterion} (the
#'   compound value \eqn{\Psi_\alpha}), \code{max_d} (maximum compound
#'   sensitivity; \eqn{\le} \code{eps0} certifies optimality), \code{converged},
#'   \code{psi} (per-component \eqn{\Psi_j}), \code{psi_star},
#'   \code{component_criterion} and \code{component_criterion_star} (the
#'   per-component values of this design and of each component's own optimum
#'   on the scale \code{\link{optimal_design}} reports: \eqn{\log\det
#'   \Sigma_j / v_j} for D, \eqn{\mathrm{tr}\,\Sigma_j / v_j} for A, i.e.
#'   \eqn{-\log\Psi_j} and \eqn{1/\Psi_j}; smaller is better),
#'   \code{efficiency_lower_bound} (per component, a guaranteed lower bound
#'   relative to the TRUE component optimum: the ratio \eqn{\Psi_j/\Psi_j^{*}}
#'   to the derived reference design times that reference's own certified
#'   bound, Theorems 4.5 and 4.6 of Becker and Yang, \emph{Post Hoc Control
#'   Group Selection via Constrained Optimal Design}), \code{reference_bound}
#'   (those bounds), \code{criterion_bound} (the compound design's own bound
#'   \eqn{\Psi_\alpha/(\Psi_\alpha + \mathrm{max\_d})}), \code{alpha},
#'   \code{at} (the weights actually used), \code{information} (a list of the
#'   per-component information matrices \eqn{M_j}), \code{iterations},
#'   \code{times}, \code{grid_sizes} and \code{total_time}.
#' @seealso \code{\link{optimal_design}}, \code{\link{compound_sensitivity}}.
#' @examples
#' \dontrun{
#' # one design for two models at once: a logistic main-effects model
#' # (D-optimal) and a logistic model with an interaction (D-optimal)
#' f1 <- function(x, th) { q <- c(1, x[1], x[2]); e <- sum(q * th)
#'                         (exp(e / 2) / (1 + exp(e))) * q }
#' f2 <- function(x, th) { q <- c(1, x[1], x[2], x[1] * x[2]); e <- sum(q * th)
#'                         (exp(e / 2) / (1 + exp(e))) * q }
#' res <- compound_design(
#'   components = list(list(info_vector = f1, theta = c(0.5, 1, -1), p = 0),
#'                     list(info_vector = f2, theta = c(0.5, 1, -1, 0.5), p = 0)),
#'   alpha      = c(0.5, 0.5),
#'   design_box = list(c(-2, 2), c(-2, 2)),
#'   step_sequence = c(0.5, 0.1))
#' print(res)
#' }
#' @export
compound_design <- function(components, alpha = NULL,
                            design_box = NULL, step_sequence = NULL,
                            candidate_set = NULL, factor_levels = NULL,
                            efficiency = TRUE, psi_star = NULL,
                            reference_bound = NULL,
                            xi0_points = NULL, xi0_weights = numeric(0),
                            n0 = 0, n1 = 1,
                            max_iter = 100L, eps0 = 1e-6,
                            accept_tol = 1e-9, verbose = FALSE) {
  t_start <- proc.time()[3]

  if (!is.list(components) || length(components) < 1L)
    stop("'components' must be a non-empty list of component specifications.",
         call. = FALSE)
  J <- length(components)

  alpha <- if (is.null(alpha)) rep(1 / J, J) else as.numeric(alpha)
  if (length(alpha) == 1L) alpha <- rep(alpha, J)
  if (length(alpha) != J)
    stop(sprintf("'alpha' must have one weight per component (%d).", J),
         call. = FALSE)
  if (any(!is.finite(alpha)) || any(alpha < 0) || sum(alpha) <= 0)
    stop("'alpha' must be nonnegative, finite, and not all zero.", call. = FALSE)
  alpha <- alpha / sum(alpha)

  use_set <- !is.null(candidate_set)
  if (!use_set && is.null(design_box))
    stop("Supply 'candidate_set', or 'design_box' (plus 'step_sequence').",
         call. = FALSE)

  # ---- resolve the models and the shared design-space metadata ------------
  su <- .cmp_setup(components, design_box, candidate_set, factor_levels,
                   xi0_points, xi0_weights, n0, n1)
  comps <- su$comps; meta <- su$meta
  is_factor <- su$is_factor; nlevels <- su$nlevels

  # ---- the design space, as the solver will see it ------------------------
  if (use_set) {
    X <- as.matrix(candidate_set); storage.mode(X) <- "double"
    runner <- function(at, vb)
      .cmp_solve_set(comps, at, X, max_iter = max_iter, eps0 = eps0, verbose = vb)
  } else {
    nstage <- if (is.list(step_sequence)) length(step_sequence)
              else if (is.matrix(step_sequence)) nrow(step_sequence)
              else length(step_sequence)
    if (nstage == 0L) {
      if (any(!is_factor))
        stop("'step_sequence' must contain at least one grid step for ",
             "continuous covariates.", call. = FALSE)
      stage_by <- list(rep(1, length(is_factor)))
    } else {
      stage_by <- .normalize_step_sequence(step_sequence, is_factor)
      if (all(is_factor)) stage_by <- stage_by[1]
    }
    runner <- function(at, vb)
      .cmp_solve_box(comps, at, meta$lo, meta$hi, stage_by, is_factor, nlevels,
                     max_iter, eps0, accept_tol, vb)
  }

  # ---- reference values: one single-component run each --------------------
  # The reference DESIGNS are kept, not just their values: they are what makes
  # the cross-efficiency table below possible, and they cost nothing extra.
  ref_designs <- NULL
  ref_bound   <- rep(NA_real_, J)   # certified bound of each reference design
  if (isTRUE(efficiency)) {
    if (is.null(psi_star)) {
      psi_star <- numeric(J)
      ref_designs <- vector("list", J)
      for (j in seq_len(J)) {
        if (verbose) cat(sprintf("reference solve for %s ...\n", comps[[j]]$name))
        one <- comps[j]
        rj <- if (use_set)
                .cmp_solve_set(one, 1, X, max_iter = max_iter, eps0 = eps0)
              else .cmp_solve_box(one, 1, meta$lo, meta$hi, stage_by, is_factor,
                                  nlevels, max_iter, eps0, accept_tol, FALSE)
        psi_j <- rj$psi[1]; conv_j <- isTRUE(rj$converged); md_j <- rj$sensitivity
        ref_j <- list(support = rj$support, weights = rj$weights)
        bound_j <- .cmp_ref_bound(comps[[j]]$p, psi_j, rj$sensitivity)
        # cross-check with the single-criterion engine and keep the better design
        alt <- .cmp_reference_single(comps[[j]], design_box, step_sequence,
                                     candidate_set, factor_levels,
                                     xi0_points, xi0_weights, n0, n1,
                                     max_iter, eps0, accept_tol)
        if (!is.null(alt)) {
          pa <- tryCatch(.cmp_eval(one, 1, alt$support, alt$weights)$psi[1],
                         error = function(e) NA_real_)
          if (is.finite(pa) && (pa > psi_j * (1 + 1e-9) || (!conv_j && alt$converged))) {
            psi_j <- pa; ref_j <- list(support = alt$support, weights = alt$weights)
            conv_j <- alt$converged; md_j <- alt$max_d
            bound_j <- if (is.null(alt$efficiency_bound)) NA_real_ else alt$efficiency_bound
          }
        }
        psi_star[j]      <- psi_j
        ref_designs[[j]] <- ref_j
        ref_bound[j]     <- bound_j
        if (!conv_j)
          warning(sprintf(paste0("compound_design(): the reference solve for %s ",
                                 "did not converge (max_d = %.3e); Psi* and the ",
                                 "efficiencies reported for that component are ",
                                 "unreliable, and its effective weight is not ",
                                 "alpha/Psi* as intended."),
                          comps[[j]]$name, md_j), call. = FALSE)
      }
    } else {
      psi_star <- as.numeric(psi_star)
      if (length(psi_star) != J)
        stop(sprintf("'psi_star' must have one value per component (%d).", J),
             call. = FALSE)
      if (!is.null(reference_bound)) {
        reference_bound <- as.numeric(reference_bound)
        if (length(reference_bound) != J)
          stop(sprintf("'reference_bound' must have one value per component (%d).", J),
               call. = FALSE)
        ref_bound <- reference_bound
      }
    }
    if (any(!is.finite(psi_star)) || any(psi_star <= 0))
      stop("every reference value 'psi_star' must be finite and positive; a ",
           "component that cannot be estimated anywhere in the design space ",
           "must be removed.", call. = FALSE)
    at <- alpha / psi_star
  } else {
    psi_star <- rep(NA_real_, J)
    at <- alpha
  }

  # ---- the compound solve -------------------------------------------------
  if (verbose) cat("compound solve ...\n")
  res <- runner(at, verbose)

  Mlist <- lapply(seq_len(J), function(j) {
    cc <- comps[[j]]
    M <- cc$infor0 + .opt_infor_from_support(res$support, res$weights, cc$b,
                                             cc$info_mode, cc$info_vector,
                                             cc$info_matrix, cc$theta_use, cc$k)
    dimnames(M) <- list(cc$coef_names, cc$coef_names)
    M
  })
  eff <- if (isTRUE(efficiency)) res$psi / psi_star else rep(NA_real_, J)
  nms <- vapply(comps, function(cc) cc$name, character(1))
  # an efficiency above 1 is impossible against a true optimum: it means the
  # reference value for that component is too small (its solve was not optimal)
  bad <- which(is.finite(eff) & eff > 1 + 1e-6)
  if (length(bad))
    warning(sprintf(paste0("compound_design(): efficiency above 1 for %s -- the ",
                           "reference optimum Psi* of that component is not the ",
                           "true optimum, so its efficiencies are unreliable. ",
                           "Try a finer or longer step_sequence, or supply ",
                           "'psi_star' from a converged single-criterion solve."),
                    paste(nms[bad], collapse = ", ")), call. = FALSE)

  # ---- cross-efficiency: every design scored under every component --------
  # Row j is the design optimal for component j alone; the last row is this
  # compound design.  Diagonal entries are 1 by construction.  This is the
  # table that shows what a single-component design costs you elsewhere, and
  # hence what the compound design buys.
  cross <- NULL
  if (!is.null(ref_designs) && !any(vapply(ref_designs, is.null, logical(1)))) {
    rows <- lapply(ref_designs, function(d)
      .cmp_eval(comps, at, d$support, d$weights)$psi / psi_star)
    rows[[J + 1L]] <- res$psi / psi_star
    cross <- do.call(rbind, rows)
    dimnames(cross) <- list(c(paste("optimal for", nms), "THIS DESIGN"), nms)
  }

  pvec <- vapply(comps, function(cc) cc$p, numeric(1))
  out <- list(support = res$support, weights = res$weights,
              criterion = res$value, max_d = res$sensitivity,
              converged = isTRUE(res$converged),
              psi = stats::setNames(res$psi, nms),
              psi_star = stats::setNames(psi_star, nms),
              # the same values on the scale optimal_design() reports
              component_criterion      = stats::setNames(.cmp_single_scale(res$psi, pvec), nms),
              component_criterion_star = stats::setNames(.cmp_single_scale(psi_star, pvec), nms),
              # ONE efficiency per component, relative to the TRUE component
              # optimum: the ratio to the derived reference design times that
              # reference's own certified bound (Becker & Yang, Thm 4.5 / 4.6)
              efficiency_lower_bound = stats::setNames(.cmp_certify(eff, ref_bound), nms),
              reference_bound        = stats::setNames(ref_bound, nms),
              # the compound design's own bound Psi_alpha / Psi_alpha* >=
              # Psi_alpha / (Psi_alpha + max_d); used by compound_exact_design()
              criterion_bound = res$value / (res$value + max(res$sensitivity, 0)),
              alpha = stats::setNames(alpha, nms),
              at = stats::setNames(at, nms),
              cross_efficiency = if (is.null(cross)) NULL
                                 else sweep(pmin(cross, 1), 2,
                                            ifelse(is.finite(ref_bound), ref_bound, 1), `*`),
              reference_designs = if (is.null(ref_designs)) NULL
                                  else stats::setNames(ref_designs, nms),
              information = stats::setNames(Mlist, nms),
              p = vapply(comps, function(cc) cc$p, numeric(1)),
              iterations = res$iter,
              is_factor = is_factor,
              efficiency_weighted = isTRUE(efficiency),
              times = if (is.null(res$times)) NA_real_ else res$times,
              grid_sizes = if (is.null(res$grid_sizes)) nrow(res$support)
                           else res$grid_sizes,
              total_time = proc.time()[3] - t_start)
  class(out) <- "compound_design"

  if (!out$converged)
    warning(sprintf(paste0("compound_design(): the returned design did NOT ",
                           "converge (max_d = %.3e > eps0 = %g); it is not ",
                           "optimal. Increase max_iter, or coarsen the grid."),
                    out$max_d, eps0), call. = FALSE)
  out
}

#' Compound criterion value of a given design.
#'
#' Evaluates the compound criterion \eqn{\Psi_\alpha} at \emph{any} design ---
#' one you wrote down yourself, one taken from the literature, a standard
#' factorial, or a \code{\link{compound_design}} result --- together with the
#' per-component values \eqn{\Psi_j}, the component efficiencies and the
#' per-component information matrices. Optionally it also scans a candidate set
#' and reports the maximum compound sensitivity, which certifies optimality.
#'
#' This is the compound counterpart of
#' \code{\link{verify_optimality}(\dots, criterion_only = TRUE)}.
#'
#' @param support \eqn{m \times d} matrix of support points (one per row), OR a
#'   \code{"compound_design"} result (its support and weights are used).
#' @param weights length-\eqn{m} vector of weights (nonnegative, summing to 1);
#'   ignored when \code{support} is a design result.
#' @param components the component list, exactly as in
#'   \code{\link{compound_design}}.
#' @param alpha component weights; default equal. Ignored when \code{support} is
#'   a \code{"compound_design"} result and \code{alpha} is not given, in which
#'   case the result's own weights are reused.
#' @param reference_bound optional certified efficiency bounds of the designs
#'   the \code{psi_star} came from (one per component); taken from a
#'   \code{compound_design} object when one is passed as \code{support}, and
#'   computed alongside \code{psi_star} when the reference solves are run here.
#'   The reported \code{efficiency_lower_bound} is the ratio to the reference
#'   times this bound; without it the plain ratio is reported.
#' @param psi_star reference values \eqn{\Psi_j^{*}}. Required for
#'   \code{efficiency = TRUE} unless a design space (\code{candidate_set}, or
#'   \code{design_box} + \code{step}) is supplied, in which case they are
#'   computed there; reused automatically from a \code{"compound_design"}
#'   result.
#' @param efficiency weight efficiencies (\code{TRUE}, default) or the raw
#'   \eqn{\Psi_j} (\code{FALSE}).
#' @param candidate_set,design_box,step a design space to scan for the maximum
#'   sensitivity, and over which any missing \code{psi_star} are computed.
#'   Omit both to get the criterion value only.
#' @param factor_levels factor specification for \code{candidate_set}.
#' @param xi0_points,xi0_weights,n0,n1 an existing design, as in
#'   \code{\link{compound_design}}.
#' @param tol threshold on the maximum sensitivity below which the design is
#'   declared optimal.
#' @param max_points safety cap (default \code{1e6}) on the number of design
#'   points the \code{design_box} + \code{step} grid may generate, exactly as in
#'   \code{\link{verify_optimality}}. The sensitivity scan is what makes a fine
#'   \code{step} expensive here, and the grid grows as the reciprocal of the
#'   step raised to the number of continuous covariates. If the grid would
#'   exceed the cap you are asked (interactive session) whether to abort, build
#'   it anyway, or compute the criterion only; non-interactively you are stopped
#'   with a message to coarsen \code{step}, raise \code{max_points}, or set
#'   \code{criterion_only = TRUE}.
#' @param criterion_only if \code{TRUE}, skip the sensitivity scan entirely: no
#'   grid is built and the design is scored on its own support, so the result
#'   carries \code{criterion}, \code{psi} and \code{efficiency_lower_bound} but none of
#'   \code{sensitivity}, \code{max_d}, \code{maximiser} or \code{is_optimal}.
#'   Needs \code{psi_star} (or \code{efficiency = FALSE}), since without a
#'   design space the reference values cannot be computed.
#' @param ... ignored.
#' @return a list with \code{criterion} (\eqn{\Psi_\alpha}), \code{psi}
#'   (per-component \eqn{\Psi_j}), \code{component_criterion} and
#'   \code{component_criterion_star} (the same values on the scale
#'   \code{\link{optimal_design}} reports, see \code{\link{compound_design}}),
#'   \code{efficiency_lower_bound} (per component: the ratio to the reference
#'   times the reference design's certified bound, a guaranteed lower bound
#'   relative to the TRUE component optimum; the plain ratio when no bound is
#'   available), \code{reference_bound}, \code{psi_star},
#'   \code{alpha}, \code{at}, \code{information} (per-component \eqn{M_j}) and,
#'   when a design space was supplied, \code{sensitivity} (one value per
#'   candidate), \code{max_d}, \code{maximiser}, \code{is_optimal} and
#'   \code{euler} (the Euler-identity residual, numerically zero).
#' @details The support need not lie on the scanning grid: it is evaluated at
#'   the exact points supplied, so an off-grid design is scored correctly.
#' @seealso \code{\link{compound_design}}, \code{\link{verify_optimality}}.
#' @examples
#' \dontrun{
#' f <- function(x, th) { q <- c(1, x[1], x[2]); e <- sum(q * th)
#'                        (exp(e / 2) / (1 + exp(e))) * q }
#' cmp <- list(list(info_vector = f, theta = c(0.5, 1, -1), p = 0),
#'             list(info_vector = f, theta = c(0.5, 1, -1), p = 1,
#'                  subset = c(2, 3)))
#' X   <- as.matrix(expand.grid(seq(-2, 2, 0.1), seq(-2, 2, 0.1)))
#' # score a 2^2 factorial under the compound criterion
#' fac <- rbind(c(-2, -2), c(-2, 2), c(2, -2), c(2, 2))
#' compound_criterion(fac, rep(0.25, 4), cmp, candidate_set = X)
#' }
#' @export
compound_criterion <- function(support, weights = NULL, components,
                               alpha = NULL, psi_star = NULL,
                               reference_bound = NULL,
                               efficiency = TRUE,
                               candidate_set = NULL, design_box = NULL,
                               step = NULL, factor_levels = NULL,
                               xi0_points = NULL, xi0_weights = numeric(0),
                               n0 = 0, n1 = 1, tol = 1e-6,
                               max_points = 1e6, criterion_only = FALSE, ...) {
  # accept a compound_design() result directly
  if (inherits(support, "compound_design")) {
    if (is.null(weights))  weights  <- support$weights
    if (is.null(alpha))    alpha    <- support$alpha
    if (is.null(psi_star) && isTRUE(support$efficiency_weighted))
      psi_star <- support$psi_star
    if (is.null(reference_bound) && !is.null(support$reference_bound))
      reference_bound <- support$reference_bound
    support <- support$support
  }
  if (is.null(weights)) stop("'weights' is required.", call. = FALSE)
  support <- as.matrix(support); storage.mode(support) <- "double"
  weights <- as.numeric(weights)
  if (any(weights < 0)) stop("'weights' must be nonnegative.", call. = FALSE)
  if (abs(sum(weights) - 1) > 1e-8) {
    warning("weights do not sum to 1; normalising.", call. = FALSE)
    weights <- weights / sum(weights)
  }

  if (!is.list(components) || length(components) < 1L)
    stop("'components' must be a non-empty list.", call. = FALSE)
  J <- length(components)
  alpha <- if (is.null(alpha)) rep(1 / J, J) else as.numeric(alpha)
  if (length(alpha) == 1L) alpha <- rep(alpha, J)
  if (length(alpha) != J)
    stop(sprintf("'alpha' must have one weight per component (%d).", J),
         call. = FALSE)
  if (any(!is.finite(alpha)) || any(alpha < 0) || sum(alpha) <= 0)
    stop("'alpha' must be nonnegative, finite, and not all zero.", call. = FALSE)
  alpha <- alpha / sum(alpha)

  # the design space to scan (optional).  Guard an enormous grid the way
  # verify_optimality() does -- count the points WITHOUT building them, since
  # it is the scan, not the criterion, that a fine step makes expensive.
  box_path <- is.null(candidate_set) && !is.null(design_box) && !is.null(step)
  if (box_path && !isTRUE(criterion_only)) {
    m0   <- .parse_design_box(design_box)
    by   <- .expand_stage_step(as.numeric(step), m0$is_factor,
                               sum(!m0$is_factor), length(design_box))
    npts <- prod(ifelse(m0$is_factor, m0$nlevels,
                        round((m0$hi - m0$lo) / by) + 1))
    if (npts > max_points) {
      msg <- sprintf(paste0("compound_criterion(): the design_box + step grid ",
                            "has %.0f design points (> max_points = %.0f)."),
                     npts, max_points)
      if (interactive()) {
        choice <- utils::menu(
          c("Abort",
            "Proceed anyway (build the full grid and check optimality)",
            "Criterion only (skip the sensitivity scan; no grid is built)"),
          title = paste(msg, "What would you like to do?"))
        if (choice == 2L) warning(paste(msg, "Proceeding anyway."), call. = FALSE)
        else if (choice == 3L) criterion_only <- TRUE
        else stop(msg, call. = FALSE)
      } else
        stop(paste(msg, "Increase 'step' (a coarser grid) to reduce the number",
                   "of design points, raise 'max_points', or set",
                   "criterion_only = TRUE to skip the sensitivity scan."),
             call. = FALSE)
    }
  }
  X <- if (!is.null(candidate_set)) {
         as.matrix(candidate_set)
       } else if (box_path && !isTRUE(criterion_only)) {
         candidate_grid(design_box, step)
       } else NULL
  if (!is.null(X)) storage.mode(X) <- "double"

  # resolve the components against whatever describes the covariates
  ref_box <- if (is.null(candidate_set)) design_box else NULL
  ref_set <- if (is.null(candidate_set)) NULL else as.matrix(candidate_set)
  if (is.null(ref_box) && is.null(ref_set)) ref_set <- support
  comps <- lapply(seq_len(J), function(j)
    .cmp_resolve_one(components[[j]], j, ref_box, ref_set, factor_levels))
  probe <- as.numeric(support[1, ])
  samp  <- support[unique(round(seq(1, nrow(support),
                                    length.out = min(7L, nrow(support))))), ,
                   drop = FALSE]
  comps <- lapply(comps, .cmp_finalize, probe = probe,
                  xi0_points = xi0_points, xi0_weights = xi0_weights,
                  n0 = n0, n1 = n1, samp = samp)

  # reference values
  if (isTRUE(efficiency)) {
    if (is.null(psi_star)) {
      if (is.null(X))
        stop("efficiency = TRUE needs the reference values: supply 'psi_star', ",
             "or a design space ('candidate_set', or 'design_box' + 'step') ",
             "over which they can be computed",
             if (isTRUE(criterion_only))
               " -- a criterion_only run builds no design space, so it needs ",
             if (isTRUE(criterion_only)) "'psi_star'",
             ".", call. = FALSE)
      refs <- lapply(seq_len(J), function(j) .cmp_solve_set(comps[j], 1, X))
      psi_star <- vapply(refs, function(r) r$psi[1], numeric(1))
      if (is.null(reference_bound))
        reference_bound <- vapply(seq_len(J), function(j)
          .cmp_ref_bound(comps[[j]]$p, refs[[j]]$psi[1], refs[[j]]$sensitivity),
          numeric(1))
    }
    psi_star <- as.numeric(psi_star)
    if (length(psi_star) != J)
      stop(sprintf("'psi_star' must have one value per component (%d).", J),
           call. = FALSE)
    if (any(!is.finite(psi_star)) || any(psi_star <= 0))
      stop("every 'psi_star' must be finite and positive.", call. = FALSE)
    at <- alpha / psi_star
  } else {
    psi_star <- rep(NA_real_, J)
    at <- alpha
  }

  ev  <- .cmp_eval(comps, at, support, weights, X)
  nms <- vapply(comps, function(cc) cc$name, character(1))

  Mlist <- lapply(seq_len(J), function(j) {
    cc <- comps[[j]]
    M <- cc$infor0 + .opt_infor_from_support(support, weights, cc$b,
                                             cc$info_mode, cc$info_vector,
                                             cc$info_matrix, cc$theta_use, cc$k)
    dimnames(M) <- list(cc$coef_names, cc$coef_names)
    M
  })

  pvec <- vapply(comps, function(cc) cc$p, numeric(1))
  out <- list(criterion = ev$criterion,
              psi = stats::setNames(ev$psi, nms),
              # the same values on the scale optimal_design() reports
              component_criterion      = stats::setNames(.cmp_single_scale(ev$psi, pvec), nms),
              component_criterion_star = stats::setNames(.cmp_single_scale(psi_star, pvec), nms),
              # ONE efficiency per component: the ratio to the reference times the
              # reference's certified bound when available (Becker & Yang, Thm
              # 4.5 / 4.6), else the plain ratio
              efficiency_lower_bound = stats::setNames(
                if (isTRUE(efficiency))
                  .cmp_certify(ev$psi / psi_star,
                               if (is.null(reference_bound)) rep(NA_real_, J)
                               else as.numeric(reference_bound))
                else rep(NA_real_, J), nms),
              reference_bound = stats::setNames(
                if (is.null(reference_bound)) rep(NA_real_, J) else as.numeric(reference_bound), nms),
              psi_star = stats::setNames(psi_star, nms),
              alpha = stats::setNames(alpha, nms),
              at = stats::setNames(at, nms),
              information = stats::setNames(Mlist, nms),
              support = support, weights = weights)
  if (!is.null(X)) {
    out$sensitivity <- ev$sensitivity
    out$max_d       <- ev$max_d
    out$maximiser   <- X[ev$index, ]
    out$is_optimal  <- ev$max_d <= tol
    out$euler       <- ev$euler
  }
  out
}

#' Compound sensitivity function of a given design.
#'
#' Evaluates the compound sensitivity function over a candidate set --- the
#' quantity the compound equivalence theorem bounds by zero. Useful for checking
#' a design the package did not compute, or for verifying a result on a finer
#' grid than the one it was found on.
#'
#' A thin wrapper on \code{\link{compound_criterion}}, kept for readability when
#' the sensitivity is what you are after.
#'
#' @param object a \code{"compound_design"} result, or an \eqn{m \times d}
#'   support matrix (then supply \code{weights}).
#' @param components the component list, as in \code{\link{compound_design}}.
#' @param candidate_set the candidate set to scan.
#' @param weights weights, when \code{object} is a bare support matrix.
#' @param ... further arguments passed to \code{\link{compound_criterion}}
#'   (\code{alpha}, \code{psi_star}, \code{efficiency}, \code{xi0_points}, ...).
#' @return a list with \code{sensitivity} (one value per candidate),
#'   \code{max_d}, \code{index} (the maximising candidate), \code{criterion} and
#'   \code{euler} (the Euler identity residual, which should be ~0).
#' @seealso \code{\link{compound_criterion}}, \code{\link{compound_design}}.
#' @export
compound_sensitivity <- function(object, components, candidate_set,
                                 weights = NULL, ...) {
  r <- compound_criterion(object, weights = weights, components = components,
                          candidate_set = candidate_set, ...)
  if (is.null(r$sensitivity))
    stop("a 'candidate_set' is required to evaluate the sensitivity.",
         call. = FALSE)
  list(sensitivity = r$sensitivity, max_d = r$max_d,
       index = which.max(r$sensitivity), criterion = r$criterion,
       euler = r$euler)
}

#' @param x a \code{"compound_design"} result.
#' @param ... ignored.
#' @rdname compound_design
#' @export
print.compound_design <- function(x, ...) {
  cat("Compound optimal design\n")
  cat(sprintf("  components : %d   (%s-weighted)\n", length(x$alpha),
              if (isTRUE(x$efficiency_weighted)) "efficiency" else "raw"))
  cat(sprintf("  criterion  : Psi_alpha = %.8f%s\n", x$criterion,
              if (isTRUE(x$efficiency_weighted))
                paste0("   (weighted average of the component efficiencies: the best any",
                       "\n                single design can reach here; it would equal 1 only if",
                       "\n                one design were optimal for every component at once)")
              else ""))
  cat(sprintf("  max_d      : %.3e   %s\n", x$max_d,
              if (isTRUE(x$converged)) "(optimal)" else "(NOT converged)"))
  cat("\n  per-component criterion values, on the scale optimal_design() reports",
      "\n  (D: log det Sigma / v;  A: tr(Sigma) / v;  smaller is better);",
      "\n  efficiency_lower_bound = guaranteed lower bound relative to that component's",
      "\n  TRUE optimum\n")
  tab <- data.frame(alpha = round(x$alpha, 4),
                    p = ifelse(x$p == 0, "D", "A"),
                    criterion = signif(x$component_criterion, 7))
  if (isTRUE(x$efficiency_weighted)) {
    tab$optimal                <- signif(x$component_criterion_star, 7)
    tab$efficiency_lower_bound <- round(x$efficiency_lower_bound, 6)
  }
  print(tab)

  if (!is.null(x$cross_efficiency)) {
    cat("\n  efficiency of each design under every component",
        "\n  (rows: design;  columns: component;  lower bounds relative to each",
        "\n   component's true optimum)\n")
    print(round(x$cross_efficiency, 3))
  }

  cat("\n  design (support points and weights):\n")
  sup <- x$support
  d <- data.frame(sup, weight = round(x$weights, 6))
  names(d) <- c(paste0("x", seq_len(ncol(sup))), "weight")
  print(d, row.names = FALSE)
  invisible(x)
}
