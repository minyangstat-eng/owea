# ===========================================================================
# verify.R -- check whether a GIVEN approximate design is optimal.
#
# Given a design (support + weights), a criterion (p) and a design space
# (design_box + step, or a candidate_set), verify_optimality() reports the
# maximum of the sensitivity (directional-derivative) function over the space
# (0 at the optimum, by the equivalence theorem), the criterion value, and the
# resulting information matrix -- after validating the support lies in the space
# and the weights are nonnegative and sum to 1.  It reuses the same internals as
# optimal_design(); the max-sensitivity computation mirrors solve.R's global_md.
# ===========================================================================

#' Verify the optimality of a given approximate design.
#'
#' Checks whether a supplied approximate design (support points and weights) is
#' optimal for a given criterion over a given design space, using the general
#' equivalence theorem. It accepts the same model and quantity-of-interest
#' arguments as \code{\link{optimal_design}}, plus the design to check and the
#' design space (a \code{design_box} + \code{step}, as in
#' \code{\link{candidate_grid}}, or an explicit \code{candidate_set}).
#'
#' The support/weights are typically the \code{support} and \code{weights}
#' returned by \code{\link{optimal_design}}. The design is validated first: the
#' support points must lie in the design space (continuous coordinates within the
#' box, factor coordinates valid levels; or, for a \code{candidate_set}, each
#' support point must be one of the candidate points) and the weights must be
#' nonnegative and sum to 1.
#'
#' @param support \eqn{m \times d} matrix of support points (one per row), OR an
#'   \code{\link{optimal_design}}/\code{\link{exact_design}} result (its
#'   \code{support}/\code{weights} are used).
#' @param weights length-\eqn{m} vector of weights (nonnegative, summing to 1);
#'   ignored when \code{support} is a design result.
#' @param design_box,step design space as a box plus a grid step, exactly as in
#'   \code{\link{candidate_grid}} (a \code{c(lo, hi)} pair is continuous; a
#'   single integer \code{L >= 2} is a factor). Supply this OR \code{candidate_set}.
#' @param candidate_set \eqn{n \times d} matrix of candidate points defining a
#'   discrete design space. Supply this OR \code{design_box} + \code{step}.
#' @param info_vector,info_matrix,theta,link,f,x,fx,xx,ff,intercept,coding,ncat
#'   the model, exactly as in \code{\link{optimal_design}}.
#' @param p criterion: \code{0} = D-optimality or \code{1} = A-optimality (only
#'   these two are supported).
#' @param wb,subset,grad_g the quantity of interest, exactly as in
#'   \code{\link{optimal_design}}.
#' @param xi0_points,xi0_weights,n0,n1 an existing design to combine with (as in
#'   \code{\link{optimal_design}}); \code{n0 = 0} (default) is single-stage.
#' @param factor_levels factor specification for the \code{candidate_set} path.
#' @param tol threshold on the maximum sensitivity below which the design is
#'   declared optimal (default \code{1e-6}).
#' @param max_points safety cap (default \code{1e6}) on the number of design
#'   points the \code{design_box} + \code{step} grid may generate. If the grid
#'   would exceed it, you are asked (interactive session) whether to abort,
#'   build the grid anyway, or compute the criterion only (see
#'   \code{criterion_only}); non-interactively you are stopped with a message to
#'   coarsen \code{step}, raise \code{max_points}, or set
#'   \code{criterion_only = TRUE}.
#' @param criterion_only if \code{TRUE}, skip the max-sensitivity
#'   (equivalence-theorem) scan entirely: the design is still validated exactly
#'   as usual (weights nonnegative and summing to 1; support points checked to
#'   lie among the design points, with a warning if not) and the criterion value
#'   and information matrix are returned, but the \code{design_box} +
#'   \code{step} grid is never built, so this works for any grid size.
#'   \code{max_sensitivity} is \code{NA}, \code{is_optimal$value} is \code{NA}
#'   (optimality is NOT assessed) and \code{maximiser} is \code{NULL}.
#' @details The design space is the finite set \code{candidate_grid(design_box,
#'   step)} (or the \code{candidate_set}). Support points that are not among those
#'   design points are allowed but flagged with a warning; choose \code{step} so
#'   the design's support fall on the grid, or pass a \code{candidate_set}.
#' @return a list with \code{is_optimal} -- itself a list with \code{value}
#'   (logical, \code{max_sensitivity <= tol}) and \code{note} (a message stating
#'   the \code{tol} used) -- plus \code{max_sensitivity} (the maximum directional
#'   derivative over the design space; \eqn{\le 0}, i.e. \eqn{\le}\code{tol}, at
#'   the optimum), \code{criterion} (the normalised Phi_p value),
#'   \code{information} (the resulting \eqn{k \times k} information matrix;
#'   combined with any existing design), \code{maximiser} (the design-space point
#'   attaining the maximum sensitivity), and echoes of \code{support},
#'   \code{weights}, \code{p}, \code{tol}. Under \code{criterion_only} the
#'   sensitivity fields are placeholders: \code{is_optimal$value} and
#'   \code{max_sensitivity} are \code{NA} and \code{maximiser} is \code{NULL}.
#' @seealso \code{\link{optimal_design}}, \code{\link{candidate_grid}},
#'   \code{\link{design_information}}.
#' @export
verify_optimality <- function(support, weights = NULL,
                              design_box = NULL, step = NULL, candidate_set = NULL,
                              info_vector = NULL, info_matrix = NULL, theta = NULL,
                              p = 0L,
                              link = NULL, f = NULL, x = NULL, fx = NULL, xx = NULL,
                              ff = NULL, intercept = TRUE, coding = "zero-sum",
                              ncat = NULL,
                              wb = NULL, subset = NULL, grad_g = NULL,
                              xi0_points = NULL, xi0_weights = numeric(0),
                              n0 = 0, n1 = 1, factor_levels = NULL, tol = 1e-6,
                              max_points = 1e6, criterion_only = FALSE) {
  p <- .check_criterion(p)

  # accept an optimal_design()/exact_design() result directly
  if (is.list(support) && !is.null(support$support) && !is.null(support$weights)) {
    if (is.null(weights)) weights <- support$weights
    support <- support$support
  }
  if (is.null(weights)) stop("'weights' is required.", call. = FALSE)
  support <- as.matrix(support); storage.mode(support) <- "double"
  weights <- as.numeric(weights)

  # ---- resolve the model (same three ways as optimal_design) -------------
  spec <- .resolve_model_spec(link, f, x, fx, xx, intercept, coding,
                              design_box, candidate_set, factor_levels,
                              info_vector, info_matrix, theta, ff = ff, ncat = ncat)
  info_vector <- spec$info_vector; info_matrix <- spec$info_matrix
  theta <- spec$theta
  if (spec$spec_given && is.null(factor_levels)) factor_levels <- spec$factor_levels
  coef_names <- spec$coef_names

  if (is.null(info_matrix) && is.null(info_vector))
    stop("Supply 'info_matrix', 'info_vector', or a model spec ('link' + terms).",
         call. = FALSE)
  if (!is.null(info_matrix) && !is.null(info_vector))
    stop("Supply only one of 'info_matrix' / 'info_vector'.", call. = FALSE)
  info_mode <- if (is.null(info_vector)) 1L else 0L
  if (info_mode == 0L && .info_vec_needs_theta(info_vector) && is.null(theta))
    stop("'theta' is required when 'info_vector' is function(x, theta).",
         call. = FALSE)
  info_vector <- .normalize_info_vector(info_vector, theta)
  info_matrix <- .normalize_info_matrix(info_matrix, theta)

  # ---- design space X + factor metadata ----------------------------------
  use_set <- !is.null(candidate_set)
  if (use_set) {
    X <- as.matrix(candidate_set); storage.mode(X) <- "double"
    meta <- .factor_levels_to_meta(factor_levels, ncol(X))
    is_factor <- meta$is_factor; nlevels <- meta$nlevels
    lo <- NULL; hi <- NULL
  } else {
    if (is.null(design_box) || is.null(step))
      stop("Supply a design space: 'candidate_set', or 'design_box' + 'step'.",
           call. = FALSE)
    m <- .parse_design_box(design_box)
    is_factor <- m$is_factor; nlevels <- m$nlevels; lo <- m$lo; hi <- m$hi
    # guard against an enormous grid (count the points without building them,
    # as optimal_design's check_global does).  criterion_only never builds it.
    by   <- .expand_stage_step(as.numeric(step), is_factor, sum(!is_factor),
                               length(design_box))
    npts <- prod(ifelse(is_factor, nlevels, round((hi - lo) / by) + 1))
    if (npts > max_points && !criterion_only) {
      msg <- sprintf(paste0("verify_optimality(): the design_box + step grid has ",
                            "%.0f design points (> max_points = %.0f)."),
                     npts, max_points)
      if (interactive()) {
        choice <- utils::menu(
          c("Abort",
            "Proceed anyway (build the full grid and check optimality)",
            "Criterion only (skip the max-sensitivity check; no grid is built)"),
          title = paste(msg, "What would you like to do?"))
        if (choice == 2L) warning(paste(msg, "Proceeding anyway."), call. = FALSE)
        else if (choice == 3L) criterion_only <- TRUE
        else stop(msg, call. = FALSE)
      } else
        stop(paste(msg, "Increase 'step' (a coarser grid) to reduce the number ",
                   "of design points, raise 'max_points', or set ",
                   "criterion_only = TRUE to skip the max-sensitivity check."),
             call. = FALSE)
    }
    if (!criterion_only) X <- candidate_grid(design_box, step)
  }
  have_X <- use_set || !criterion_only        # criterion_only box path: no grid
  if (have_X) storage.mode(X) <- "double"

  # ---- infer k and theta_use (mirrors optimal_design) --------------------
  # without a grid, probe at the first grid point it WOULD have (factor level 1,
  # continuous lo), which is exactly X[1, ] when the grid is built
  probe <- if (have_X) as.numeric(X[1, ]) else as.numeric(ifelse(is_factor, 1, lo))
  k <- if (info_mode == 1L) nrow(as.matrix(info_matrix(probe)))
       else length(as.numeric(info_vector(probe, theta)))
  theta_use <- if (is.null(theta)) rep(0.0, k) else as.numeric(theta)

  # ---- validate the given design -----------------------------------------
  d_space <- if (have_X) ncol(X) else length(design_box)
  if (ncol(support) != d_space)
    stop(sprintf("'support' has %d column(s) but the design space has %d.",
                 ncol(support), d_space), call. = FALSE)
  if (length(weights) != nrow(support))
    stop(sprintf("length(weights) (%d) must equal the number of support points (%d).",
                 length(weights), nrow(support)), call. = FALSE)
  if (any(!is.finite(weights)) || any(weights < -1e-10))
    stop("weights must be finite and nonnegative.", call. = FALSE)
  if (abs(sum(weights) - 1) > 1e-8)
    stop(sprintf("weights must sum to 1 (they sum to %.10g).", sum(weights)),
         call. = FALSE)
  .validate_factor_columns(support, is_factor, nlevels, "support")
  # Are the support points actual design points (grid points / candidate_set
  # members)?  If not, warn but proceed -- the design is still verified as given.
  # A box grid is a Cartesian product, so its nearest point is found coordinate
  # by coordinate (.offgrid_dmax) -- no scan over the N grid rows; an arbitrary
  # candidate_set has no such structure, so there the rows are searched.
  dmax <- if (use_set) {
    idx <- .nearest_idx(X, support)
    vapply(seq_len(nrow(support)),
           function(i) max(abs(support[i, ] - X[idx[i], ])), numeric(1))
  } else .offgrid_dmax(support, lo, hi, by, is_factor, nlevels)
  bad  <- which(dmax > 1e-6)
  if (length(bad)) {
    where <- if (use_set) "the candidate_set"
             else "the design points generated by design_box + step"
    warning(sprintf(paste0("%d support point(s) (e.g. #%d) are not among %s; verifying ",
                           "the design as given anyway."), length(bad), bad[1], where),
            call. = FALSE)
  }

  # ---- assemble the C++-engine pieces (as in solve.R) --------------------
  wb_use <- .wb_from(wb, subset, grad_g, theta_use, k)
  i0 <- .make_infor0(xi0_points, xi0_weights, n0, n1,
                     info_mode, info_vector, info_matrix, theta_use, k)
  infor0 <- i0$infor0; b <- i0$b

  # ---- max sensitivity over the design space (equivalence theorem) -------
  # (skipped under criterion_only: it is the only part that needs the grid)
  oi <- .opt_infor_from_support(support, weights, b, info_mode, info_vector,
                                info_matrix, theta_use, k)
  ve <- if (criterion_only) NULL else {
    scaled <- .scale_info(.build_info_data(X, info_mode, info_vector,
                                           info_matrix, theta_use)$info_data,
                          info_mode, b)
    verify_equiv_cpp(as.integer(p), as.matrix(wb_use), info_mode, scaled,
                     oi, as.matrix(infor0))
  }

  # ---- criterion value + combined information matrix ---------------------
  M <- infor0 + oi
  if (!is.null(coef_names) && length(coef_names) == k)
    dimnames(M) <- list(coef_names, coef_names)
  crit <- criterion_cpp(as.integer(p), 1L, 1.0, 1L,
                        matrix(as.numeric(M), ncol = 1L),
                        as.matrix(wb_use), matrix(0.0, k, k))

  is_opt <- if (criterion_only)
    list(value = NA,
         note  = paste0("criterion_only: the max-sensitivity check was ",
                        "skipped, so optimality was NOT assessed."))
  else
    list(value = isTRUE(ve$max_d <= tol),
         note  = sprintf(paste0("'value' is TRUE when max_sensitivity <= tol; ",
                                "tol = %g was used."), tol))
  list(is_optimal      = is_opt,
       max_sensitivity = if (criterion_only) NA_real_ else ve$max_d,
       criterion       = crit,
       information     = M,
       maximiser       = if (criterion_only) NULL else X[ve$index, ],
       p = p, support = support, weights = weights, tol = tol)
}
