# ===========================================================================
# problem.R -- the DesignProblem container for a fixed candidate set.
# ===========================================================================

#' Construct a design problem on a fixed candidate set.
#'
#' The model enters through ONE per-point function: supply either
#' \code{info_matrix(x)} returning the \eqn{k \times k} Fisher information
#' matrix at \eqn{x} (the general path), or \code{info_vector(x, theta)}
#' returning the length-k information vector \eqn{f} (per-point information
#' \eqn{f f'}, the fast path).
#'
#' @param X \eqn{n \times N} numeric matrix of candidate design points (one per
#'   row).
#' @param theta numeric parameter vector. Required only when the model function
#'   is written as \code{function(x, theta)}, or when using \code{grad_g} or an
#'   existing design; otherwise optional.
#' @param info_matrix information-matrix model, either \code{function(x)} or
#'   \code{function(x, theta)}, returning the \eqn{k \times k} matrix (general
#'   path). Supply this OR \code{info_vector}.
#' @param p criterion: 0 = D-optimality or 1 = A-optimality (only these two are
#'   supported).
#' @param subset integer vector of parameter indices of interest (partial
#'   parameters); a convenience alternative to \code{grad_g}.
#' @param grad_g function \code{theta -> v x k} matrix \eqn{dg/d\theta^T} for a
#'   general differentiable function \eqn{g(\theta)}; default identity.
#' @param xi0_points \eqn{n_0 \times N} support of an already-run design
#'   (multistage); \code{NULL} for a single-stage design.
#' @param xi0_weights weights of that existing design.
#' @param n0,n1 sample sizes of the existing and the new stage.
#' @param info_vector information-vector model, either \code{function(x)} or
#'   \code{function(x, theta)}, returning the length-k vector \eqn{f} (fast
#'   path). Supply this OR \code{info_matrix}.
#' @return an object of class \code{"DesignProblem"}.
#' @export
DesignProblem <- function(X, theta = NULL, info_matrix = NULL, p = 0L,
                          subset = NULL, grad_g = NULL,
                          xi0_points = NULL, xi0_weights = numeric(0),
                          n0 = 0, n1 = 1, info_vector = NULL) {
  if (is.null(info_matrix) && is.null(info_vector))
    stop("Supply 'info_matrix' or 'info_vector'.", call. = FALSE)
  if (!is.null(info_matrix) && !is.null(info_vector))
    stop("Supply only one of 'info_matrix' / 'info_vector'.", call. = FALSE)
  p <- .check_criterion(p)
  info_mode <- if (is.null(info_vector)) 1L else 0L

  X <- as.matrix(X); storage.mode(X) <- "double"
  if (info_mode == 0L && .info_vec_needs_theta(info_vector) && is.null(theta))
    stop("'theta' is required when 'info_vector' is function(x, theta).",
         call. = FALSE)
  theta <- if (is.null(theta)) NULL else as.numeric(theta)

  # accept info_vector(x) / (x, theta) and info_matrix(x) / (x, theta)
  info_vector <- .normalize_info_vector(info_vector, theta)
  info_matrix <- .normalize_info_matrix(info_matrix, theta)

  bd <- .build_info_data(X, info_mode, info_vector, info_matrix, theta)
  info_raw <- bd$info_data; k <- bd$k
  if (is.null(theta)) theta <- rep(0.0, k)

  wb <- .wb_from(NULL, subset, grad_g, theta, k)
  i0 <- .make_infor0(xi0_points, xi0_weights, n0, n1,
                     info_mode, info_vector, info_matrix, theta, k)
  info_data <- .scale_info(info_raw, info_mode, i0$b)

  r_pt <- .point_info_rank(info_mode, info_vector, info_matrix, theta, X)
  min_support <- .min_support_rule(i0$infor0, k, wb, r_pt)

  prob <- list(p = as.integer(p), wb = wb, info_mode = info_mode,
               info_data = info_data, infor0 = i0$infor0, b = i0$b,
               X = X, k = k, N = ncol(X), n = nrow(X), theta = theta,
               info_vector = info_vector, info_matrix = info_matrix,
               n0 = as.numeric(n0), n1 = as.numeric(n1),
               min_support = min_support)
  class(prob) <- "DesignProblem"
  prob
}

# Information matrix at a single point x (k x k, unscaled).
info_at <- function(prob, x)
  .info_mat_at(x, prob$info_mode, prob$info_vector, prob$info_matrix, prob$theta)
