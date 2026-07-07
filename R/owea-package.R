#' owea: Optimal Weights Exchange Algorithm for approximate optimal designs
#'
#' Finds D- and A-optimal approximate designs for linear, nonlinear and
#' generalized linear models with the Optimal Weights Exchange Algorithm (OWEA)
#' of Yang, Biedermann & Tang (2013, JASA).
#'
#' The model enters the algorithm through a single per-point function. You may
#' supply EITHER
#' \itemize{
#'   \item an \strong{information vector} \code{info_vector(x, theta)} returning
#'     the length-\eqn{k} vector \eqn{f(x,\theta)} whose outer product
#'     \eqn{f f'} is the per-point information matrix (the fast path); or
#'   \item an \strong{information matrix} \code{info_matrix(x)} returning the
#'     full \eqn{k \times k} Fisher information matrix at \eqn{x} (the general
#'     path, for models whose per-point information is not rank one).
#' }
#' A single C++ core (RcppArmadillo) handles both representations.
#'
#' Main entry points:
#' \describe{
#'   \item{\code{\link{optimal_design}}}{continuous design box + step sequence
#'     (the recommended high-level entry; accepts either representation,
#'     existing designs, and optional merging).}
#'   \item{\code{\link{owea}} / \code{\link{DesignProblem}}}{a fixed finite
#'     candidate set.}
#'   \item{\code{\link{appro_opt}} / \code{\link{appro_opt_seq}}}{the
#'     information-vector API (fast path) kept for compatibility.}
#' }
#'
#' @references Yang, M., Biedermann, S. & Tang, E. (2013). On Optimal Designs
#'   for Nonlinear Models: A General and Efficient Algorithm. \emph{JASA}
#'   108(504), 1411-1420.
#'
#' @useDynLib owea, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
