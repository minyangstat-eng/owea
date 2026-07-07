# ===========================================================================
# criterion.R -- criterion helpers.
#
# The package reports the NORMALISED Phi_p criterion everywhere:
#   pp == 0 :  log|Sigma| / v
#   pp  > 0 : ( tr(Sigma^p) / v )^(1/p) ,   Sigma = wb M^{-1} wb'.
# ===========================================================================

# Integer matrix power; p = 0 gives the identity.
matpow <- function(A, p) {
  if (p == 0) return(diag(nrow(A)))
  R <- A
  for (i in seq_len(p - 1L)) R <- R %*% A
  R
}

tr <- function(A) sum(diag(A))

# Validate the criterion index: only D (p = 0) and A (p = 1) are supported.
# Returns the coerced integer, or errors. Called by every public entry point
# that takes a criterion (p / pp).
.check_criterion <- function(p) {
  pin <- suppressWarnings(as.integer(p))
  if (length(pin) != 1L || is.na(pin) || !(pin %in% c(0L, 1L)))
    stop("Only D-optimality (p = 0) and A-optimality (p = 1) are supported; got p = ",
         paste(p, collapse = ", "), ".", call. = FALSE)
  pin
}

#' Normalised Phi_p criterion value of an approximate design.
#'
#' @param pp criterion: 0 = D-optimality or 1 = A-optimality (only these two are
#'   supported).
#' @param index 1-based support indices (columns of \code{infor_vec_all}).
#' @param weight corresponding weights.
#' @param infor_vec_all the \eqn{k \times N} information-vector matrix.
#' @param wb the \eqn{v \times k} selection matrix.
#' @return the normalised Phi_p criterion value.
#' @export
phi_value <- function(pp, index, weight, infor_vec_all, wb) {
  pp <- .check_criterion(pp)
  infor_vec_all <- as.matrix(infor_vec_all)
  k <- nrow(infor_vec_all)
  criterion_cpp(as.integer(pp), as.integer(index), as.numeric(weight),
                0L, infor_vec_all, as.matrix(wb), matrix(0.0, k, k))
}
